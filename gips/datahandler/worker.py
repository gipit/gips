"""Worker functions for GIPS scheduled tasks.

These functions are called inside worker processes.
"""
import os, shutil

from django.db import transaction, IntegrityError

from gips import utils
from gips.core import SpatialExtent, TemporalExtent
from gips.datahandler import api
from . import queue
from gips.datahandler.logger import Logger
from gips.inventory import DataInventory, ProjectInventory
from gips.inventory import dbinv, orm
from gips.scripts.spatial_wrapper import aggregate

def _de_u(u):
    return u.encode('ascii', 'ignore')


def query (job_id):
    """Determine which assets to fetch and products to process"""
    with transaction.atomic():
        job = dbinv.models.Job.objects.get(pk=job_id)
        if job.status != 'initializing':
            err_msg = "Expected Job status to be 'initializing', but got {}"
            utils.verbose_out(err_msg.format(job.status), 1)
            return job # not sure this is useful for anything
        job.status = 'scheduled'
        job.save()

    api.query_service(
        job.variable.driver,
        eval(job.spatial),
        eval(job.temporal),
        [job.variable.product])

    job.refresh_from_db()
    if job.status != 'scheduled':
        # sanity check; have to keep going, but whine about it
        err_msg = "Expected Job status to be 'initializing, but got {}"
        utils.verbose_out(err_msg.format(job.status), 1)
    job.status = 'in-progress'
    job.save()

    return job


def fetch(asset_id):
    """Fetch the asset file specified."""
    # Notify everyone that this process is doing the fetch.
    Logger().log("begin")
    with transaction.atomic():
        # TODO confirm transaction.atomic prevents DB writes while it's active
        asset = dbinv.models.Asset.objects.get(pk=asset_id)
        if asset.status != 'scheduled':
            # TODO log/msg about giving up here
            return asset
        # TODO: this status change needs to be made earlier or could be scheduled twice
        asset.status = 'in-progress'
        asset.save()

    Logger().log("fetching {} {} {} {}".format(asset.driver, asset.asset, asset.tile, asset.date))
    # locate the correct fetch function & execute it, then archive the file and update the DB
    DataClass  = utils.import_data_class(asset.driver)
    AssetClass = DataClass.Asset
    filenames  = AssetClass.fetch(asset.asset, asset.tile, asset.date)
    if len(filenames) == 0:
        asset.status = 'retry'
        asset.save()
        return asset  # TODO: seems odd that this function returns
                      #       asset...discuss with tolson and fisk.

    a_obj = AssetClass._archivefile(filenames[0])[0] # for now neglect the 'update' case

    Logger().log('fetched')
    # update DB now that the work is done; no need for an atomic transaction
    # because the only critical action is Model.save(), which is already atomic
    asset.refresh_from_db()
    if asset.status != 'in-progress':
        # sanity check; have to keep going but do whine about it
        err_msg = "Expected Asset status to be 'in-progress' but got '{}'"
        utils.verbose_out(err_msg.format(asset.status), 1)
    asset.sensor = a_obj.sensor
    asset.name   = a_obj.archived_filename
    asset.status = 'complete'
    asset.save()

    Logger().log('done')
    if len(filenames) > 1:
        err_msg = 'Expected to fetch one asset but fetched {} instead'.format(len(filenames))
        raise ValueError(err_msg, filenames)
    return asset


def process(product_id):
    """Produce the specified product file."""
    # Notify everyone that this process is doing the fetch.
    with transaction.atomic():
        product = dbinv.models.Product.objects.get(pk=product_id)
        if product.status != 'scheduled':
            return product
        product.status = 'in-progress'
        product.save()

    DataClass = utils.import_data_class(_de_u(product.driver))
    data = DataClass(_de_u(product.tile), product.date) # should autopopulate with the right assets
    data.process([_de_u(product.product)])

    # sanity check
    product.refresh_from_db()
    if product.status != 'complete': # check process()'s work
        # sanity check; have to keep going but do whine about it
        err_msg = "Expected product status to be 'complete' but got '{}'"
        utils.verbose_out(err_msg.format(product.status), 1)

    return product


def _export(driver, spatial_spec, temporal_spec, products, outdir,
            notld=True, suffix='', nproc=1, start_ext=None, end_ext=None,
            **kwargs):
    """Performs a 'gips_project' with the given paramenters."""

    DataClass = utils.import_data_class(driver)

    # TODO trash and recreate output dir

    extents = SpatialExtent.factory(DataClass, **spatial_spec)
    if start_ext is not None and end_ext is not None:
        extents = extents[start_ext:end_ext]

    # create tld: SITENAME--KEY_DATATYPE_SUFFIX
    if notld:
        tld = outdir
    else:
        key = spatial_spec.get('key', '')
        if key:
            key = '--' + key
        suffix = '_' + suffix if suffix else ''
        res = kwargs.get('res', '')
        if res:
            res = '_{}x{}'.format(res[0], res[1])
        bname = (
            extents[0].site.LayerName() +
            key + res + '_' + driver + suffix
        )
        tld = os.path.join(outdir, bname)

    t_extent = TemporalExtent(**temporal_spec)

    def _mosaic(extent):
        inv = DataInventory(DataClass, extent, t_extent, products, **kwargs)
        datadir = os.path.join(tld, extent.site.Value())
        if inv.numfiles > 0:
            inv.mosaic(datadir=datadir, **kwargs)
        else:
            utils.verbose_out(
                'No data found for {} within temporal extent {}'
                .format(str(t_extent), str(t_extent)),
                2,
            )

    if nproc == 1:
        map(_mosaic, extents)
    else:
        pool = Pool(processes=nproc)
        pool.map(_mosaic, extents)
        
    # TODO nothing meaningful to return?


def _aggregate(job, outdir, nproc=1):
    aggregate(job, outdir,  nproc)


def export_and_aggregate(job_id, start_ext, end_ext,
                         nprocs=1, outdir=None, cleanup=False, **mosaic_kwargs):
    """Does the same things as gips_project + zonal summary."""
    with transaction.atomic():
        job = dbinv.models.Job.objects.get(
            pk=job_id,
        )
        task = job.postprocessjobs_set.get(
            args=repr((job_id, start_ext, end_ext))
        )
        if task.status != 'scheduled':
            # TODO log/msg about giving up here
            return task # not sure this is useful for anything
        task.status = 'in-progress'
        task.save()

    # setup output dir
    if outdir is None:
        outdir = os.path.join(utils.settings().EXPORT_DIR,
                              queue.get_job_name(),
                              str(job_id),)
    # poor man's binary semaphore since mkdir is atomic;
    # exception on a priori existence is what we want here
    os.makedirs(outdir)

    # force alltouch to be true for the geokit API
    # TODO: GIPS API should expose spatial options
    mosaic_kwargs['alltouch'] = True
    _export(
        task.job.variable.driver,
        eval(task.job.spatial),
        eval(task.job.temporal),
        [_de_u(task.job.variable.product)],
        outdir,
        start_ext=start_ext,
        end_ext=end_ext,
        **mosaic_kwargs
    )
    _aggregate(job, outdir, nprocs)

    # bookkeeping & cleanup
    task.refresh_from_db()
    if task.status != 'in-progress':
        # sanity check; have to keep going, but whine about it
        err_msg = "Expected Job status to be 'post-processing', but got {}"
        utils.verbose_out(err_msg.format(job.status), 1)
    task.status = 'complete'
    task.save()

    if cleanup:
        shutil.rmtree(outdir)
