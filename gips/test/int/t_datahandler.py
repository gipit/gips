
import os, shutil
import datetime

import pytest
import envoy

from gips.inventory import dbinv, orm
from gips.test.sys import util, expected

from gips.datahandler import worker

@pytest.mark.django_db()
def t_fetch():
    """Simple test for worker.fetch.  Warning:  Operates destructively."""
    from gips.datahandler.worker import fetch
    import gippy

    # config
    gippy.Options.SetVerbose(4) # substantial verbosity for testing purposes

    os.environ['GIPS_ORM'] = 'true'
    key_kwargs = dict(driver='modis',
                      asset ='MCD43A2',
                      tile  ='h12v04',
                      date  =datetime.date(2012, 12, 1))
    fetch_args = tuple(key_kwargs[k] for k in ('driver', 'asset', 'tile', 'date'))
    api_kwargs = dict(sensor='dontcare',
                      name  ='unspecified',
                      status='scheduled',
                      **key_kwargs)
    expected = os.path.join(
            util.DATA_REPO_ROOT,
            # taken from the modis fetch test; hash isn't used but is:  531008224
            'modis/tiles/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf')

    # setup
    orm.setup()
    if os.path.exists(expected): # remove the asset if it's present
        os.remove(expected)
    dbinv.update_or_add_asset(**api_kwargs) # set the asset to 'scheduled' status

    # test
    returned_asset = fetch(*fetch_args)
    queried_asset = dbinv.models.Asset.objects.get(**key_kwargs)
    assert (expected == returned_asset.name == queried_asset.name
            and 'complete' == returned_asset.status == queried_asset.status
            and os.path.isfile(expected))


@pytest.mark.django_db()
def t_process():
    """Simple test for worker.process.  Warning:  Operates destructively."""
    from gips.datahandler.worker import process
    import gippy

    # check dependency:  if it's not there the test fails (to get it run the prior fetch test)
    dependee_fname = os.path.join(
            util.DATA_REPO_ROOT,
            'modis/tiles/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf')
    assert os.path.isfile(dependee_fname), "Couldn't satisfy asset dependency, giving up"
    dependee_model = dbinv.models.Asset(driver='modis',     asset ='MCD43A2',
                                        tile  ='h12v04',    date  =datetime.date(2012, 12, 1),
                                        sensor='dontcare',  name  =dependee_fname,
                                        status='complete')
    dependee_model.save()

    # config
    gippy.Options.SetVerbose(4) # substantial verbosity for testing purposes

    os.environ['GIPS_ORM'] = 'true'
    key_kwargs = dict(driver ='modis',
                      product='quality',
                      tile   ='h12v04',
                      date   =datetime.date(2012, 12, 1))
    process_args = tuple(key_kwargs[k] for k in ('driver', 'product', 'tile', 'date'))
    api_kwargs = dict(sensor='dontcare',
                      name  ='unspecified',
                      status='scheduled',
                      **key_kwargs)
    fname = os.path.join(
            util.DATA_REPO_ROOT,
            'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_quality.tif')

    # setup
    orm.setup()
    if os.path.lexists(fname): # remove the prod if it's present
        os.remove(fname)
    dbinv.update_or_add_product(**api_kwargs) # set the prod to 'scheduled' status

    # test
    returned_prod = process(*process_args)
    queried_prod = dbinv.models.Product.objects.get(**key_kwargs)

    expected = {'fname': fname,              'status': 'complete'}
    returned = {'fname': returned_prod.name, 'status': returned_prod.status}
    queried  = {'fname': queried_prod.name,  'status': queried_prod.status}
    assert (expected == returned == queried) and os.path.lexists(fname)


def t_aggregate_helper():
    pass


from django.core.management import call_command
#from gips.test.int.t_torque import disable_db_blocker
@pytest.fixture
def disable_db_blocker(django_db_blocker):
    with django_db_blocker.unblock():
        orm.setup()
        call_command('migrate', interactive=False)
        yield

def t_export_helper(disable_db_blocker):
    """Test gips_project modis with warping."""
    import gippy
    gippy.Options.SetVerbose(4) # substantial verbosity for testing purposes

    os.environ['GIPS_ORM'] = 'true'

    # setup - check for prereqs then clear test directory
    prerequisites = [os.path.join(util.DATA_REPO_ROOT, f) for f in (
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_fsnow.tif',
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_indices.tif',
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_quality.tif',
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_snow.tif',
        'modis/tiles/h12v04/2012336/h12v04_2012336_MOD_clouds.tif',
        'modis/tiles/h12v04/2012336/h12v04_2012336_MOD-MYD_obstime.tif',
        'modis/tiles/h12v04/2012336/h12v04_2012336_MOD-MYD_temp.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_fsnow.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_indices.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_quality.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_snow.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_clouds.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD-MYD_obstime.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD-MYD_temp.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_temp8td.tif',
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_temp8tn.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_fsnow.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_indices.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_quality.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_snow.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MOD_clouds.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MOD-MYD_obstime.tif',
        'modis/tiles/h12v04/2012338/h12v04_2012338_MOD-MYD_temp.tif',
    )]
    for p in prerequisites:
        try:
            dbinv.models.Product.objects.get(name=p) # raises on failure, which is desired here
            assert os.path.lexists(p)
        except:
            print ("Prerequisite not met, try: "
                   "`py.test --sys --clear-repo -k 'modis and process' -s --ll debug`")
            raise
    if os.path.exists(util.OUTPUT_DIR): # remove the output dir if it's present
        shutil.rmtree(util.OUTPUT_DIR)

    # set up args for the call
    # some of the same arguments used by system test for gips config
    kwargs = {'alltouch': False, 'crop': False, 'interpolation': 0, 'overwrite': False,
              'res': [100.0, 100.0], 'tree': False}

    spatial_spec = {'ptile': 0, 'pcov': 0, 'tiles': None, 'where': '', 'key': '',
                    'site': util.NH_SHP_PATH}

    # run the test
    worker._export(driver='modis', spatial_spec=spatial_spec,
                   temporal_spec={'days': None, 'dates': '2012-12-01,2012-12-03'},
                   products=None, outdir=util.OUTPUT_DIR, notld=True,
                   **kwargs)
    expected = set((
        '2012336_MCD_fsnow.tif',
        '2012336_MCD_indices.tif',
        '2012336_MCD_quality.tif',
        '2012336_MCD_snow.tif',
        '2012336_MOD_clouds.tif',
        '2012336_MOD-MYD_obstime.tif',
        '2012336_MOD-MYD_temp.tif',
        '2012337_MCD_fsnow.tif',
        '2012337_MCD_indices.tif',
        '2012337_MCD_quality.tif',
        '2012337_MCD_snow.tif',
        '2012337_MOD_clouds.tif',
        '2012337_MOD-MYD_obstime.tif',
        '2012337_MOD-MYD_temp.tif',
        '2012337_MOD_temp8td.tif',
        '2012337_MOD_temp8tn.tif',
        '2012338_MCD_fsnow.tif',
        '2012338_MCD_indices.tif',
        '2012338_MCD_quality.tif',
        '2012338_MCD_snow.tif',
        '2012338_MOD_clouds.tif',
        '2012338_MOD-MYD_obstime.tif',
        '2012338_MOD-MYD_temp.tif',
    ))
    actual = set(os.listdir(os.path.join(util.OUTPUT_DIR, '0')))
    assert expected == actual


def t_aggregate():
    pass
