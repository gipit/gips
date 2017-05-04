#!/usr/bin/env python

"""This is "the cron job" that facilitates work moving through the system.

It's intended to run periodically.  When it does, it looks for work that
is ready to be performed, and submits that work to work queues.  Each
'schedule_something' function handles one part of the work of a Job,
which always has four parts, which must run in sequence:  query, fetch,
process, and post-process.  The post-process step is also known as
export-and-aggregate.
"""

from django.core.exceptions import ObjectDoesNotExist
from django.db import transaction
from django.db import models

import fcntl
import math
from pprint import pprint, pformat

from gips import utils
from gips.core import SpatialExtent, TemporalExtent
from gips.inventory import dbinv, orm
from gips.datahandler import api
from gips.datahandler.logger import Logger
from gips.datahandler import queue


class Lock:
    def __init__(self, filename):
        self.filename = filename
        # This will create it if it does not exist already
        self.handle = open(filename, 'w')

    # Bitwise OR fcntl.LOCK_NB if you need a non-blocking lock
    def acquire(self):
        fcntl.flock(self.handle, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def release(self):
        fcntl.flock(self.handle, fcntl.LOCK_UN)

    def __del__(self):
        self.handle.close()


def schedule_query ():
    """Find new (ie unworked) Jobs and submit them to a work queue.

    The first step is determining assets that need to be fetched and
    products that need to be processed.  Will be called by a cron'd
    scheduler job.
    """
    orm.setup()
    jfilter = dbinv.models.Job.objects.filter
    # locate work-ready Jobs and reserve them
    with transaction.atomic():
        jobs = jfilter(status='requested')
        if not jobs.exists():
            return
        job_pks = [j.pk for j in jobs] # this has to go first, else QuerySet weirdness
        jobs.update(status='initializing')

    # submit the work to the queue, but unreserve the Jobs if anything breaks
    try:
        call_signatures = [[pk] for pk in job_pks]
        outcomes = queue.submit('query', call_signatures, 1)
    except Exception as e:
        Logger().log("Error submitting query for jobs, IDs {}:\n{}".format(job_pks, e))
        # possibly better to go to 'failed' instead?
        jfilter(pk__in=job_pks, status='initializing').update(status='requested')
        raise
    ids = [{'torque_id': o[0], 'db_id': o[1][0][0]} for o in outcomes]
    Logger().log("Submitted query for job(s):\n{}".format(pformat(ids)))


def schedule_fetch (driver):
    """Find unworked fetch tasks and submit them to a work queue."""
    orm.setup()

    with transaction.atomic():
        # check if another fetch is currently running
        active = dbinv.models.Asset.objects.filter(
            driver=driver,
            status__in=['scheduled', 'in-progress', 'retry']
        )

        # check if jobs still alive
        active_jobs = active.order_by('sched_id').distinct('sched_id')
        for j in active_jobs:
            if j.sched_id and queue.is_job_alive(j.sched_id):
                Logger().log(
                    'Previous {} fetch job still running. Moving on...'
                    .format(driver)
                )
                return

        # if anything leftover, clean it up
        needs_cleanup = active.filter(
            assetstatuschange__status='requested',
        ).annotate(
            num_retries=models.Count('assetstatuschange__status'),
        )
        if needs_cleanup.exists():
            # something went wrong. not 'complete', but no running/scheduled job
            Logger().log('Attempting to reschedule incomplete tasks')

            retries = []
            fails = []
            for a in needs_cleanup:
                if a.num_retries >= 3:
                    fails.append(a)
                    # if it has been retried 3 times, give up
                    a.status = 'failed'
                else:
                    retries.append(a)
                    a.status = 'requested'
                    a.sched_id = ''
                a.save()

            if fails:
                Logger().log(
                    "Fetch failed three times - give up\n{}"
                    .format(pformat([(f.id, f.driver, f.asset, f.tile, f.date)
                                     for f in fails]))
                )
            if retries:
                Logger().log(
                    "Retrying fetch\n{}"
                    .format(pformat([(r.id, r.driver, r.asset, r.tile, r.date)
                                     for r in retries]))
                )

    # TODO: these should be configured elsewhere
    num_jobs = 10
    per_job = 10
    asfu = dbinv.models.Asset.objects.select_for_update
    # locate fetch-ready Assets and reserve them
    with transaction.atomic():
        assets = asfu().filter(
            driver=driver,
            status='requested',
        ).order_by('id')[:num_jobs*per_job]
        if not assets.exists():
            return
        asset_pks = [a.pk for a in assets]
        dbinv.models.update_status(assets, 'scheduled')

    # submit the work to the queue, but unreserve the Assets if anything breaks
    try:
        call_signatures = [[pk] for pk in asset_pks]
        outcomes = queue.submit('fetch', call_signatures, per_job, chain=True)
        print outcomes
    except Exception as e:
        Logger().log("Error submitting asset for jobs, IDs {}:\n{}".format(asset_pks, e))
        # possibly better to go to 'failed' instead?
        with transaction.atomic():
            assets = asfu().filter(pk__in=asset_pks, status='scheduled')
            dbinv.models.update_status(assets, 'requested')
        raise

    # finally, bookkeeping:  Remember which fetch job is working on which asset
    # glean the mapping of job IDs to asset IDs from the submit() return value
    ids = [{'job_id': o[0], 'assets': [a[0] for a in o[1]]} for o in outcomes]
    with transaction.atomic():
        for i in ids:
            a = asfu().filter(id__in=i['assets']).update(sched_id=i['job_id'])
            print "len", a
    Logger().log('scheduled fetch job(s):\n{}'.format(pformat(ids)))


def schedule_process ():
    """Find ready processing tasks and submit them to a work queue.

    To be ready means the task hasn't been started already, and all its
    prerequisite assets are in place.
    """
    orm.setup()
    psfu = dbinv.models.Product.objects.select_for_update
    with transaction.atomic():
        products = psfu().filter(
            status='requested'
        ).exclude(
            ~models.Q(assetdependency__asset__status='complete')
        )
        if not products.exists():
            return
        dbinv.models.update_status(products, 'scheduled')
        product_pks = [p.pk for p in products]

    # submit the work to the queue, but unreserve the Products if anything breaks
    try:
        call_signatures = [[pk] for pk in product_pks]
        outcomes = queue.submit('process', call_signatures, 5)
        print outcomes
    except Exception as e:
        Logger().log("Error submitting products for processing, IDs {}:\n{}".format(product_pks, e))
        # possibly better to go to 'failed' instead?
        with transaction.atomic():
            products = psfu().filter(pk__in=product_pks, status='scheduled')
            dbinv.models.update_status(products, 'requested')
        raise

    # finally, bookkeeping:  Remember which job is working on which product
    # glean the mapping of job IDs to product IDs from the submit() return value
    ids = [{'torque_id': o[0], 'products': [p[0] for p in o[1]]} for o in outcomes]
    with transaction.atomic():
        for i in ids:
            products = dbinv.models.Product.objects.filter(id__in=i['products'])
            products.update(sched_id=i['torque_id'])
    Logger().log('scheduled process job(s):\n{}'.format(pformat(ids)))


def schedule_export_and_aggregate ():
    """Find ready post-processing tasks and submit them to a work queue.

    To be ready means that the task is unworked and its prerequisites
    are in place, in this case the products created by the processing
    step.
    """
    orm.setup()
    # check if any in-progress jobs have finished fetch/process
    jobs = dbinv.models.Job.objects.filter(
        status='in-progress',
    )
    for job in jobs:
        status = api.processing_status(
            job.variable.driver.encode('ascii', 'ignore'),
            eval(job.spatial),
            eval(job.temporal),
            [job.variable.product.encode('ascii', 'ignore')],
        )
        if (
                status['requested'] == 0
                and status['scheduled'] == 0
                and status['in-progress'] == 0
                and status['retry'] == 0
        ):
            Logger().log('job {} has finished fetch and process'.format(job.pk))
            # nothing being worked on. must be done pre-processing
            DataClass = utils.import_data_class(job.variable.driver.encode('ascii', 'ignore'))
            spatial_spec = eval(job.spatial)
            extents = SpatialExtent.factory(DataClass, **spatial_spec)

            num_ext = len(extents)
            print "num_ext", num_ext
            # TODO: make this configurable
            ext_per = 10 if num_ext > 15 else num_ext
            make_args = lambda i: [
                job.pk,
                i*ext_per,
                (i+1)*ext_per if (i+1)*ext_per < num_ext else num_ext
            ]
            job_args = [make_args(i)
                        for i in range(int(math.ceil(num_ext / ext_per)))]
            # create the PostProcessJobs rows in the DB before calling the remote worker
            with transaction.atomic():
                for j in job_args:
                    pp = dbinv.models.PostProcessJobs.objects.update_or_create(
                        defaults = {'status': 'scheduled'},
                        job = job,
                        args = repr(tuple(j)),
                    )
            # call the remote worker to work the PPJ(s?) then update the Job status
            with transaction.atomic():
                outcomes = queue.submit('export_and_aggregate', job_args, 1)
                for o in outcomes:
                    pp = dbinv.models.PostProcessJobs.objects.select_for_update().get(
                        job = job,
                        args = repr(tuple(o[1][0])),
                    )
                    pp.sched_id = o[0]
                    pp.save()
                job.refresh_from_db()
                job.status = 'post-processing'
                job.save()

                Logger().log('scheduled post-process job(s):\n{}'.format(pformat(outcomes)))

    # check whether any post-processing jobs have finished
    with transaction.atomic():
        jobs = dbinv.models.Job.objects.select_for_update().filter(
            status='post-processing',
        )
        for job in jobs:
            tasks = job.postprocessjobs_set.all()
            failed = False
            complete = True
            for task in tasks:
                if task.status != 'complete':
                    complete = False
                if task.status in ('scheduled', 'in-progress'):
                    # make sure jobs are still queued / running
                    if not queue.is_job_alive(task.sched_id):
                        Logger().log("task {} for job {} failed to complete"
                                     .format(task.sched_id, job.pk))
                        task.status = 'failed'
                        task.save()
                        failed = True
                    else:
                        # something is still working. move on
                        break

            if failed:
                job.status = 'failed'
                job.save()
            elif complete:
                job.status = 'complete'
                job.save()


def print_outcomes(kind, outcomes):
    if outcomes is not None and len(outcomes) > 0: # then there was something to do:
        print kind, "job(s) located and submitted:"
        [pprint(outcome) for outcome in outcomes]


def main ():
    lock = Lock('/tmp/geokitd.lock')
    try:
        lock.acquire()
    except IOError:
        Logger().log("previous scheduler process still running")
        exit(1)
    
    outcomes = schedule_query()
    print_outcomes('Query', outcomes)

    for driver in utils.settings().REPOS.keys():
        outcomes = schedule_fetch(driver)
        print_outcomes('Fetch', outcomes)

    outcomes = schedule_process()
    print_outcomes('Process', outcomes)

    outcomes = schedule_export_and_aggregate()
    print_outcomes('Export & Aggregate', outcomes)

    lock.release()

if __name__ == '__main__':
    main()
