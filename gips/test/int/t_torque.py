import os
import datetime
import time

import pytest
import envoy

from django.core.management import call_command

import gippy
from gips.inventory import dbinv, orm
from gips.test.sys import util
from gips.datahandler.torque import submit

@pytest.fixture
def disable_db_blocker(django_db_blocker):
    with django_db_blocker.unblock():
        orm.setup()
        call_command('migrate', interactive=False)
        yield

def t_single_fetch(disable_db_blocker):
    """Simple test for torque.submit('fetch'...).  Warning:  Destructive."""

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
    if os.path.exists(expected): # remove the asset if it's present
        os.remove(expected)
    dbinv.update_or_add_asset(**api_kwargs) # set the asset to 'requested' status

    # test
    submission_outcome = submit('fetch', [fetch_args])
    #(job_exit_status, job_stdout, job_stderr) = job_outcomes[0]
    print "qsub called; exit status, stdout, stderr: {}, {}, {}".format(*submission_outcome[0])
    # wait on outcome - once per second check status
    queried_asset = dbinv.models.Asset.objects.get(**key_kwargs)
    for i in range(1, 46):
        queried_asset.refresh_from_db()
        if queried_asset.status in ('complete', 'failed'):
            break
        time.sleep(1)
        if i % 5 == 0:
            print "Wait time: {}s".format(i)

    assert (True # TODO check outcome of job in torque via pbs_python call
            and (expected, 'complete') == (queried_asset.name, queried_asset.status)
            and os.path.isfile(expected))


def t_submit_three_fetch(disable_db_blocker):
    """Simple test for torque.submit('fetch'...).  Warning:  Destructive."""

    # config
    gippy.Options.SetVerbose(4) # substantial verbosity for testing purposes
    os.environ['GIPS_ORM'] = 'true'

    expected = []
    for file_path in (
            # taken from the modis fetch test; hash isn't used but is:  531008224
            'modis/tiles/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'modis/tiles/h12v04/2012337/MCD43A2.A2012337.h12v04.006.2016112013509.hdf',
            'modis/tiles/h12v04/2012338/MCD43A2.A2012338.h12v04.006.2016112020013.hdf'):
        expected.append(os.path.join(util.DATA_REPO_ROOT, file_path))

    # assemble args & such
    fetch_args = []
    key_kwargs = []
    api_kwargs = []
    for day in (1, 2, 3):
        today = datetime.date(2012, 12, day)
        todays_key = dict(driver='modis', asset='MCD43A2', tile='h12v04', date=today)
        fetch_args.append(['modis', 'MCD43A2', 'h12v04', today])
        key_kwargs.append(todays_key)
        api_kwargs.append(dict(sensor='dontcare', name='unspecified', status='scheduled',
                               **todays_key))

    # setup
    [os.remove(e) for e in expected if os.path.exists(e)]   # remove the asset if it's present
    [dbinv.update_or_add_asset(**ak) for ak in api_kwargs]  # set the asset to 'scheduled' status

    # test
    outcomes = submit('fetch', fetch_args, 2) # should be two batches: [two fetches, one fetch]
    for o in outcomes:
        print "qsub called; exit status, stdout, stderr: {}, {}, {}".format(*o)
    # wait on outcome - once per second check status
    queried_assets = [dbinv.models.Asset.objects.get(**kk) for kk in key_kwargs]
    for i in range(1, 61):
        [qa.refresh_from_db() for qa in queried_assets]
        if all([qa.status in ('complete', 'failed') for qa in queried_assets]):
            break # done waiting when all three are in a final state
        time.sleep(1)
        if i % 5 == 0:
            print "Wait time: {}s".format(i)

    assert {
        'status':       ['complete', 'complete', 'complete'],
        'filenames':    set(expected),
        'files_exist':  [True, True, True]
    } == {
        'status':       [qa.status for qa in queried_assets],
        'filenames':    set(qa.name for qa in queried_assets),
        'files_exist':  [os.path.isfile(e) for e in expected],
    }
