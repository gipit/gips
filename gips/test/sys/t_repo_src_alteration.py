import os
import logging
from datetime import datetime

import pytest

from .util import *

logger = logging.getLogger(__name__)

# skip everything unless --src-altering
pytestmark = pytest.mark.skipif(not pytest.config.getoption("src_altering"),
                                reason="--src-altering is required for this test")

# only one day for most of these tests (to save time), and note --fetch:
STD_MODIS_ARGS = ('modis', '-s', NH_SHP_PATH, '-d', '2012-12-01,2012-12-01', '-v', '4', '--fetch')

# will need to support changing the driver on a per-test basis; one idea is to have the test set
# the driver on the fixture, then have the fixture inspect it when needed
driver = 'modis'


@pytest.yield_fixture
def careful_repo_env(request, expected):
    """Carefully provides access to the data repo.

    It won't let the test run if existing files in the repo would be
    overwritten by files created during the test."""
    logger.debug("starting careful_repo_env setup")
    gtfe = GipsTestFileEnv(DATA_REPO_ROOT, start_clear=False)
    intersection = set(gtfe._find_files().keys()) & set(expected.created.keys())
    if intersection:
        msg = '{} output files found before test; test cannot proceed'
        raise RuntimeError(msg.format(len(intersection)), list(intersection))
    logger.debug("careful_repo_env yielding to test")
    yield gtfe
    gtfe.remove_created()
    driver = request.function.func_name.split('_')[1]
    gtfe.run('gips_inventory', driver, '--rectify')


is_wednesday = datetime.today().date().weekday() == 2

@pytest.mark.skipif(is_wednesday, reason="Data isn't available on Wednesdays")
def t_modis_inv_fetch(careful_repo_env, expected):
    """Test gips_inventory --fetch; actually contacts data provider."""
    actual = careful_repo_env.run('gips_inventory', *STD_MODIS_ARGS)
    assert expected.created == actual.created


@pytest.mark.timeout(20)
@pytest.mark.skipif(not is_wednesday, reason="This test is only meaningful on Wednesdays")
# TODO use careful version once expectations are known:
# def t_modis_inv_fetch_on_wednesday(careful_repo_env):
def t_modis_inv_fetch_on_wednesday(repo_env):
    """Test gips_inventory --fetch's graceful failure during planned outages.

    Actually contacts data provider, but only on Wednesdays."""
    # TODO some files may be generated; put those in expected/ when they've been identified
    expected = GipsProcResult()
    actual = repo_env.run('gips_inventory', *STD_MODIS_ARGS)
    assert expected.created == actual.created


@pytest.mark.xfail(reason="waiting on a fix for issue #86", run=False)
def t_modis_archive(careful_repo_env, repo_env, expected):
    """Test gips_archive modis using faked source/asset files."""
    files = ( # list of fake files
        'MCD43A4.A2012337.h12v04.006.2016112013509.hdf',
        'MOD11A2.A2012337.h12v04.005.2012346152330.hdf',
        'MYD10A1.A2012337.h12v04.005.2012340112013.hdf',
        'MOD11A1.A2012337.h12v04.005.2012339204007.hdf',
        'MYD11A1.A2012337.h12v04.005.2012341072847.hdf',
        'MOD10A1.A2012337.h12v04.005.2012340033542.hdf',
        'MCD43A2.A2012337.h12v04.006.2016112013509.hdf',
        'MOD10A1.A2012336.h12v04.005.2012339213007.hdf',
        'MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
        'MYD11A1.A2012336.h12v04.005.2012341040543.hdf',
        'MYD10A1.A2012336.h12v04.005.2012340031954.hdf',
        'MOD11A1.A2012336.h12v04.005.2012339180517.hdf',
        'MCD43A4.A2012336.h12v04.006.2016112010833.hdf',
        'MCD43A2.A2012338.h12v04.006.2016112020013.hdf',
        'MYD11A1.A2012338.h12v04.005.2012341075802.hdf',
        'MOD10A1.A2012338.h12v04.005.2012341091201.hdf',
        'MOD11A1.A2012338.h12v04.005.2012341041222.hdf',
        'MYD10A1.A2012338.h12v04.005.2012340142152.hdf',
        'MCD43A4.A2012338.h12v04.006.2016112020013.hdf', # this one should test --update
    )
    # put the faked assets into place in the stage, and a fake stale asset into the archive
    for f in files:
        careful_repo_env.writefile(os.path.join('modis/stage', f),
                                   content='This file is named ' + f + '!')
    stale_file_path = os.path.join(
            DATA_REPO_ROOT,
            'modis/tiles/h12v04/2012338/MCD43A4.A2012338.h12v04.005.2012120921015.hdf')
    careful_repo_env.writefile(stale_file_path, content='This STALE file is named ' + f + '!')
    careful_repo_env.run('gips_inventory', 'modis', '--rectify') # put the stale file in the DB

    # run the test
    stage_dir = os.path.join(DATA_REPO_ROOT, 'modis/stage')
    actual = careful_repo_env.run('gips_archive', 'modis', '--update', cwd=stage_dir)
    inv_actual = repo_env.run('gips_inventory', 'modis')

    assert (expected._post_archive_inv_stdout == inv_actual.stdout and
            expected.created == actual.created and
            expected.deleted == actual.deleted)


def t_sentinel2_inv_fetch(careful_repo_env, expected):
    """Test gips_inventory sentinel2 --fetch; actually contacts data provider."""
    args = ('sentinel2', '-s', DURHAM_SHP_PATH, '-d', '2017-003', '-v', '4', '--fetch')
    actual = careful_repo_env.run('gips_inventory', *args)
    assert expected.created == actual.created
