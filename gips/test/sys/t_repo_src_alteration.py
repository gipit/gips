import logging
from datetime import datetime

import pytest

from .util import *

logger = logging.getLogger(__name__)


@pytest.yield_fixture
def careful_repo_env(request, expected):
    """Carefully provides access to the data repo.

    It won't let the test run if existing files in the repo would be
    overwritten by files created during the test."""
    gtfe = GipsTestFileEnv(DATA_REPO_ROOT, start_clear=False)
    intersection = set(gtfe._find_files().keys()) & set(expected.created.keys())
    if intersection:
        msg = '{} output files found before test; test cannot proceed'
        raise RuntimeError(msg.format(len(intersection)), list(intersection))
    yield gtfe
    gtfe.remove_created()


src_altering = pytest.mark.skipif(not pytest.config.getoption("src_altering"),
                                  reason="--src-altering is required for this test")

is_wednesday = datetime.today().date().weekday() == 2

@src_altering
@pytest.mark.skipif(is_wednesday, reason="Data isn't available on Wednesdays")
def t_modis_inv_fetch(careful_repo_env, expected):
    """Test gips_inventory --fetch; actually contacts data provider."""
    # only get data for one day to save time
    args = ('modis', '-s', NH_SHP_PATH,
            '-d', '2012-12-01,2012-12-01', '-v', '4', '--fetch')
    actual = careful_repo_env.run('gips_inventory', *args)
    assert expected == actual


@pytest.mark.timeout(20)
@pytest.mark.skipif(not is_wednesday, reason="This test is only meaningful on Wednesdays")
# TODO use careful version once expectations are known:
# def t_modis_inv_fetch_on_wednesday(careful_repo_env):
def t_modis_inv_fetch_on_wednesday(repo_env):
    """Test gips_inventory --fetch's graceful failure during planned outages.

    Actually contacts data provider, but only on Wednesdays."""
    # only get data for one day to save time
    args = ('modis', '-s', NH_SHP_PATH,
            '-d', '2012-12-01,2012-12-01', '-v', '4', '--fetch')
    # TODO some files may be generated; put those in expected/ when they've been identified
    expected = GipsProcResult()
    actual = repo_env.run('gips_inventory', *args)
    assert expected == actual
