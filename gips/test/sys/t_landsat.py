import logging

import envoy
import pytest

from .util import *

logger = logging.getLogger(__name__)

# changing this will require changes in expected/
STD_ARGS = ('landsat', '-s', NH_SHP_PATH, '-d', '2015-352', '-v', '4')

product_args = tuple('-p acca bqashadow ref-toa ndvi-toa rad-toa'.split())

STD_PROD_ARGS = STD_ARGS + product_args

@pytest.fixture
def setup_landsat_data(pytestconfig):
    """Use gips_inventory to ensure presence of MODIS data in the data repo."""
    driver = 'landsat'
    if not pytestconfig.getoption('setup_repo'):
        logger.info("Skipping repo setup per lack of option.")
        return
    cmd_str = 'gips_inventory ' + ' '.join(STD_ARGS) + ' --fetch'
    logger.info("Downloading data via `{}`".format(cmd_str))
    outcome = envoy.run(cmd_str)
    logger.info("{} data download complete.".format(driver))
    if outcome.status_code != 0:
        raise RuntimeError("{} data setup via `gips_inventory` failed".format(driver),
                           outcome.std_out, outcome.std_err, outcome)


def t_info(repo_env, expected):
    """Test `gips_info modis` and confirm recorded output is given."""
    actual = repo_env.run('gips_info', 'landsat')
    assert expected == actual


def t_inventory(setup_landsat_data, repo_env, expected):
    """Test `gips_inventory` to confirm it notices the emplaced data file."""
    actual = repo_env.run('gips_inventory', *STD_ARGS)
    assert expected == actual

@slow
def t_process(setup_landsat_data, repo_env, expected):
    """Test gips_process on landsat data."""
    actual = repo_env.run('gips_process', *STD_PROD_ARGS)
    assert expected == actual


def t_project(setup_landsat_data, clean_repo_env, output_tfe, expected):
    """Test gips_project landsat with warping."""
    args = STD_PROD_ARGS + ('--res', '30', '30', '--outdir', OUTPUT_DIR, '--notld')
    actual = output_tfe.run('gips_project', *args)
    assert expected == actual


def t_project_no_warp(setup_landsat_data, clean_repo_env, output_tfe, expected):
    """Test gips_project modis without warping."""
    args = STD_PROD_ARGS + ('--outdir', OUTPUT_DIR, '--notld')
    actual = output_tfe.run('gips_project', *args)
    assert expected == actual


def t_tiles(setup_landsat_data, clean_repo_env, output_tfe, expected):
    """Test gips_tiles modis with warping."""
    args = STD_PROD_ARGS + ('--outdir', OUTPUT_DIR, '--notld')
    actual = output_tfe.run('gips_tiles', *args)
    assert expected == actual

@slow
def t_tiles_copy(setup_landsat_data, clean_repo_env, output_tfe, expected):
    """Test gips_tiles modis with copying."""
    # doesn't quite use STD_ARGS
    # -p ref-toa ndvi-toa rad-toa
    args = ('landsat', '-t', '012030', '-d', '2015-352', '-v', '4',
            '--outdir', OUTPUT_DIR, '--notld') + product_args
    actual = output_tfe.run('gips_tiles', *args)
    assert expected == actual


def t_stats(setup_landsat_data, clean_repo_env, output_tfe, expected):
    """Test gips_stats on projected files."""
    # generate data needed for stats computation
    args = STD_PROD_ARGS + ('--res', '30', '30', '--outdir', OUTPUT_DIR, '--notld')
    prep_run = output_tfe.run('gips_project', *args)
    assert prep_run.exit_status == 0 # confirm it worked; not really in the test

    # compute stats
    gtfe = GipsTestFileEnv(OUTPUT_DIR, start_clear=False)
    actual = gtfe.run('gips_stats', OUTPUT_DIR)

    # check for correct stats content
    assert expected == actual
