import logging
from datetime import datetime

import pytest
import envoy

from .util import *

logger = logging.getLogger(__name__)

pytestmark = sys # skip everything unless --sys

driver = 'sentinel2'

# changing this will require changes in expected/
# date chosen for relative lack of cloud cover -------------vvv
STD_ARGS = ('sentinel2', '-s', DURHAM_SHP_PATH, '-d', '2017-010', '-v', '4')
# note that two dates doubles runtime and gips_process is approx 5 min / day runtime
#STD_ARGS = ('sentinel2', '-s', DURHAM_SHP_PATH, '-d', '2017-010,2017-030', '-v', '4')
# for now sentinel2 only supports top-of-atmo products
PROD_ARGS = list(STD_ARGS) + ['-p',
    'ref-toa', # standard
    'ndvi-toa', 'evi-toa', 'lswi-toa', 'ndsi-toa', 'satvi-toa', 'msavi2-toa', 'vari-toa', # indices
    'ndti-toa', 'crc-toa', 'crcm-toa', 'isti-toa', # tillage indices
    'bi-toa', 'sti-toa', 'brgt-toa', # these are broken atm but still undergo testing
]

@pytest.fixture
def setup_fixture(pytestconfig):
    """Use gips_inventory to ensure presence of data in the data repo."""
    if not pytestconfig.getoption('setup_repo'):
        logger.debug("Skipping repo setup per lack of option.")
        return
    logger.info("Downloading Sentinel 2 data . . .")
    cmd_str = 'gips_inventory ' + ' '.join(STD_ARGS) + ' --fetch'
    outcome = envoy.run(cmd_str)
    logger.info("Data download complete.")
    if outcome.status_code != 0:
        raise RuntimeError("Data setup via `gips_inventory` failed",
                           outcome.std_out, outcome.std_err, outcome)


def t_inventory(setup_fixture, repo_env, expected):
    """Test `gips_inventory` and confirm recorded output is given."""
    actual = repo_env.run('gips_inventory', *STD_ARGS)
    assert expected == actual


@slow
def t_process(setup_fixture, repo_env, expected):
    """Test `gips_process`, confirming products are created."""
    process_actual = repo_env.run('gips_process', *PROD_ARGS)
    inventory_actual = envoy.run('gips_inventory ' + ' '.join(STD_ARGS))
    assert (expected == process_actual
            and expected._inv_stdout == inventory_actual.std_out)


def t_info(repo_env, expected):
    """Test `gips_info` and confirm recorded output is given."""
    actual = repo_env.run('gips_info', driver)
    assert expected == actual


def t_project(setup_fixture, clean_repo_env, output_tfe, expected):
    """Test gips_project with warping."""
    args = PROD_ARGS + ['--res', '100', '100', '--outdir', OUTPUT_DIR, '--notld']
    actual = output_tfe.run('gips_project', *args)
    assert expected == actual


def t_stats(setup_fixture, clean_repo_env, output_tfe, expected):
    """Test gips_stats on projected files."""
    # generate data needed for stats computation
    args = PROD_ARGS + ['--res', '30', '30', '--outdir', OUTPUT_DIR, '--notld']
    prep_run = output_tfe.run('gips_project', *args)
    assert prep_run.exit_status == 0 # confirm it worked; not really in the test

    # compute stats
    gtfe = GipsTestFileEnv(OUTPUT_DIR, start_clear=False)
    actual = gtfe.run('gips_stats', OUTPUT_DIR)

    # check for correct stats content
    assert expected == actual

@slow
@acolite
@pytest.mark.skip(reason="Overflows in zlib.crc32")
def t_process_acolite(setup_fixture, repo_env, expected):
    """Generate acolite products."""
    aco_args = list(STD_ARGS) + \
               '-p rhow fai oc2chl oc3chl spm655 turbidity acoflags'.split()
    actual = repo_env.run('gips_process', *aco_args)
    assert expected == actual
