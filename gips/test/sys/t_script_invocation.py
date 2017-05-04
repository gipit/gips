import logging

import pytest
import envoy

from .util import *

logger = logging.getLogger(__name__)

pytestmark = sys # skip everything unless --sys


# TODO exit statuses of all scripts, not just inventory
def t_inventory_nonzero_exit(repo_env):
    """Confirm bad argument passing results in nonzero exit status."""
    actual = repo_env.run('gips_inventory', 'landsat', '--fetch', '--rectify')
    assert actual.exit_status != 0 
