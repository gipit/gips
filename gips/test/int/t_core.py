"""Integration tests for core functions, such as those found in gips.core and gips.data.core."""

import os

from gips.data.landsat import landsat

def t_make_temp_proc_dir():
    """Test Data.make_temp_proc_dir:

    Confirm creation of tempdir in stage/ with the right name, and
    destruction of same on context exit.
    """
    lsd = landsat.landsatData(search=False)
    with lsd.make_temp_proc_dir() as tempdir:
        tempdir_existed = os.path.exists(tempdir)

    # confirm was here, now it's gone, and had the right name
    dir_prefix = lsd.Repository.path('stage') + '/proc'
    assert tempdir_existed and not os.path.exists(tempdir) and tempdir.startswith(dir_prefix)
