import datetime
import fnmatch

import pytest

from gips.data.modis import modisAsset
from gips.data import core

# TODO these aren't covered (because demo filenames couldn't be located in the
# central archive): MYD11A2 MOD10A2 MYD10A2 
known_good_filenames = { # cracked in half for convenience
    'MCD43A4': 'A2012336.h12v04.006.2016112010833.hdf',
    'MCD43A2': 'A2012337.h12v04.006.2016112013509.hdf',
    'MOD09Q1': 'A2000289.h08v04.005.2008205220247.hdf',
    'MOD10A1': 'A2002021.h08v04.005.2008289111219.hdf',
    'MYD10A1': 'A2002221.h08v04.005.2007173090625.hdf',
    'MOD11A1': 'A2006144.h08v04.005.2008115065841.hdf',
    'MYD11A1': 'A2006150.h08v04.005.2008116225141.hdf',
    'MOD11A2': 'A2000233.h08v04.005.2007212193347.hdf',
    'MCD12Q1': 'A2004001.h02v06.051.2014287173613.hdf',
}

bad_filename_tails = (
    'A2003266.h12v09.005.hdf',  # missing generation stamp per issue #79
    'bork.hdf',                 # not enough period-delimited substrings
    '_2004001._02_06.051.2014287173613.hdf', # try to trip it up with wrong characters
)

# generator for taking the product of asset types and bad filename tails
at_x_bft = ((at, bft) for at in known_good_filenames.keys() for bft in bad_filename_tails)

@pytest.mark.parametrize('asset_type, bad_filename_tail', at_x_bft)
def t_discover_filename_globs(mocker, asset_type, bad_filename_tail):
    repo_prefix = '/dontcare/'
    bad_filename = repo_prefix + asset_type + '.' + bad_filename_tail

    # use known repo value instead of looking in config
    data_path = mocker.patch.object(modisAsset.Repository, 'data_path')
    data_path.return_value = repo_prefix

    # simulate glob behavior on an emulated directory listing
    fake_dir_listing = []
    for at, tail in known_good_filenames.items():
        if asset_type == at:
            fake_dir_listing.append(bad_filename)
        else:
            fake_dir_listing.append(repo_prefix + at + '.' + tail)

    glob = mocker.patch.object(core.glob, 'glob')
    glob.side_effect = lambda path: [fn for fn in fake_dir_listing if fnmatch.fnmatchcase(fn, path)]

    found = modisAsset.discover('hi mom!', datetime.date(9999, 9, 9))
    # order of found Assets depends on unpredictable dict iteration, hence set()
    actual   = set(asset.filename for asset in found)
    expected = set(l for l in fake_dir_listing if asset_type not in l)
    assert expected == actual
