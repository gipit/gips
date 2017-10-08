import datetime
import fnmatch
import os

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
    """Run Asset.discover() to confirm only valid filenames are accepted."""
    # fake out config so it goes to FS inventory
    m_use_orm = mocker.patch('gips.data.core.orm.use_orm', return_value=False)
    # to confirm it wasn't called
    m_asset_search = mocker.patch('gips.data.core.dbinv.asset_search')

    # set up mocked repo path and install where needed
    repo_prefix = '/dontcare/'
    bad_filename = repo_prefix + asset_type + '.' + bad_filename_tail
    data_path = mocker.patch.object(modisAsset.Repository, 'data_path')
    data_path.return_value = repo_prefix

    # simulate glob behavior on an emulated directory listing
    fake_dir_listing = []
    for at, tail in known_good_filenames.items():
        if asset_type == at:
            fake_dir_listing.append(bad_filename)
        else:
            fake_dir_listing.append(repo_prefix + at + '.' + tail)

    isdir = mocker.patch.object(core.os.path, 'isdir')
    isdir.return_value = True
    isfile = mocker.patch.object(core.os.path, 'isfile')
    isfile.return_value = True

    listdir = mocker.patch.object(core.os, 'listdir')
    listdir.side_effect = lambda path: [os.path.basename(fn) for fn in fake_dir_listing]

    # run the test
    found = modisAsset.discover('hi mom!', datetime.date(9999, 9, 9))

    # order of found Assets depends on unpredictable dict iteration, hence set()
    actual   = set(asset.filename for asset in found)
    expected = set(l for l in fake_dir_listing if asset_type not in l)
    assert expected == actual


@pytest.fixture
def basic_asset_db(db):
    # This data isn't entirely valid but is correct enough for simple tests.  Also unicode isn't super
    # relevant here; it's an artifact of cutpasting
    path_prefix = modisAsset.Repository.data_path()
    assets = [
        { # basic row in table
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 1),
            'driver': u'modis',
            'sensor': u'MCD',
            'tile':   u'h12v04'
        },
        { # duplicate tile string to confirm nonrepetition in outcome
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 2),
            'driver': u'modis',
            'sensor': u'MCD',
            'tile':   u'h12v04'
        },
        { # second tile to confirm multiple tiles will be returned
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 2),
            'driver': u'modis',
            'sensor': u'MCD',
            'tile':   u'h13v05'
        },
        { # different driver to confirm no cross-driver pollution of results
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 1),
            'driver': u'landsat',
            'sensor': u'MCD',
            'tile':   u'trolololo'
        },
    ]
    from gips.inventory.dbinv import models
    [models.Asset(**f).save() for f in assets]
    return assets
    # TODO put this fn into utils and refactor wrt t_dbinv_api.py
    # TODO add in some disparate asset types for Asset.discover to chew on (it can return multiple
    # values if called with asset=None)


@pytest.mark.django_db
def t_discover(basic_asset_db):
    """Confirm Asset.discover works (with the inventory DB).

    Set up the inventory db with some assets then confirm
    Asset.discover's return value is correct.
    """
    disc_out = modisAsset.discover('h12v04', datetime.date(2012, 12, 1))
    actual = disc_out[0]
    expected = modisAsset(basic_asset_db[0]['name'])

    attribs = ('filename', 'asset', 'date', 'tile')
    expected_d = {a: getattr(expected, a) for a in attribs}
    actual_d   = {a: getattr(actual, a)   for a in attribs}

    assert len(disc_out) == 1 and expected_d == actual_d
