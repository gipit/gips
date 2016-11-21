"""Unit tests for core functions, such as those found in gips.core and gips.data.core."""

import sys
from datetime import datetime

import pytest
import mock

import gips
from gips import core
from gips.data.landsat.landsat import landsatRepository, landsatData
from gips.inventory import dbinv

def t_version_override(mocker):
    """Test gips.__init__.detect_version() for correct override of __version__."""
    env = mocker.patch.object(gips.os, 'environ')
    # os.environ.get is called by libs as well as detect_version(); fortunately no harm seems to
    # come from giving them bad results.

    # no override requested
    env.get.side_effect = lambda key, default=None: default # key not found
    version_a = gips.detect_version()

    # override requested
    env.get.side_effect = lambda key, default=None: 'fancy-new-version'
    version_b = gips.detect_version()

    env.get.assert_has_calls([ # assert two identical calls
        mock.call('GIPS_OVERRIDE_VERSION', gips.version.__version__) for _ in range(2)
    ])
    assert (version_a, version_b) == (gips.version.__version__, 'fancy-new-version')


def t_repository_find_tiles_normal_case(mocker):
    """Test Repository.find_tiles using landsatRepository as a guinea pig."""
    m_list_tiles = mocker.patch('gips.data.core.dbinv.list_tiles')
    expected = [u'tile1', u'tile2', u'tile3'] # list_tiles returns unicode
    m_list_tiles.return_value = expected
    actual = landsatRepository.find_tiles()
    assert expected == actual


def t_repository_find_tiles_error_case(mocker):
    """Confirm Repository.find_tiles quits on error."""
    m_list_tiles = mocker.patch('gips.data.core.dbinv.list_tiles')
    m_list_tiles.side_effect = RuntimeError('AAAAAAAAAAH!') # intentionally break list_tiles

    with pytest.raises(RuntimeError):
        landsatRepository.find_tiles()


def t_repository_find_dates_normal_case(mocker):
    """Test Repository.find_dates using landsatRepository as a guinea pig."""
    m_list_dates = mocker.patch('gips.data.core.dbinv.list_dates')
    dt = datetime
    expected = [dt(1900, 1, 1), dt(1950, 10, 10), dt(2000, 12, 12)]
    m_list_dates.return_value = expected
    actual = landsatRepository.find_dates('some-tile')
    assert expected == actual


@pytest.mark.skip(reason="Letting exception bubble up for now; if that changes un-skip this test.")
def t_repository_find_dates_error_case(mocker):
    """Confirm Repository.find_dates quits on error."""
    m_list_dates = mocker.patch('gips.data.core.dbinv.list_dates')
    m_list_dates.side_effect = Exception('AAAAAAAAAAH!') # intentionally break list_dates

    # confirm call was still a success via the righ code path
    with pytest.raises(SystemExit):
        landsatRepository.find_dates('some-tile')


@pytest.mark.parametrize('add_to_db', (False, True))
def t_data_add_file(mocker, add_to_db):
    """Basic test for Data.AddFile; calls it once then tests its state."""
    mocker.patch('gips.data.core.orm.use_orm', return_value=True)
    m_uoap = mocker.patch('gips.data.core.dbinv.update_or_add_product')
    t_sensor    = 'test-sensor'
    t_product   = 'test-product'
    t_filename  = 'test-filename.tif'
    lsd = landsatData('test-tile', datetime(1955, 11, 5), search=False)

    lsd.AddFile(t_sensor, t_product, t_filename, add_to_db)

    # alas, cleanest to make multiple assertions for this test
    if add_to_db:
        m_uoap.assert_called_once_with(driver='landsat', product=t_product, sensor=t_sensor,
                                       tile=lsd.id, date=lsd.date, name=t_filename)
    else:
        m_uoap.assert_not_called()
    assert (lsd.filenames == {(t_sensor, t_product): t_filename}
            and lsd.sensors == {t_product: t_sensor})


@pytest.mark.django_db
def t_data_add_file_repeat(mocker):
    """Confirm that calling Data.AddFile twice results in overwrite.

    Thus, confirm it's possible to replace file entries with new versions.
    """
    mocker.patch('gips.data.core.orm.use_orm', return_value=True)
    t_tile      = 'test-tile'
    t_date      = datetime(1955, 11, 5)
    t_sensor    = 'test-sensor'
    t_product   = 'test-product'
    t_filename  = 'test-filename.tif'
    t_new_filename = 'test-new-filename.tif'
    lsd = landsatData(search=False) # heh
    lsd.id = t_tile
    lsd.date = t_date
    lsd.AddFile(t_sensor, t_product, t_filename)
    lsd.AddFile(t_sensor, t_product, t_new_filename) # should overwrite the old one
    # confirm overwrite happens both in the Data object and in the database
    assert (t_new_filename == lsd.filenames[(t_sensor, t_product)] and
            t_new_filename == dbinv.models.Product.objects.get(tile=t_tile, date=t_date).name)


@pytest.mark.parametrize('search', (False, True))
def t_data_init_search(mocker, search):
    """Confirm Data.__init__ searches the FS only when told to.

    Do this by instantiating landsatData."""
    # set up mocks:  prevent it from doing I/O and use for assertions below
    mocker.patch.object(landsatData.Asset, 'discover')
    mocker.patch.object(landsatData, 'ParseAndAddFiles')

    lsd = landsatData(tile='t', date=datetime(9999, 1, 1), search=search) # call-under-test

    # assert normal activity & entry of search block
    assert (lsd.id == 't'
            and lsd.date == datetime(9999, 1, 1)
            and '' not in (lsd.path, lsd.basename)
            and lsd.assets == lsd.filenames == lsd.sensors == {}
            and lsd.ParseAndAddFiles.called == lsd.Asset.discover.called == search)


@pytest.fixture
def df_mocks(mocker):
    """Mocks for testing Data.fetch below."""
    mocker.patch.object(landsatData.Asset, 'discover', return_value=False)
    return mocker.patch.object(landsatData.Asset, 'fetch')

# useful constant for the following tests
df_args = (['rad', 'ndvi', 'bqashadow'], ['012030'], core.TemporalExtent('2012-12-01'))

def t_data_fetch_base_case(df_mocks):
    """Test base case of data.core.Data.fetch.

    It should return data about the fetch on success, and shouldn't
    raise an exception."""
    assert landsatData.fetch(*df_args) == [('DN', '012030', datetime(2012, 12, 1, 0, 0))]


def t_data_fetch_error_case(df_mocks):
    """Test error case of data.core.Data.fetch.

    It should return [], and shouldn't raise an exception."""
    df_mocks.side_effect = RuntimeError('aaah!')
    assert landsatData.fetch(*df_args) == []
