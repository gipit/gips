from datetime import date

import pytest

from ...data import chirps

@pytest.mark.parametrize("asset_fp, expected_date", [
    ('/whatever/africa-daily-chirps-v2.0.1995.05.03.tif.gz', date(1995, 5, 3)),
    ('/whatever/africa-daily-chirps-v2.0.1996.05.18.tif.gz', date(1996, 5, 18)),
    ('/whatever/africa-daily-chirps-v2.0.1997.05.01.tif.gz', date(1997, 5, 1)),
])
def t_chirpsAsset_constructor(asset_fp, expected_date):
    """Test the asset constructor using three representative assets."""
    a = chirps.chirpsAsset(asset_fp)
    expected = (asset_fp,   expected_date, 'africa', 'chirps', 'africa-daily')
    actual   = (a.filename, a.date,        a.tile,   a.sensor, a.asset)
    assert expected == actual

@pytest.mark.parametrize("asset_fp", [
    '/whatever/',
    '/whatever/whatever-daily-chirps-v2.0.1996.05.18.tif.gz',
    '/whatever/africa-daily-chirps-v2.0.1997.055.01.tif.gz',
    '/whatever/africa-daily-chirps-v2.0.1997.05.01.tif.gz/hmm',
])
def t_chirpsAsset_constructor_error(asset_fp):
    """Test asset constructor error cases, mostly the regex that parses the filename."""
    with pytest.raises(RuntimeError):
        a = chirps.chirpsAsset(asset_fp)

@pytest.fixture
def nlst_april_1992_mock(mocker):
    """Mock ftp connection to return results from a month in 1992."""
    m_ftp_connect = mocker.patch.object(chirps.chirpsAsset, 'ftp_connect')
    m_ftp_connect.return_value.nlst.return_value = [
        'chirps-v2.0.1992.04.{:0>2}.tif.gz'.format(i) for i in range(1, 31)]
    return m_ftp_connect

@pytest.mark.parametrize('asset, tile, date, expected_fn', [
    ('africa-daily', 'africa', date(1992, 4, 1),  'chirps-v2.0.1992.04.01.tif.gz'),
    ('africa-daily', 'africa', date(1992, 4, 13), 'chirps-v2.0.1992.04.13.tif.gz'),
    ('africa-daily', 'africa', date(1992, 4, 29), 'chirps-v2.0.1992.04.29.tif.gz'),
])
def t_chirpsAsset_query_provider(nlst_april_1992_mock, asset, tile, date, expected_fn):
    """Test base case of chirps' query provider method."""
    actual = chirps.chirpsAsset.query_provider(asset, tile, date)
    assert (expected_fn, None) == actual

def t_chirpsAsset_query_provider_none_found(nlst_april_1992_mock):
    """Test no-asset-found case for chirps' query provider."""
    #with pytest.raises(ValueError):
    actual = chirps.chirpsAsset.query_provider('africa-daily', 'africa', date(1992, 5, 1))
    assert (None, None) == actual

def t_chirpsAsset_query_provider_too_many_found(nlst_april_1992_mock):
    """Test no-asset-found case for chirps' query provider."""
    nlst_april_1992_mock.return_value.nlst.return_value.append('chirps-v2.0.1992.04.01.tif.gz')
    with pytest.raises(ValueError):
        chirps.chirpsAsset.query_provider('africa-daily', 'africa', date(1992, 4, 1))
