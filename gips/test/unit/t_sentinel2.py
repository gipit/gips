import pytest

from ...data.sentinel2 import sentinel2

# pattern borrowed from t_landsat.py

old_style_fn = (
    '/data-root/sentinel2/tiles/19TCH/2016019/'
    '19TCH_S2A_OPER_PRD_MSIL1C_PDMC_20160119T230435_R011_V20160119T155447_20160119T155447.zip')
new_style_fn = (
    '/data-root/sentinel2/tiles/19TCH/2017030/'
    'S2A_MSIL1C_20170130T153551_N0204_R111_T19TCH_20170130T153712.zip')

@pytest.mark.parametrize("fn, expected", [
    (old_style_fn, (old_style_fn, 'L1C', 'S2A', '19TCH', 2016, 19)),
    (new_style_fn, (new_style_fn, 'L1C', 'S2A', '19TCH', 2017, 30)),
])
def t_sentinel2Asset_constructor(mocker, fn, expected):
    """Test asset constructor against old & new-style assets."""
    mocker.patch.object(sentinel2.zipfile, 'ZipFile')
    asset = sentinel2.sentinel2Asset(fn)
    actual = (asset.filename, asset.asset, asset.sensor, asset.tile,
              asset.date.year, asset.date.timetuple().tm_yday)
    assert expected == actual
