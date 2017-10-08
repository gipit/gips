import datetime

import pytest

from ...data.daymet import daymet

# approach here is cribbed from unit/t_landsat.py

_2015002 = datetime.date(2015, 1, 2)
_afps = { # asset full paths
    'tmin': '/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_tmin.tif',
    'tmax': '/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_tmax.tif',
    'prcp': '/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_prcp.tif',
    'srad': '/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_srad.tif',
    'vp':   '/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_vp.tif',
}

# inputs and outputs for the base case - one for each asset/product type
asset_constructor_base_case_params = [
    (asset_fp, (asset_fp, ap_type, '11935', _2015002, 'daymet', {ap_type: asset_fp}))
    for ap_type, asset_fp in _afps.items()
]

@pytest.mark.parametrize("fn, expected", asset_constructor_base_case_params)
def t_daymetAsset_constructor_base_case(fn, expected):
    """Confirm well-behaved files go in correctly."""
    a = daymet.daymetAsset(fn)
    assert expected == (a.filename, a.asset, a.tile, a.date, a.sensor, a.products)

@pytest.mark.parametrize("fn, message", (
    ('/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_xxx.tif', 'bad asset/product type'),
    ('/repo/daymet/tiles/11935/2015002/11935_daymet_tmin.tif', 'missing a token'),
    ('/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_tmin.tif.index', 'trailing characters'),
    ('/repo/daymet/tiles/11935/2015002/foo_11935_2015002_daymet_tmin.tif', 'leading characters'),
    ('/repo/daymet/tiles/11935/2015002/11935_2015002_daymet_tmin_x.tif', 'extra token'),
))
def t_daymetAsset_constructor_error_case(fn, message):
    """Confirm malformed files cause exceptions."""
    with pytest.raises(ValueError, message=message):
        a = daymet.daymetAsset(fn)
