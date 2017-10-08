import datetime
import os

import pytest

from ...data.cdl import cdl

_2016 = datetime.date(2016, 1, 1)
_afps = {
    'cdl': '/repo/cdl/tiles/IN/2016/IN_2016_cdl_cdl.tif',
    'cdlmkii': '/repo/cdl/tiles/IN/2016/IN_2016_cdl_cdlmkii.zip',
}

def _patched_cdlData_init(self, *args, **kwargs):
    asset = cdl.cdlAsset(os.path.dirname(__file__) + '/../IN_2016_cdl_cdlmkii.zip')
    self.assets = {'cdlmkii': asset}

def _patched_metadata(self):
    return {
        'CLASS_NAME_class1': 0,
        'CLASS_NAME_class2': 1,
    }

def _patched_gdal_open(self, *args, **kwargs):
    return type('osgeo.gdal.Dataset', (object,), {'GetMetadata': _patched_metadata})()

asset_constructor_base_case_params = [
    (asset_fp, (asset_fp, ap_type, 'IN', _2016, 'cdl', {ap_type: asset_fp}))
    for ap_type, asset_fp in _afps.items()
]

@pytest.mark.parametrize("fn, expected", asset_constructor_base_case_params)
def t_cdlAsset_constructor(fn, expected):
    a = cdl.cdlAsset(fn)
    assert expected == (a.filename, a.asset, a.tile, a.date, a.sensor, a.products)

def t_cdlData_legend(mocker):
    mocker.patch('gips.data.cdl.cdlData.__init__', _patched_cdlData_init)
    mocker.patch('osgeo.gdal.Open', _patched_gdal_open)
    d = cdl.cdlData()

    expected_legend = [''] * 256
    expected_legend[0] = 'class1'
    expected_legend[1] = 'class2'
    assert d.legend() == expected_legend and \
        d.get_code('class1') == 0 and \
        d.get_cropname(0) == 'class1'
