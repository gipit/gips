import os
import datetime

from gips.data.modis.modis import modisAsset

path_prefix = modisAsset.Repository.data_path()

asset_filenames = [
    os.path.join(path_prefix, 'h12v04/2012338/MYD11A1.A2012338.h12v04.005.2012341075802.hdf'),
    os.path.join(path_prefix, 'h12v04/2012337/MOD10A1.A2012337.h12v04.005.2012340033542.hdf'),
    os.path.join(path_prefix, 'h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf'),
]

expected_assets = {
    'MYD11A1': {
        'name':   path_prefix + '/h12v04/2012338/MYD11A1.A2012338.h12v04.005.2012341075802.hdf',
        'asset':  u'MYD11A1',
        'date':   datetime.date(2012, 12, 3),
        'driver': u'modis',
        'sensor': u'MYD',
        'tile':   u'h12v04'
    },
    'MOD10A1': {
        'name':   path_prefix + '/h12v04/2012337/MOD10A1.A2012337.h12v04.005.2012340033542.hdf',
        'asset':  u'MOD10A1',
        'date':   datetime.date(2012, 12, 2),
        'driver': u'modis',
        'sensor': u'MOD',
        'tile':   u'h12v04'
    },
    'MCD43A2': {
        'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
        'asset':  u'MCD43A2',
        'date':   datetime.date(2012, 12, 1),
        'driver': u'modis',
        'sensor': u'MCD',
        'tile':   u'h12v04'
    },
}

product_filenames = [
    path_prefix + '/h12v04/2012336/h12v04_2012336_MCD_quality.tif',
    path_prefix + '/h12v04/2012337/h12v04_2012337_MOD_temp8td.tif',
    path_prefix + '/h12v04/2012338/h12v04_2012338_MCD_fsnow.tif',
]

expected_products = {
    'quality': {
        'driver':   u'modis',
        'product':  u'quality',
        'sensor':   u'MCD',
        'tile':     u'h12v04',
        'date':     datetime.date(2012, 12, 1),
        'name':     path_prefix + '/h12v04/2012336/h12v04_2012336_MCD_quality.tif',
    },
    'temp8td': {
        'driver':   u'modis',
        'product':  u'temp8td',
        'sensor':   u'MOD',
        'tile':     u'h12v04',
        'date':     datetime.date(2012, 12, 2),
        'name':     path_prefix + '/h12v04/2012337/h12v04_2012337_MOD_temp8td.tif',
    },
    'fsnow': {
        'driver':   u'modis',
        'product':  u'fsnow',
        'sensor':   u'MCD',
        'tile':     u'h12v04',
        'date':     datetime.date(2012, 12, 3),
        'name':     path_prefix + '/h12v04/2012338/h12v04_2012338_MCD_fsnow.tif',
    },
}

