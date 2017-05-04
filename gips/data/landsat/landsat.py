#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  matt.a.hanson@gmail.com
#
#    Copyright (C) 2014 Applied Geosolutions
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>
################################################################################


import os
import re
from datetime import datetime, date
import shutil
import glob
from fnmatch import fnmatchcase
import traceback
from copy import deepcopy
import commands
import tempfile
import tarfile

import numpy

import gippy
from gips import __version__ as __gips_version__
from gippy.algorithms import ACCA, Fmask, LinearTransform, Indices, AddShadowMask
from gips.data.core import Repository, Asset, Data
from gips.atmosphere import SIXS, MODTRAN
from gips.utils import VerboseOut, RemoveFiles, basename, settings
from gips import utils

from landsat_util import search, downloader

from pdb import set_trace


requirements = ['Py6S>=1.5.0']


def binmask(arr, bit):
    """ Return boolean array indicating which elements as binary have a 1 in
        a specified bit position. Input is Numpy array.
    """
    return arr & (1 << (bit - 1)) == (1 << (bit - 1))


class landsatRepository(Repository):
    """ Singleton (all class methods) to be overridden by child data classes """
    name = 'Landsat'
    description = 'Landsat 5 (TM), 7 (ETM+), 8 (OLI)'
    _tile_attribute = 'pr'

    @classmethod
    def feature2tile(cls, feature):
        tile = super(landsatRepository, cls).feature2tile(feature)
        return tile.zfill(6)


class landsatAsset(Asset):
    """ Landsat asset (original raw tar file) """
    Repository = landsatRepository

    # tassled cap coefficients for L5 and L7
    _tcapcoef = [
        [0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596],
        [-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630],
        [0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388],
        [0.0805, -0.0498, 0.1950, -0.1327, 0.5752, -0.7775],
        [-0.7252, -0.0202, 0.6683, 0.0631, -0.1494, -0.0274],
        [0.4000, -0.8172, 0.3832, 0.0602, -0.1095, 0.0985]
    ]

    # combine sensormeta with sensor
    _sensors = {
        #'LT4': {
        #    'description': 'Landsat 4',
        #},
        'LT5': {
            'description': 'Landsat 5',
            'bands': ['1', '2', '3', '4', '5', '6', '7'],
            'oldbands': ['1', '2', '3', '4', '5', '6', '7'],
            'colors': ["BLUE", "GREEN", "RED", "NIR", "SWIR1", "LWIR", "SWIR2"],
            # TODO - update bands with actual L5 values (these are L7)
            'bandlocs': [0.4825, 0.565, 0.66, 0.825, 1.65, 11.45, 2.22],
            'bandwidths': [0.065, 0.08, 0.06, 0.15, 0.2, 2.1, 0.26],
            'E': [1983, 1796, 1536, 1031, 220.0, 0, 83.44],
            'K1': [0, 0, 0, 0, 0, 607.76, 0],
            'K2': [0, 0, 0, 0, 0, 1260.56, 0],
            'tcap': _tcapcoef,
        },
        'LE7': {
            'description': 'Landsat 7',
            #bands = ['1','2','3','4','5','6_VCID_1','6_VCID_2','7','8']
            'bands': ['1', '2', '3', '4', '5', '6_VCID_1', '7'],
            'oldbands': ['1', '2', '3', '4', '5', '61', '7'],
            'colors': ["BLUE", "GREEN", "RED", "NIR", "SWIR1", "LWIR", "SWIR2"],
            'bandlocs': [0.4825, 0.565, 0.66, 0.825, 1.65, 11.45, 2.22],
            'bandwidths': [0.065, 0.08, 0.06, 0.15, 0.2, 2.1, 0.26],
            'E': [1997, 1812, 1533, 1039, 230.8, 0, 84.90],
            'K1': [0, 0, 0, 0, 0, 666.09, 0],
            'K2': [0, 0, 0, 0, 0, 1282.71, 0],
            'tcap': _tcapcoef,
        },
        'LC8': {
            'description': 'Landsat 8',
            'bands': ['1', '2', '3', '4', '5', '6', '7', '9', '10', '11'],
            'oldbands': ['1', '2', '3', '4', '5', '6', '7', '9', '10', '11'],
            'colors': ["COASTAL", "BLUE", "GREEN", "RED", "NIR",
                       "SWIR1", "SWIR2", "CIRRUS", "LWIR", "LWIR2"],
            'bandlocs': [0.443, 0.4825, 0.5625, 0.655, 0.865,
                         1.610, 2.2, 1.375, 10.8, 12.0],
            'bandwidths': [0.01, 0.0325, 0.0375, 0.025, 0.02,
                           0.05, 0.1, 0.015, 0.5, 0.5],
            'E': [2638.35, 2031.08, 1821.09, 2075.48, 1272.96,
                  246.94, 90.61, 369.36, 0, 0],
            'K1': [0, 0, 0, 0, 0,
                   0, 0, 0, 774.89, 480.89],
            'K2': [0, 0, 0, 0, 0,
                   0, 0, 0, 1321.08, 1201.14],
            'tcap': [
                [0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872],
                [-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608],
                [0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559],
                [-0.8239, 0.0849, 0.4396, -0.058, 0.2013, -0.2773],
                [-0.3294, 0.0557, 0.1056, 0.1855, -0.4349, 0.8085],
                [0.1079, -0.9023, 0.4119, 0.0575, -0.0259, 0.0252],
            ]
        },
        'LC8SR': {
            'description': 'Landsat 8 Surface Reflectance',
        }

    }

    # TODO - consider assets and sensors relationship ?
    _assets = {
        'DN': {
            'pattern': 'L' + (20 * '?') + '.tar.gz',
        },
        'SR': {
            'pattern': 'L*-SC*.tar.gz',
        }
    }

    _defaultresolution = [30.0, 30.0]

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(landsatAsset, self).__init__(filename)

        fname = os.path.basename(filename)

        VerboseOut(("fname", fname), 2)

        self.tile = fname[3:9]
        year = fname[9:13]
        doy = fname[13:16]
        self.date = datetime.strptime(year + doy, "%Y%j")

        if fnmatchcase(fname, self._assets['SR']['pattern']):
            VerboseOut('SR asset', 2)
            self.asset = 'SR'
            self.sensor = 'LC8SR'
            self.version = int(fname[20:22])
        elif fnmatchcase(fname, self._assets['DN']['pattern']):
            VerboseOut('DN asset', 2)
            self.asset = 'DN'
            self.sensor = fname[0:3]
            self.version = int(fname[19:21])
            # Landsat DN specific additions
            smeta = self._sensors[self.sensor]
            self.meta = {}
            for i, band in enumerate(smeta['colors']):
                wvlen = smeta['bandlocs'][i]
                self.meta[band] = {
                    'bandnum': i + 1,
                    'wvlen': wvlen,
                    'wvlen1': wvlen - smeta['bandwidths'][i] / 2.0,
                    'wvlen2': wvlen + smeta['bandwidths'][i] / 2.0,
                    'E': smeta['E'][i],
                    'K1': smeta['K1'][i],
                    'K2': smeta['K2'][i],
                }
            self.visbands = [col for col in smeta['colors'] if col[0:4] != "LWIR"]
            self.lwbands = [col for col in smeta['colors'] if col[0:4] == "LWIR"]
        else:
            msg = "No matching landsat asset type for '{}'".format(fname)
            raise RuntimeError(msg, filename)

        if self.sensor not in self._sensors.keys():
            raise Exception("Sensor %s not supported: %s" % (self.sensor, filename))

    @classmethod
    def fetch(cls, asset, tile, date):

        # 'SR' not fetchable at the moment
        if asset == 'SR':
            VerboseOut('SR assets are never fetchable', 4)
            return
        paths_rows = tile[:3] + "," + tile[3:]
        fdate = date.strftime('%Y-%m-%d')

        s = search.Search()
        response = s.search(paths_rows=paths_rows, start_date=fdate, end_date=fdate, cloud_max=90)
        if response['status'] == 'SUCCESS' and response['total'] > 0:
            VerboseOut('Fetching %s %s %s' % (asset, tile, fdate), 1)
            if response['total_returned'] != 1:
                raise Exception('Single date, single location, returned more than one result')
            result = response['results'][0]
            cloudpct = result['cloud']
            sceneID = result['sceneID']
            stage_dir = os.path.join(cls.Repository.path(), 'stage')
            sceneIDs = [str(sceneID)]
            d = downloader.Downloader(download_dir=stage_dir)
            d.download(sceneIDs)
            # do the following because the downloaded .bz file has owner/group
            # settings that cause the GDAL virtual filesystem access to fail
            bz_path = glob.glob(os.path.join(stage_dir, sceneID + '*'))[0]
            gz_path = os.path.splitext(bz_path)[0] + ".gz"
            cmd = "tar xvfj %s -C %s |xargs tar cvfz %s -C %s" % (bz_path, stage_dir, gz_path, stage_dir)
            VerboseOut("Reformatting bz->gz", 1)
            result = commands.getstatusoutput(cmd)
            VerboseOut("removing %s" % bz_path)
            bands_path = glob.glob(os.path.join(stage_dir, sceneID + '_*.*'))
            # clean up - the .tar.gz will get moved on archive
            os.remove(bz_path)
            for band_path in bands_path:
                os.remove(band_path)

    def updated(self, newasset):
        '''
        Compare the version for this to that of newasset.
        Return true if newasset version is greater.
        '''
        return (self.sensor == newasset.sensor and
                self.tile == newasset.tile and
                self.date == newasset.date and
                self.version < newasset.version)


class landsatData(Data):
    name = 'Landsat'
    version = '0.9.0'

    Asset = landsatAsset

    _lt5_startdate = date(1984, 3, 1)
    _lc8_startdate = date(2013, 5, 30)

    # Group products belong to ('Standard' if not specified)
    _productgroups = {
        'Index': ['bi', 'evi', 'lswi', 'msavi2', 'ndsi', 'ndvi', 'ndwi', 'satvi'],
        'Tillage': ['ndti', 'crc', 'sti', 'isti'],
        'LC8SR': ['ndvi8sr'],
        'ACOLITE': [
            'rhow', 'oc2chl', 'oc3chl', 'fai', 'spm655', 'turbidity',
            'acoflags',
            # 'rhoam',  # Dropped for the moment due to issues in ACOLITE
        ],
    }
    __toastring = 'toa: use top of the atmosphere reflectance'
    _products = {
        #'Standard':
        'rad': {
            'assets': ['DN'],
            'description': 'Surface-leaving radiance',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'ref': {
            'assets': ['DN'],
            'description': 'Surface reflectance',
            'arguments': [__toastring],
	    'startdate': _lt5_startdate,
            'latency': 1,
        },
        'temp': {
            'assets': ['DN'],
            'description': 'Brightness (apparent) temperature',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'acca': {
            'assets': ['DN'],
            'description': 'Automated Cloud Cover Assessment',
            'arguments': [
                'X: erosion kernel diameter in pixels (default: 5)',
                'Y: dilation kernel diameter in pixels (default: 10)',
                'Z: cloud height in meters (default: 4000)'
            ],
            'nargs': '*',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'fmask': {
            'assets': ['DN'],
            'description': 'Fmask cloud cover',
            'nargs': '*',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'tcap': {
            'assets': ['DN'],
            'description': 'Tassled cap transformation',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'dn': {
            'assets': ['DN'],
            'description': 'Raw digital numbers',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'volref': {
            'assets': ['DN'],
            'description': 'Volumetric water reflectance - valid for water only',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'wtemp': {
            'assets': ['DN'],
            'description': 'Water temperature (atmospherically correct) - valid for water only',
            # It's not really TOA, but the product code will take care of atm correction itself
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'bqa': {
            'assets': ['DN'],
            'description': 'LC8 band quality',
            'toa': True,
            'startdate': _lc8_startdate,
            'latency': 1,
        },
        'bqashadow': {
            'assets': ['DN'],
            'description': 'LC8 QA + Shadow Smear',
            'arguments': [
                'X: erosion kernel diameter in pixels (default: 5)',
                'Y: dilation kernel diameter in pixels (default: 10)',
                'Z: cloud height in meters (default: 4000)'
            ],
            'nargs': '*',
            'toa': True,
            'startdate': _lc8_startdate,
            'latency': 1,
        },
        #'Indices': {
        'bi': {
            'assets': ['DN'],
            'description': 'Brightness Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'evi': {
            'assets': ['DN'],
            'description': 'Enhanced Vegetation Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'lswi': {
            'assets': ['DN'],
            'description': 'Land Surface Water Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'msavi2': {
            'assets': ['DN'],
            'description': 'Modified Soil-Adjusted Vegetation Index (revised)',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'ndsi': {
            'assets': ['DN'],
            'description': 'Normalized Difference Snow Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'ndvi': {
            'assets': ['DN'],
            'description': 'Normalized Difference Vegetation Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'ndwi': {
            'assets': ['DN'],
            'description': 'Normalized Difference Water Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'satvi': {
            'assets': ['DN'],
            'description': 'Soil-Adjusted Total Vegetation Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        #'Tillage Indices': {
        'ndti': {
            'assets': ['DN'],
            'description': 'Normalized Difference Tillage Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'crc': {
            'assets': ['DN'],
            'description': 'Crop Residue Cover',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'sti': {
            'assets': ['DN'],
            'description': 'Standard Tillage Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'isti': {
            'assets': ['DN'],
            'description': 'Inverse Standard Tillage Index',
            'arguments': [__toastring],
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        # NEW!!!
        'ndvi8sr': {
            'assets': ['SR'],
            'description': 'Normalized Difference Vegetation from LC8SR',
            'startdate': _lc8_startdate,
            'latency': 1,
        },
        'landmask': {
            'assets': ['SR'],
            'description': 'Land mask from LC8SR',
            'startdate': _lc8_startdate,
            'latency': 1,
        },
        # ACOLITE products
        'rhow': {
            'assets': ['DN'],
            'description': 'Water-Leaving Radiance-Reflectance',
            'acolite-product': 'rhow_vnir',
            'acolite-key': 'RHOW',
            'gain': 0.0001,
            'offset': 0.,
            'dtype': 'int16',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        # Not sure what the issue is with this product, but it doesn't seem to
        # work as expected (multiband vis+nir product)
        # 'rhoam': {
        #     'assets': ['DN'],
        #     'description': 'Multi-Scattering Aerosol Reflectance',
        #     'acolite-product': 'rhoam_vnir',
        #     'acolite-key': 'RHOAM',
        #     'dtype': 'int16',
        #     'toa': True,
        # },
        'oc2chl': {
            'assets': ['DN'],
            'description': 'Blue-Green Ratio Chlorophyll Algorithm using bands 483 & 561',
            'acolite-product': 'CHL_OC2',
            'acolite-key': 'CHL_OC2',
            'gain': 0.0125,
            'offset': 250.,
            'dtype': 'int16',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'oc3chl': {
            'assets': ['DN'],
            'description': 'Blue-Green Ratio Chlorophyll Algorithm using bands 443, 483, & 561',
            'acolite-product': 'CHL_OC3',
            'acolite-key': 'CHL_OC3',
            'gain': 0.0125,
            'offset': 250.,
            'dtype': 'int16',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'fai': {
            'assets': ['DN'],
            'description': 'Floating Algae Index',
            'acolite-product': 'FAI',
            'acolite-key': 'FAI',
            'dtype': 'float32',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'acoflags': {
            'assets': ['DN'],
            'description': '0 = water 1 = no data 2 = land',
            'acolite-product': 'FLAGS',
            'acolite-key': 'FLAGS',
            'dtype': 'uint8',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'spm655': {
            'assets': ['DN'],
            'description': 'Suspended Sediment Concentration 655nm',
            'acolite-product': 'SPM_NECHAD_655',
            'acolite-key': 'SPM_NECHAD_655',
            'offset': 50.,
            'gain': 0.005,
            'dtype': 'int16',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
        'turbidity': {
            'assets': ['DN'],
            'description': 'Blended Turbidity',
            'acolite-product': 'T_DOGLIOTTI',
            'acolite-key': 'T_DOGLIOTTI',
            'offset': 50.,
            'gain': 0.005,
            'dtype': 'int16',
            'toa': True,
            'startdate': _lt5_startdate,
            'latency': 1,
        },
    }

    def _process_acolite(self, asset, aco_proc_dir, products):
        '''
        TODO: Move this to `gips.atmosphere`.
        TODO: Ensure this is genericized to work for S2 or Landsat.
        '''
        import netCDF4
        ACOLITEPATHS = {
            'ACO_DIR': settings().REPOS['landsat']['ACOLITE_DIR'],
            # N.B.: only seems to work when run from the ACO_DIR
            'IDLPATH': 'idl',
            'ACOLITE_BINARY': 'acolite.sav',
            # TODO: template may be the only piece that needs
            #       to be moved for driver-independence.
            'SETTINGS_TEMPLATE': os.path.join(
                os.path.dirname(__file__),
                'acolite.cfg'
            )
        }
        ACOLITE_NDV = 1.875 * 2 ** 122
        # mapping from dtype to gdal type and nodata value
        IMG_PARAMS = {
            'float32': (gippy.GDT_Float32, -32768.),
            'int16': (gippy.GDT_Int16, -32768),
            'uint8': (gippy.GDT_Byte, 1),
        }
        imeta = products.pop('meta')

        # TODO: add 'outdir' to `gips.data.core.Asset.extract` method
        # EXTRACT ASSET
        tar = tarfile.open(asset.filename)
        tar.extractall(aco_proc_dir)

        # STASH PROJECTION AND GEOTRANSFORM (in a GeoImage)
        exts = re.compile(r'.*\.((jp2)|(tif)|(TIF))$')
        tif = filter(
            lambda de: exts.match(de),
            os.listdir(aco_proc_dir)
        )[0]
        tmp = gippy.GeoImage(os.path.join(aco_proc_dir, tif))

        # PROCESS SETTINGS TEMPLATE FOR SPECIFIED PRODUCTS
        settings_path = os.path.join(aco_proc_dir, 'settings.cfg')
        template_path = ACOLITEPATHS.pop('SETTINGS_TEMPLATE')
        acolite_products = ','.join(
            [
                products[k]['acolite-product']
                for k in products
                if k != 'acoflags'  # acoflags is always internally generated
                                    # by ACOLITE,
            ]
        )
        if len(acolite_products) == 0:
            raise Exception(
                "ACOLITE: Must specify at least 1 product.\n"
                "'acoflags' cannot be generated on its own.",
            )
        with open(template_path, 'r') as aco_template:
            with open(settings_path, 'w') as aco_settings:
                for line in aco_template:
                    aco_settings.write(
                        re.sub(
                            r'GIPS_LANDSAT_PRODUCTS',
                            acolite_products,
                            line
                        )
                    )
        ACOLITEPATHS['ACOLITE_SETTINGS'] = settings_path

        # PROCESS VIA ACOLITE IDL CALL
        cmd = (
            ('cd {ACO_DIR} ; '
             '{IDLPATH} -IDL_CPU_TPOOL_NTHREADS 1 '
             '-rt={ACOLITE_BINARY} '
             '-args settings={ACOLITE_SETTINGS} '
             'run=1 '
             'output={OUTPUT} image={IMAGES}')
            .format(
                OUTPUT=aco_proc_dir,
                IMAGES=aco_proc_dir,
                **ACOLITEPATHS
            )
        )
        utils.verbose_out('Running: {}'.format(cmd), 2)
        status, output = commands.getstatusoutput(cmd)
        if status != 0:
            raise Exception(cmd, output)
        aco_nc_file = glob.glob(os.path.join(aco_proc_dir, '*_L2.nc'))[0]
        dsroot = netCDF4.Dataset(aco_nc_file)

        # EXTRACT IMAGES FROM NETCDF AND
        # COMBINE MULTI-IMAGE PRODUCTS INTO
        # A MULTI-BAND TIF, ADD METADATA, and MOVE INTO TILES
        prodout = dict()

        for key in products:
            ofname = products[key]['fname']
            aco_key = products[key]['acolite-key']
            bands = list(filter(
                lambda x: str(x) == aco_key or x.startswith(aco_key),
                dsroot.variables.keys()
            ))
            npdtype = products[key]['dtype']
            dtype, missing = IMG_PARAMS[npdtype]
            gain = products[key].get('gain', 1.0)
            offset = products[key].get('offset', 0.0)
            imgout = gippy.GeoImage(ofname, tmp, dtype, len(bands))
            # # TODO: add units to products dictionary and use here.
            # imgout.SetUnits(products[key]['units'])
            pmeta = dict()
            pmeta.update(imeta)
            pmeta = {
                mdi: products[key][mdi]
                for mdi in ['acolite-key', 'description']
            }
            pmeta['source_asset'] = os.path.basename(asset.filename)
            imgout.SetMeta(pmeta)
            for i, b in enumerate(bands):
                imgout.SetBandName(str(b), i + 1)

            for i, b in enumerate(bands):
                var = dsroot.variables[b][:]
                arr = numpy.array(var)
                if hasattr(dsroot.variables[b], '_FillValue'):
                    fill = dsroot.variables[b]._FillValue
                else:
                    fill = ACOLITE_NDV
                mask = arr != fill
                arr[numpy.invert(mask)] = missing
                # if key == 'rhow':
                #     set_trace()
                arr[mask] = ((arr[mask] - offset) / gain)
                imgout[i].Write(arr.astype(npdtype))

            prodout[key] = imgout.Filename()
            imgout = None
            imgout = gippy.GeoImage(ofname, True)
            imgout.SetGain(gain)
            imgout.SetOffset(offset)
            imgout.SetNoData(missing)
        return prodout


    def process(self, products=None, overwrite=False, **kwargs):
        """ Make sure all products have been processed """
        products = super(landsatData, self).process(products, overwrite, **kwargs)
        if len(products) == 0:
            return

        start = datetime.now()

        assets = set()
        for key, val in products.requested.items():
            assets.update(self._products[val[0]]['assets'])

        if len(assets) != 1:
            raise Exception('This driver does not support creation of products'
                            ' from different Assets at the same time')

        asset = list(assets)[0]

        # TODO: De-hack this
        # Better approach, but needs some thought, is to loop over assets
        # Ian, you are right. I just don't have enough time to do it.

        if asset == 'SR':

            datafiles = self.assets['SR'].datafiles()

            imgpaths = dict()

            for datafile in datafiles:

                key = datafile.partition('_')[2].split('.')[0]
                path = os.path.join('/vsitar/' + self.assets['SR'].filename, datafile)

                imgpaths[key] = path

            # print imgpaths

            bname = os.path.join(self.path, self.basename)

            for key, val in products.requested.items():

                if val[0] == "ndvi8sr":

                    sensor = 'LC8SR'
                    fname = '%s_%s_%s' % (bname, sensor, key)

                    img = gippy.GeoImage([imgpaths['sr_band4'], imgpaths['sr_band5']])

                    missing = float(img[0].NoDataValue())

                    red = img[0].Read().astype('float32')
                    nir = img[1].Read().astype('float32')

                    wvalid = numpy.where((red != missing) & (nir != missing) & (red + nir != 0.0))

                    red[wvalid] *= 1.E-4
                    nir[wvalid] *= 1.E-4

                    # TODO: change this so that these pixels become missing
                    red[(red != missing) & (red < 0.0)] = 0.0
                    red[red > 1.0] = 1.0
                    nir[(nir != missing) & (nir < 0.0)] = 0.0
                    nir[nir > 1.0] = 1.0

                    ndvi = missing + numpy.zeros_like(red)
                    ndvi[wvalid] = (nir[wvalid] - red[wvalid])/(nir[wvalid] + red[wvalid])

                    VerboseOut("writing " + fname, 2)
                    imgout = gippy.GeoImage(fname, img, gippy.GDT_Float32, 1)
                    imgout.SetNoData(-9999.)
                    imgout.SetOffset(0.0)
                    imgout.SetGain(1.0)
                    imgout.SetBandName('NDVI', 1)
                    imgout[0].Write(ndvi)

                if val[0] == "landmask":

                    sensor = 'LC8SR'
                    fname = '%s_%s_%s' % (bname, sensor, key)

                    img = gippy.GeoImage([imgpaths['cfmask'], imgpaths['cfmask_conf']])

                    cfmask = img[0].Read()
                    # array([  0,   1,   2,   3,   4, 255], dtype=uint8)
                    # 0 means clear! but I want 1 to mean clear

                    cfmask[cfmask > 0] = 2
                    cfmask[cfmask == 0] = 1
                    cfmask[cfmask == 2] = 0

                    VerboseOut("writing " + fname, 2)
                    imgout = gippy.GeoImage(fname, img, gippy.GDT_Byte, 1)
                    imgout.SetBandName('Land mask', 1)
                    imgout[0].Write(cfmask)


        elif asset == 'DN':

            # This block contains everything that existed in the first generation Landsat driver

            # Add the sensor for this date to the basename
            self.basename = self.basename + '_' + self.sensors[asset]

            # Read the assets
            with utils.error_handler('Error reading ' + basename(self.assets['DN'].filename)):
                img = self._readraw()

            meta = self.assets['DN'].meta
            visbands = self.assets['DN'].visbands
            lwbands = self.assets['DN'].lwbands
            md = self.meta_dict()

            # running atmosphere if any products require it
            toa = True
            for val in products.requested.values():
                toa = toa and (self._products[val[0]].get('toa', False) or 'toa' in val)
            if not toa:
                start = datetime.now()

                if not settings().REPOS[self.Repository.name.lower()]['6S']:
                    raise Exception('6S is required for atmospheric correction')
                with utils.error_handler('Problem running 6S atmospheric model'):
                    wvlens = [(meta[b]['wvlen1'], meta[b]['wvlen2']) for b in visbands]
                    geo = self.metadata['geometry']
                    atm6s = SIXS(visbands, wvlens, geo, self.metadata['datetime'], sensor=self.sensor_set[0])
                    md["AOD Source"] = str(atm6s.aod[0])
                    md["AOD Value"] = str(atm6s.aod[1])

            # Break down by group
            groups = products.groups()
            # ^--- has the info about what products the user requested

            # create non-atmospherically corrected apparent reflectance and temperature image
            reflimg = gippy.GeoImage(img)
            theta = numpy.pi * self.metadata['geometry']['solarzenith'] / 180.0
            sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (float(self.day) - 4.0) / 180.0))
            for col in self.assets['DN'].visbands:
                reflimg[col] = img[col] * (1.0 /
                        ((meta[col]['E'] * numpy.cos(theta)) / (numpy.pi * sundist * sundist)))
            for col in self.assets['DN'].lwbands:
                reflimg[col] = (((img[col].pow(-1)) * meta[col]['K1'] + 1).log().pow(-1)
                        ) * meta[col]['K2'] - 273.15

            # This is landsat, so always just one sensor for a given date
            sensor = self.sensors['DN']

            # Process standard products
            for key, val in groups['Standard'].items():
                start = datetime.now()
                # TODO - update if no atmos desired for others
                toa = self._products[val[0]].get('toa', False) or 'toa' in val
                # Create product
                with utils.error_handler(
                        'Error creating product {} for {}'
                        .format(key, basename(self.assets['DN'].filename)),
                        continuable=True):
                    fname = os.path.join(self.path, self.basename + '_' + key)
                    if val[0] == 'acca':
                        s_azim = self.metadata['geometry']['solarazimuth']
                        s_elev = 90 - self.metadata['geometry']['solarzenith']
                        erosion, dilation, cloudheight = 5, 10, 4000
                        if len(val) >= 4:
                            erosion, dilation, cloudheight = [int(v) for v in val[1:4]]
                        resset = set(
                            [(reflimg[band].Resolution().x(),
                              reflimg[band].Resolution().y())
                             for band in (self.assets['DN'].visbands +
                                          self.assets['DN'].lwbands)]
                        )
                        if len(resset) > 1:
                            raise Exception(
                                'ACCA requires all bands to have the same '
                                'spatial resolution.  Found:\n\t' + str(resset)
                            )
                        imgout = ACCA(reflimg, fname, s_elev, s_azim, erosion, dilation, cloudheight)
                    elif val[0] == 'fmask':
                        tolerance, dilation = 3, 5
                        if len(val) >= 3:
                            tolerance, dilation = [int(v) for v in val[1:3]]
                        imgout = Fmask(reflimg, fname, tolerance, dilation)
                    elif val[0] == 'rad':
                        imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(visbands))
                        for i in range(0, imgout.NumBands()):
                            imgout.SetBandName(visbands[i], i + 1)
                        imgout.SetNoData(-32768)
                        imgout.SetGain(0.1)
                        if toa:
                            for col in visbands:
                                img[col].Process(imgout[col])
                        else:
                            for col in visbands:
                                ((img[col] - atm6s.results[col][1]) / atm6s.results[col][0]).Process(imgout[col])
                        # Mask out any pixel for which any band is nodata
                        #imgout.ApplyMask(img.DataMask())
                    elif val[0] == 'ref':
                        imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(visbands))
                        for i in range(0, imgout.NumBands()):
                            imgout.SetBandName(visbands[i], i + 1)
                        imgout.SetNoData(-32768)
                        imgout.SetGain(0.0001)
                        if toa:
                            for c in visbands:
                                reflimg[c].Process(imgout[c])
                        else:
                            for c in visbands:
                                (((img[c] - atm6s.results[c][1]) / atm6s.results[c][0]) * (1.0 / atm6s.results[c][2])).Process(imgout[c])
                        # Mask out any pixel for which any band is nodata
                        #imgout.ApplyMask(img.DataMask())
                    elif val[0] == 'tcap':
                        tmpimg = gippy.GeoImage(reflimg)
                        tmpimg.PruneBands(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'])
                        arr = numpy.array(self.Asset._sensors[self.sensor_set[0]]['tcap']).astype('float32')
                        imgout = LinearTransform(tmpimg, fname, arr)
                        outbands = ['Brightness', 'Greenness', 'Wetness', 'TCT4', 'TCT5', 'TCT6']
                        for i in range(0, imgout.NumBands()):
                            imgout.SetBandName(outbands[i], i + 1)
                    elif val[0] == 'temp':
                        imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(lwbands))
                        for i in range(0, imgout.NumBands()):
                            imgout.SetBandName(lwbands[i], i + 1)
                        imgout.SetNoData(-32768)
                        imgout.SetGain(0.1)
                        [reflimg[col].Process(imgout[col]) for col in lwbands]
                    elif val[0] == 'dn':
                        rawimg = self._readraw()
                        rawimg.SetGain(1.0)
                        rawimg.SetOffset(0.0)
                        imgout = rawimg.Process(fname)
                        rawimg = None
                    elif val[0] == 'volref':
                        bands = deepcopy(visbands)
                        bands.remove("SWIR1")
                        imgout = gippy.GeoImage(fname, reflimg, gippy.GDT_Int16, len(bands))
                        [imgout.SetBandName(band, i + 1) for i, band in enumerate(bands)]
                        imgout.SetNoData(-32768)
                        imgout.SetGain(0.0001)
                        r = 0.54    # Water-air reflection
                        p = 0.03    # Internal Fresnel reflectance
                        pp = 0.54   # Water-air Fresnel reflectance
                        n = 1.34    # Refractive index of water
                        Q = 1.0     # Downwelled irradiance / upwelled radiance
                        A = ((1 - p) * (1 - pp)) / (n * n)
                        srband = reflimg['SWIR1'].Read()
                        nodatainds = srband == reflimg['SWIR1'].NoDataValue()
                        for band in bands:
                            bimg = reflimg[band].Read()
                            diffimg = bimg - srband
                            diffimg = diffimg / (A + r * Q * diffimg)
                            diffimg[bimg == reflimg[band].NoDataValue()] = imgout[band].NoDataValue()
                            diffimg[nodatainds] = imgout[band].NoDataValue()
                            imgout[band].Write(diffimg)
                    elif val[0] == 'wtemp':
                        imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(lwbands))
                        [imgout.SetBandName(lwbands[i], i + 1) for i in range(0, imgout.NumBands())]
                        imgout.SetNoData(-32768)
                        imgout.SetGain(0.1)
                        tmpimg = gippy.GeoImage(img)
                        for col in lwbands:
                            band = tmpimg[col]
                            m = meta[col]
                            lat = self.metadata['geometry']['lat']
                            lon = self.metadata['geometry']['lon']
                            dt = self.metadata['datetime']
                            atmos = MODTRAN(m['bandnum'], m['wvlen1'], m['wvlen2'], dt, lat, lon, True)
                            e = 0.95
                            band = (tmpimg[col] - (atmos.output[1] + (1 - e) * atmos.output[2])) / (atmos.output[0] * e)
                            band = (((band.pow(-1)) * meta[col]['K1'] + 1).log().pow(-1)) * meta[col]['K2'] - 273.15
                            band.Process(imgout[col])

                    elif val[0] == 'bqa':
                        if 'LC8' not in self.sensor_set:
                            continue
                        imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, 7)
                        qaimg = self._readqa()
                        qadata = qaimg.Read()
                        notfilled = ~binmask(qadata, 1)
                        notdropped = ~binmask(qadata, 2)
                        notterrain = ~binmask(qadata, 3)
                        notcirrus = ~binmask(qadata, 14) & binmask(qadata, 13)
                        notcloud = ~binmask(qadata, 16) & binmask(qadata, 15)
                        allgood = notfilled * notdropped * notterrain * notcirrus * notcloud
                        imgout[0].Write(allgood.astype('int16'))
                        imgout[1].Write(notfilled.astype('int16'))
                        imgout[2].Write(notdropped.astype('int16'))
                        imgout[3].Write(notterrain.astype('int16'))
                        imgout[4].Write(notsnow.astype('int16'))
                        imgout[5].Write(notcirrus.astype('int16'))
                        imgout[6].Write(notcloud.astype('int16'))

                    elif val[0] == 'bqashadow':
                        if 'LC8' not in self.sensor_set:
                            continue
                        imgout = gippy.GeoImage(fname, img, gippy.GDT_UInt16, 1)
                        imgout[0].SetNoData(0)
                        qaimg = self._readqa()
                        qadata = qaimg.Read()
                        fill = binmask(qadata, 1)
                        dropped = binmask(qadata, 2)
                        terrain = binmask(qadata, 3)
                        cirrus = binmask(qadata, 14)
                        othercloud = binmask(qadata, 16)
                        cloud = (cirrus + othercloud) + 2 * (fill + dropped + terrain)
                        abfn = fname + '-intermediate'
                        abimg = gippy.GeoImage(abfn, img, gippy.GDT_UInt16, 1)
                        abimg[0].SetNoData(2)
                        abimg[0].Write(cloud.astype(numpy.uint16))
                        abimg.Process()
                        abimg = None
                        abimg = gippy.GeoImage(abfn + '.tif')

                        s_azim = self.metadata['geometry']['solarazimuth']
                        s_elev = 90 - self.metadata['geometry']['solarzenith']
                        erosion, dilation, cloudheight = 5, 10, 4000
                        if len(val) >= 4:
                            erosion, dilation, cloudheight = [int(v) for v in val[1:4]]
                        imgout = AddShadowMask(
                            abimg, imgout, 0, s_elev, s_azim, erosion,
                            dilation, cloudheight, {'notes': 'dev-version'}
                        )
                        imgout.Process()
                        abimg = None
                        os.remove(abfn + '.tif')
                    fname = imgout.Filename()
                    imgout.SetMeta(md)
                    imgout = None
                    self.AddFile(sensor, key, fname)
                    VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname), datetime.now() - start), 1)

            # Process Indices
            indices0 = dict(groups['Index'], **groups['Tillage'])
            if len(indices0) > 0:
                start = datetime.now()
                indices = {}
                indices_toa = {}
                for key, val in indices0.items():
                    if 'toa' in val:
                        indices_toa[key] = val
                    else:
                        indices[key] = val
                # Run TOA
                if len(indices_toa) > 0:
                    fnames = [os.path.join(self.path, self.basename + '_' + key) for key in indices_toa]
                    prodout = Indices(reflimg, dict(zip([p[0] for p in indices_toa.values()], fnames)), md)
                    prodout = dict(zip(indices_toa.keys(), prodout.values()))
                    [self.AddFile(sensor, key, fname) for key, fname in prodout.items()]
                # Run atmospherically corrected
                if len(indices) > 0:
                    fnames = [os.path.join(self.path, self.basename + '_' + key) for key in indices]
                    for col in visbands:
                        img[col] = ((img[col] - atm6s.results[col][1]) / atm6s.results[col][0]
                                ) * (1.0 / atm6s.results[col][2])
                    prodout = Indices(img, dict(zip([p[0] for p in indices.values()], fnames)), md)
                    prodout = dict(zip(indices.keys(), prodout.values()))
                    [self.AddFile(sensor, key, fname) for key, fname in prodout.items()]
                VerboseOut(' -> %s: processed %s in %s' % (
                        self.basename, indices0.keys(), datetime.now() - start), 1)
            img = None
            # cleanup directory
            try:
                if settings().REPOS[self.Repository.name.lower()]['extract']:
                    for bname in self.assets['DN'].datafiles():
                        if bname[-7:] != 'MTL.txt':
                            files = glob.glob(os.path.join(self.path, bname) + '*')
                            RemoveFiles(files)
                shutil.rmtree(os.path.join(self.path, 'modtran'))
            except:
                # TODO error-handling-fix: continuable handler
                # VerboseOut(traceback.format_exc(), 4)
                pass

            if groups['ACOLITE']:
                start = datetime.now()
                # TEMPDIR FOR PROCESSING
                aco_proc_dir = tempfile.mkdtemp(
                    prefix='aco_proc_',
                    dir=os.path.join(self.Repository.path(), 'stage')
                )
                with utils.error_handler(
                        'Error creating ACOLITE products {} for {}'
                        .format(
                            groups['ACOLITE'].keys(),
                            basename(self.assets['DN'].filename)
                        ),
                        continuable=True):
                    # amd is 'meta' (common to all products) and product info dicts
                    amd = {
                        'meta': md.copy()
                    }
                    for p in groups['ACOLITE']:
                        amd[p] = {
                            'fname': os.path.join(
                                self.path, self.basename + '_' + p + '.tif'
                            )
                        }
                        amd[p].update(self._products[p])
                        amd[p].pop('assets')
                    #set_trace()
                    prodout = self._process_acolite(
                        asset=self.assets['DN'],
                        aco_proc_dir=aco_proc_dir,
                        products=amd,
                    )
                    endtime = datetime.now()
                    for k, fn in prodout.items():
                        self.AddFile(sensor, k, fn)
                    VerboseOut(
                        ' -> {}: processed {} in {}'
                        .format(self.basename, prodout.keys(), endtime - start),
                        1
                    )
                shutil.rmtree(aco_proc_dir)
                ## end ACOLITE

    def filter(self, pclouds=100, sensors=None, **kwargs):
        """Check if Data object passes filter.

        User can't enter pclouds, but can pass in --sensors.  kwargs
        isn't used.
        """
        if pclouds < 100:
            self.meta()
            if self.metadata['clouds'] > pclouds:
                return False
        if sensors:
            if type(sensors) is str:
                sensors = [sensors]
            sensors = set(sensors)
            # ideally, the data class would be trimmed by
            if not sensors.intersection(self.sensor_set):
                return False
        return True

    def meta(self):
        """ Read in Landsat MTL (metadata) file """

        # test if metadata already read in, if so, return

        datafiles = self.assets['DN'].datafiles()

        # locate MTL file and save it to disk if it isn't saved already
        mtlfilename = [f for f in datafiles if 'MTL.txt' in f][0]
        if os.path.exists(mtlfilename) and os.stat(mtlfilename).st_size == 0:
            os.remove(mtlfilename)
        if not os.path.exists(mtlfilename):
            mtlfilename = self.assets['DN'].extract([mtlfilename])[0]
        # Read MTL file
        with utils.error_handler('Error reading metadata file ' + mtlfilename):
            text = open(mtlfilename, 'r').read()
        if len(text) < 10:
            raise Exception('MTL file is too short. {}'.format(mtlfilename))

        sensor = self.assets['DN'].sensor
        smeta = self.assets['DN']._sensors[sensor]

        # Process MTL text - replace old metadata tags with new
        # NOTE This is not comprehensive, there may be others
        text = text.replace('ACQUISITION_DATE', 'DATE_ACQUIRED')
        text = text.replace('SCENE_CENTER_SCAN_TIME', 'SCENE_CENTER_TIME')

        for (ob, nb) in zip(smeta['oldbands'], smeta['bands']):
            text = re.sub(r'\WLMIN_BAND' + ob, 'RADIANCE_MINIMUM_BAND_' + nb, text)
            text = re.sub(r'\WLMAX_BAND' + ob, 'RADIANCE_MAXIMUM_BAND_' + nb, text)
            text = re.sub(r'\WQCALMIN_BAND' + ob, 'QUANTIZE_CAL_MIN_BAND_' + nb, text)
            text = re.sub(r'\WQCALMAX_BAND' + ob, 'QUANTIZE_CAL_MAX_BAND_' + nb, text)
            text = re.sub(r'\WBAND' + ob + '_FILE_NAME', 'FILE_NAME_BAND_' + nb, text)
        for l in ('LAT', 'LON', 'MAPX', 'MAPY'):
            for c in ('UL', 'UR', 'LL', 'LR'):
                text = text.replace('PRODUCT_' + c + '_CORNER_' + l, 'CORNER_' + c + '_' + l + '_PRODUCT')
        text = text.replace('\x00', '')
        # Remove junk
        lines = text.split('\n')
        mtl = dict()
        for l in lines:
            meta = l.replace('\"', "").strip().split('=')
            if len(meta) > 1:
                key = meta[0].strip()
                item = meta[1].strip()
                if key != "GROUP" and key != "END_GROUP":
                    mtl[key] = item

        # Extract useful metadata
        lats = (float(mtl['CORNER_UL_LAT_PRODUCT']), float(mtl['CORNER_UR_LAT_PRODUCT']),
                float(mtl['CORNER_LL_LAT_PRODUCT']), float(mtl['CORNER_LR_LAT_PRODUCT']))
        lons = (float(mtl['CORNER_UL_LON_PRODUCT']), float(mtl['CORNER_UR_LON_PRODUCT']),
                float(mtl['CORNER_LL_LON_PRODUCT']), float(mtl['CORNER_LR_LON_PRODUCT']))
        lat = (min(lats) + max(lats)) / 2.0
        lon = (min(lons) + max(lons)) / 2.0
        dt = datetime.strptime(mtl['DATE_ACQUIRED'] + ' ' + mtl['SCENE_CENTER_TIME'][:-2], '%Y-%m-%d %H:%M:%S.%f')
        try:
            clouds = float(mtl['CLOUD_COVER'])
        except:
            # TODO error-handling-fix: no try-except needed
            clouds = 0

        filenames = []
        gain = []
        offset = []
        dynrange = []
        for i, b in enumerate(smeta['bands']):
            minval = int(float(mtl['QUANTIZE_CAL_MIN_BAND_' + b]))
            maxval = int(float(mtl['QUANTIZE_CAL_MAX_BAND_' + b]))
            minrad = float(mtl['RADIANCE_MINIMUM_BAND_' + b])
            maxrad = float(mtl['RADIANCE_MAXIMUM_BAND_' + b])
            gain.append((maxrad - minrad) / (maxval - minval))
            offset.append(minrad)
            dynrange.append((minval, maxval))
            filenames.append(mtl['FILE_NAME_BAND_' + b].strip('\"'))

        _geometry = {
            'solarzenith': (90.0 - float(mtl['SUN_ELEVATION'])),
            'solarazimuth': float(mtl['SUN_AZIMUTH']),
            'zenith': 0.0,
            'azimuth': 180.0,
            'lat': lat,
            'lon': lon,
        }

        try:
            qafilename = [f for f in datafiles if '_BQA.TIF' in f][0]
        except Exception:
            # TODO error-handling-fix: no try-except needed
            qafilename = None

        self.metadata = {
            'filenames': filenames,
            'gain': gain,
            'offset': offset,
            'dynrange': dynrange,
            'geometry': _geometry,
            'datetime': dt,
            'clouds': clouds,
            'qafilename': qafilename
        }
        #self.metadata.update(smeta)

    @classmethod
    def meta_dict(cls):
        meta = super(landsatData, cls).meta_dict()
        meta['GIPS-landsat Version'] = cls.version
        return meta

    def _readqa(self):
        # make sure metadata is loaded
        if not hasattr(self, 'metadata'):
            self.meta()
        if settings().REPOS[self.Repository.name.lower()]['extract']:
            # Extract files
            qadatafile = self.assets['DN'].extract([self.metadata['qafilename']])
        else:
            # Use tar.gz directly using GDAL's virtual filesystem
            qadatafile = os.path.join('/vsitar/' + self.assets['DN'].filename, self.metadata['qafilename'])
        qaimg = gippy.GeoImage(qadatafile)
        return qaimg


    def _readraw(self):
        """ Read in Landsat bands using original tar.gz file """
        start = datetime.now()
        # make sure metadata is loaded
        if not hasattr(self, 'metadata'):
            self.meta()

        if settings().REPOS[self.Repository.name.lower()]['extract']:
            # Extract all files
            datafiles = self.assets['DN'].extract(self.metadata['filenames'])
        else:
            # Use tar.gz directly using GDAL's virtual filesystem
            datafiles = [os.path.join('/vsitar/' + self.assets['DN'].filename, f)
                    for f in self.metadata['filenames']]

        image = gippy.GeoImage(datafiles)
        image.SetNoData(0)

        # TODO - set appropriate metadata
        #for key,val in meta.iteritems():
        #    image.SetMeta(key,str(val))

        # Geometry used for calculating incident irradiance
        # colors = self.assets['DN']._sensors[self.sensor_set[0]]['colors']

        sensor = self.assets['DN'].sensor
        colors = self.assets['DN']._sensors[sensor]['colors']

        for bi in range(0, len(self.metadata['filenames'])):
            image.SetBandName(colors[bi], bi + 1)
            # need to do this or can we index correctly?
            band = image[bi]
            gain = self.metadata['gain'][bi]
            band.SetGain(gain)
            band.SetOffset(self.metadata['offset'][bi])
            dynrange = self.metadata['dynrange'][bi]
            # #band.SetDynamicRange(dynrange[0], dynrange[1])
            # dynrange[0] was used internally to for conversion to radiance
            # from DN in GeoRaster.Read:
            #   img = Gain() * (img-_minDC) + Offset();  # (1)
            # and with the removal of _minDC and _maxDC it is now:
            #   img = Gain() * img + Offset();           # (2)
            # And 1 can be re-written as:
            #   img = Gain() * img - Gain() *
            #                       _minDC + Offset();   # (3)
            #       = Gain() * img + Offset
            #                      - _min * Gain() ;     # (4)
            # So, since the gippy now has line (2), we can add
            # the final term of (4) [as below] to keep that functionality.
            image[bi] = band - dynrange[0] * gain
            # I verified this by example.  With old gippy.GeoRaster:
            #     In [8]: a.min()
            #     Out[8]: -64.927711
            # with new version.
            #     In [20]: ascale.min() - 1*0.01298554277169103
            #     Out[20]: -64.927711800095906


        VerboseOut('%s: read in %s' % (image.Basename(), datetime.now() - start), 2)
        return image
