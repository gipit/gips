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
from datetime import datetime
import shutil
import numpy
import glob
import traceback
from copy import deepcopy

import gippy
from gippy.algorithms import ACCA, Fmask, LinearTransform, Indices
from gips.data.core import Repository, Asset, Data
from gips.atmosphere import SIXS, MODTRAN
from gips.utils import VerboseOut, RemoveFiles, basename, settings

requirements = ['Py6S>=1.5.0']


class LandsatRepository(Repository):
    """ Singleton (all class methods) to be overridden by child data classes """
    name = 'Landsat'
    description = 'Landsat 5 (TM), 7 (ETM+), 8 (OLI)'

    @classmethod
    def feature2tile(cls, feature):
        tile = super(LandsatRepository, cls).feature2tile(feature)
        return tile.zfill(6)


class LandsatAsset(Asset):
    """ Landsat asset (original raw tar file) """
    Repository = LandsatRepository

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
            'colors': ["COASTAL", "BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2", "CIRRUS", "LWIR", "LWIR2"],
            'bandlocs': [0.443, 0.4825, 0.5625, 0.655, 0.865, 1.610, 2.2, 1.375, 10.8, 12.0],
            'bandwidths': [0.01, 0.0325, 0.0375, 0.025, 0.02, 0.05, 0.1, 0.015, 0.5, 0.5],
            'E': [2638.35, 2031.08, 1821.09, 2075.48, 1272.96, 246.94, 90.61, 369.36, 0, 0],
            'K1': [0, 0, 0, 0, 0, 0, 0, 0, 774.89, 480.89],
            'K2': [0, 0, 0, 0, 0, 0, 0, 0, 1321.08, 1201.14],
            'tcap': [
                [0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872],
                [-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608],
                [0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559],
                [-0.8239, 0.0849, 0.4396, -0.058, 0.2013, -0.2773],
                [-0.3294, 0.0557, 0.1056, 0.1855, -0.4349, 0.8085],
                [0.1079, -0.9023, 0.4119, 0.0575, -0.0259, 0.0252],
            ]
        }
    }

    # TODO - consider assets and sensors relationship ?
    _assets = {
        '': {
            'pattern': 'L*.tar.gz'
            #_pattern = r'^L[TEC][4578].*\.tar\.gz$'
        }
    }

    _defaultresolution = [30.0, 30.0]

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(LandsatAsset, self).__init__(filename)
        fname = os.path.basename(filename)
        self.sensor = fname[0:3]
        self.tile = fname[3:9]
        year = fname[9:13]
        doy = fname[13:16]
        self.date = datetime.strptime(year + doy, "%Y%j")
        if self.sensor not in self._sensors.keys():
            raise Exception("Sensor %s not supported: %s" % (self.sensor, filename))
        # Landsat specific additions
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


class LandsatData(Data):
    name = 'Landsat'
    version = '0.9.0'

    Asset = LandsatAsset

    _prodpattern = '*.tif'
    # Group products belong to ('Standard' if not specified)
    _productgroups = {
        'Index': ['bi', 'evi', 'lswi', 'msavi2', 'ndsi', 'ndvi', 'ndwi', 'satvi'],
        'Tillage': ['ndti', 'crc', 'sti', 'isti'],
    }
    __toastring = 'toa: use top of the atmosphere reflectance'
    _products = {
        #'Standard':
        'rad': {
            'description': 'Surface-leaving radiance',
            'arguments': [__toastring]
        },
        'ref': {
            'description': 'Surface reflectance',
            'arguments': [__toastring]
        },
        'temp': {
            'description': 'Brightness (apparent) temperature',
            'toa': True
        },
        'acca': {
            'description': 'Automated Cloud Cover Assessment',
            'arguments': [
                'X: erosion kernel diameter in pixels (default: 5)',
                'Y: dilation kernel diameter in pixels (default: 10)',
                'Z: cloud height in meters (default: 4000)'
            ],
            'nargs': '*',
            'toa': True
        },
        'fmask': {
            'description': 'Fmask cloud cover',
            'nargs': '*',
            'toa': True
        },
        'tcap': {
            'description': 'Tassled cap transformation',
            'toa': True
        },
        'dn': {
            'description': 'Raw digital numbers',
            'toa': True
        },
        'volref': {
            'description': 'Volumetric water reflectance - valid for water only',
            'arguments': [__toastring]
        },
        'wtemp': {
            'description': 'Water temperature (atmospherically correct) - valid for water only',
            # It's not really TOA, but the product code will take care of atm correction itself
            'toa': True
        },
        #'Indices': {
        'bi': {
            'description': 'Brightness Index',
            'arguments': [__toastring]
        },
        'evi': {
            'description': 'Enhanced Vegetation Index',
            'arguments': [__toastring]
        },
        'lswi': {
            'description': 'Land Surface Water Index',
            'arguments': [__toastring]
        },
        'msavi2': {
            'description': 'Modified Soil-Adjusted Vegetation Index (revised)',
            'arguments': [__toastring]
        },
        'ndsi': {
            'description': 'Normalized Difference Snow Index',
            'arguments': [__toastring]
        },
        'ndvi': {
            'description': 'Normalized Difference Vegetation Index',
            'arguments': [__toastring]
        },
        'ndwi': {
            'description': 'Normalized Difference Water Index',
            'arguments': [__toastring]
        },
        'satvi': {
            'description': 'Soil-Adjusted Total Vegetation Index',
            'arguments': [__toastring]
        },
        #'Tillage Indices': {
        'ndti': {
            'description': 'Normalized Difference Tillage Index',
            'arguments': [__toastring]
        },
        'crc': {
            'description': 'Crop Residue Cover',
            'arguments': [__toastring]
        },
        'sti': {
            'description': 'Standard Tillage Index',
            'arguments': [__toastring]
        },
        'isti': {
            'description': 'Inverse Standard Tillage Index',
            'arguments': [__toastring]
        },
    }

    def process(self, products=None, overwrite=False, **kwargs):
        """ Make sure all products have been processed """
        products = super(LandsatData, self).process(products, overwrite, **kwargs)
        if len(products) == 0:
            return

        start = datetime.now()

        # Add the sensor for this date to the basename
        self.basename = self.basename + '_' + self.sensor_set[0]

        # Read the assets
        try:
            img = self._readraw()
        except Exception, e:
            VerboseOut(traceback.format_exc(), 5)
            raise Exception('Error reading %s: %s' % (basename(self.assets[''].filename), e))

        meta = self.assets[''].meta
        visbands = self.assets[''].visbands
        lwbands = self.assets[''].lwbands
        md = self.meta_dict()

        # running atmosphere if any products require it
        toa = True
        for val in products.requested.values():
            toa = toa and (self._products[val[0]].get('toa', False) or 'toa' in val)
        if not toa:
            start = datetime.now()
            if not settings().REPOS[self.Repository.name]['6S']:
                raise Exception('6S is required for atmospheric correction')
            try:
                wvlens = [(meta[b]['wvlen1'], meta[b]['wvlen2']) for b in visbands]
                geo = self.metadata['geometry']
                atm6s = SIXS(visbands, wvlens, geo, self.metadata['datetime'], sensor=self.sensor_set[0])
                md["AOD Source"] = str(atm6s.aod[0])
                md["AOD Value"] = str(atm6s.aod[1])
            except Exception, e:
                VerboseOut(traceback.format_exc(), 4)
                raise Exception('Problem running 6S atmospheric model: %s' % e)

        # Break down by group
        groups = products.groups()

        # create non-atmospherically corrected apparent reflectance and temperature image
        reflimg = gippy.GeoImage(img)
        theta = numpy.pi * self.metadata['geometry']['solarzenith'] / 180.0
        sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (float(self.day) - 4.0) / 180.0))
        for col in self.assets[''].visbands:
            reflimg[col] = img[col] * (1.0 / ((meta[col]['E'] * numpy.cos(theta)) / (numpy.pi * sundist * sundist)))
        for col in self.assets[''].lwbands:
            reflimg[col] = (((img[col].pow(-1)) * meta[col]['K1'] + 1).log().pow(-1)) * meta[col]['K2'] - 273.15

        # This is landsat, so always just one sensor for a given date
        sensor = self.sensors['']

        # Process standard products
        for key, val in groups['Standard'].items():
            start = datetime.now()
            # TODO - update if no atmos desired for others
            toa = self._products[val[0]].get('toa', False) or 'toa' in val
            # Create product
            try:
                fname = os.path.join(self.path, self.basename + '_' + key)
                if val[0] == 'acca':
                    s_azim = self.metadata['geometry']['solarazimuth']
                    s_elev = 90 - self.metadata['geometry']['solarzenith']
                    try:
                        erosion = int(val[1]) if len(val) > 1 else 5
                        dilation = int(val[2]) if len(val) > 2 else 10
                        cloudheight = int(val[3]) if len(val) > 3 else 4000
                    except:
                        erosion = 5
                        dilation = 10
                        cloudheight = 4000
                    imgout = ACCA(reflimg, fname, s_elev, s_azim, erosion, dilation, cloudheight)
                elif val[0] == 'fmask':
                    try:
                        tolerance = int(val[1]) if len(val) > 1 else 3
                        dilation = int(val[2]) if len(val) > 2 else 5
                    except:
                        tolerance = 3
                        dilation = 5
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
                fname = imgout.Filename()
                imgout.SetMeta(md)
                imgout = None
                self.AddFile(sensor, key, fname)
                VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname), datetime.now() - start), 1)
            except Exception, e:
                VerboseOut('Error creating product %s for %s: %s' % (key, basename(self.assets[''].filename), e), 2)
                VerboseOut(traceback.format_exc(), 3)

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
                    img[col] = ((img[col] - atm6s.results[col][1]) / atm6s.results[col][0]) * (1.0 / atm6s.results[col][2])
                prodout = Indices(img, dict(zip([p[0] for p in indices.values()], fnames)), md)
                prodout = dict(zip(indices.keys(), prodout.values()))
                [self.AddFile(sensor, key, fname) for key, fname in prodout.items()]
            VerboseOut(' -> %s: processed %s in %s' % (self.basename, indices0.keys(), datetime.now() - start), 1)
        img = None

        # cleanup directory
        try:
            if settings().REPOS[self.Repository.name]['extract']:
                for bname in self.assets[''].datafiles():
                    if bname[-7:] != 'MTL.txt':
                        files = glob.glob(os.path.join(self.path, bname) + '*')
                        RemoveFiles(files)
            shutil.rmtree(os.path.join(self.path, 'modtran'))
        except:
            #VerboseOut(traceback.format_exc(), 4)
            pass

    def filter(self, pclouds=100, **kwargs):
        """ Check if tile passes filter """
        if pclouds < 100:
            self.meta()
            if self.metadata['clouds'] > pclouds:
                return False
        return True

    def meta(self):
        """ Read in Landsat MTL (metadata) file """

        # test if metadata already read in, if so, return

        datafiles = self.assets[''].datafiles()
        mtlfilename = [f for f in datafiles if 'MTL.txt' in f][0]
        if not os.path.exists(mtlfilename):
            mtlfilename = self.assets[''].extract([mtlfilename])[0]
        # Read MTL file
        try:
            text = open(mtlfilename, 'r').read()
        except IOError as e:
            raise Exception('({})'.format(e))

        smeta = self.assets['']._sensors[self.sensor_set[0]]

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

        self.metadata = {
            'filenames': filenames,
            'gain': gain,
            'offset': offset,
            'dynrange': dynrange,
            'geometry': _geometry,
            'datetime': dt,
            'clouds': clouds
        }
        #self.metadata.update(smeta)

    @classmethod
    def meta_dict(cls):
        meta = super(LandsatData, cls).meta_dict()
        meta['GIPS-landsat Version'] = cls.version
        return meta

    def _readraw(self):
        """ Read in Landsat bands using original tar.gz file """
        start = datetime.now()
        # make sure metadata is loaded
        self.meta()

        if settings().REPOS[self.Repository.name]['extract']:
            # Extract all files
            datafiles = self.assets[''].extract(self.metadata['filenames'])
        else:
            # Use tar.gz directly using GDAL's virtual filesystem
            datafiles = [os.path.join('/vsitar/' + self.assets[''].filename, f) for f in self.metadata['filenames']]

        image = gippy.GeoImage(datafiles)
        image.SetNoData(0)

        # TODO - set appropriate metadata
        #for key,val in meta.iteritems():
        #    image.SetMeta(key,str(val))

        # Geometry used for calculating incident irradiance
        colors = self.assets['']._sensors[self.sensor_set[0]]['colors']
        for bi in range(0, len(self.metadata['filenames'])):
            image.SetBandName(colors[bi], bi + 1)
            # need to do this or can we index correctly?
            band = image[bi]
            band.SetGain(self.metadata['gain'][bi])
            band.SetOffset(self.metadata['offset'][bi])
            dynrange = self.metadata['dynrange'][bi]
            band.SetDynamicRange(dynrange[0], dynrange[1])
            image[bi] = band

        VerboseOut('%s: read in %s' % (image.Basename(), datetime.now() - start), 2)
        return image

    @classmethod
    def extra_arguments(cls):
        return {
            '--%clouds': {
                'dest': 'pclouds',
                'help': 'Threshold of max %% cloud cover',
                'default': 100,
                'type': int
            },
        }
