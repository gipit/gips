#!/usr/bin/env python
################################################################################
#    GIPPY: Geospatial Image Processing library for Python
#
#    Copyright (C) 2014 Matthew A Hanson
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
import glob
import re
from datetime import datetime
import shutil
import numpy
from collections import OrderedDict

import gippy
from gippy.atmosphere import atmosphere
from gippy.data.core import Asset, Tile, Data
from gippy.utils import VerboseOut

import traceback
from pdb import set_trace


class LandsatAsset(Asset):
    """ Landsat asset (original raw tar file) """
    _rootpath = '/titan/data/landsat'

    # combine sensormeta with sensor
    _sensors = {
        'LT4': {
            'description': 'Landsat 4',
        },
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
            'K2': [0, 0, 0, 0, 0, 1260.56, 0]
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
        },
        'LC8': {
            'description': 'Landsat 8',
            'bands': ['1', '2', '3', '4', '5', '6', '7', '9'],  # ,'10','11']
            'oldbands': ['1', '2', '3', '4', '5', '6', '7', '9'],  # ,'10','11']
            'colors': ["COASTAL", "BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2", "CIRRUS"],     # ,"LWIR1","LWIR2"]
            'bandlocs': [0.443, 0.4825, 0.5625, 0.655, 0.865, 1.610, 2.2, 1.375],    # , 10.8, 12.0]
            'bandwidths': [0.01, 0.0325, 0.0375, 0.025, 0.02, 0.05, 0.1, 0.015],        # , 0.5, 0.5]
            'E': [2638.35, 2031.08, 1821.09, 2075.48, 1272.96, 246.94, 90.61, 369.36],   # , 0, 0]
            'K1': [0, 0, 0, 0, 0, 0, 0, 0],  # 774.89, 480.89]
            'K2': [0, 0, 0, 0, 0, 0, 0, 0],  # 1 321.08, 1201.14]
        }
    }

    # TODO - consider assets and sensors relationship ?
    _assets = {
        '': {
            'pattern': 'L*.tar.gz'
            #_pattern = r'^L[TEC][4578].*\.tar\.gz$'
        }
    }

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(LandsatAsset, self).__init__(filename)
        self.basename = self.basename[:-12]
        self.sensor = self.basename[0:3]
        self.tile = self.basename[3:9]
        year = self.basename[9:13]
        doy = self.basename[13:16]
        self.date = datetime.strptime(year+doy, "%Y%j")


class LandsatTile(Tile):

    Asset = LandsatAsset

    _prodpattern = '*.tif'
    _products = OrderedDict([
        ('rgb', {
            'description': 'RGB image for viewing (quick processing)',
            'function': 'RGB',
            'atmcorr': False,
        }),
        ('rad', {
            'description': 'Surface leaving radiance',
            'function': 'Rad',
            'atmcorr': True,
        }),
        ('toarad', {
            'description': 'Top of atmosphere radiance',
            'function': 'Rad',
            'atmcorr': False,
        }),
        ('ref', {
            'description': 'Surface reflectance',
            'function': 'Ref',
            'atmcorr': True,
        }),
        ('toaref', {
            'description': 'Top of atmosphere reflectance',
            'function': 'Ref',
            'atmcorr': False,
        }),
        #('ind', {
        #    'description': 'Atmospherically corrected common indices (NDVI,EVI,LSWI,NDSI,BI)',
        #    'function': 'Indices',
        #    'atmcorr': True,
        #}),
        #('toaind', {
        #    'description': 'Top of atmosphere common indices (NDVI,EVI,LSWI,NDSI,BI)',
        #    'function': 'Indices',
        #    'atmcorr': False,
        #}),
        ('ndvi', {
            'description': 'Atmospherically corrected NDVI',
            'function': 'NDVI',
            'atmcorr': True,
        }),
        ('evi', {
            'description': 'Atmospherically corrected EVI',
            'function': 'EVI',
            'atmcorr': True,
        }),
        ('lswi', {
            'description': 'Atmospherically corrected LSWI',
            'function': 'LSWI',
            'atmcorr': True,
        }),
        ('ndsi', {
            'description': 'Atmospherically corrected NDSI',
            'function': 'NDSI',
            'atmcorr': True,
        }),
        ('ndti', {
            'description': 'Atmospherically corrected NDTI',
            'function': 'NDTI',
            'atmcorr': True,
        }),
        ('satvi', {
            'description': 'Atmospherically corrected SATVI',
            'function': 'SATVI',
            'atmcorr': True,
        }),
        ('toabi', {
            'description': 'TOA BI',
            'function': 'BI',
            'atmcorr': False,
        }),
        ('toandvi', {
            'description': 'TOA NDVI',
            'function': 'NDVI',
            'atmcorr': False,
        }),

        # Tillage
        ('crc', {
            'description': 'Crop Residue Cover',
            'function': 'CRC',
            'atmcorr': True,
        }),
        ('crcm', {
            'description': 'Crop Residue Cover',
            'function': 'CRCm',
            'atmcorr': True,
        }),
        ('isti', {
            'description': 'Inverse Standard Tillage Index',
            'function': 'iSTI',
            'atmcorr': True,
        }),
        ('sti', {
            'description': 'Standard Tillage Index',
            'function': 'STI',
            'atmcorr': True,
        }),

        # Other
        ('acca', {
            'description': 'ACCA Cloud Algorithm',
            'function': 'ACCA',
            'atmcorr': False,
        }),
    ])

    def process(self, products):
        """ Make sure all products have been pre-processed """
        start = datetime.now()
        bname = os.path.basename(self.assets[''].filename)
        #try:
        img = self._readraw()
        #except Exception, e:
            #VerboseOut(traceback.format_exc(), 4)
        #    raise Exception('Error reading %s %s' % (bname, e))

        runatm = False
        for p in products:
            if self._products[p]['atmcorr']:
                runatm = True
        if runatm:
            atmospheres = [atmosphere(i, self.metadata) for i in range(1, img.NumBands()+1)]
        VerboseOut('%s: read in %s' % (bname, datetime.now() - start))

        for p, fout in products.items():
            try:
                start = datetime.now()
                if (self._products[p]['atmcorr']):
                    for i in range(0, img.NumBands()):
                        b = img[i]
                        b.SetAtmosphere(atmospheres[i])
                        img[i] = b
                        VerboseOut('Band %i: atmospherically correcting' % (i+1), 3)
                else:
                    img.ClearAtmosphere()
                try:
                    # If product is ACCA, get params from suffix, if possible
                    mat = re.match('.*acca(?:_(?P<eds>(\d+)_(\d+)_(\d+))|(?:_(?P<ed>(\d+)_(\d+)))|(?:_(?P<e>(\d+))))*', fout)
                    if mat:
                        if mat.group('eds'):
                            erosion = int(mat.group(2))
                            dilation = int(mat.group(3))
                            cloudheight = int(mat.group(4))
                        elif mat.group('ed'):
                            erosion = int(mat.group(6))
                            dilation = int(mat.group(7))
                            cloudheight = 4000
                        elif mat.group('e'):
                            erosion = int(mat.group(9))
                            dilation = 10
                            cloudheight = 4000
                        else:
                            erosion = 5
                            dilation = 10
                            cloudheight = 4000
                        s_azim = self.metadata['geometry']['solarazimuth']
                        s_elev = 90 - self.metadata['geometry']['solarzenith']
                        VerboseOut('s_elev, s_azim, erosion, dilation, cloudheight = ' +
                                   str((s_elev, s_azim, erosion, dilation, cloudheight)))
                        fcall = 'gippy.%s(img, fout, s_elev, s_azim, erosion, dilation, cloudheight)' % self._products[p]['function']
                    else:
                        fcall = 'gippy.%s(img, fout)' % self._products[p]['function']
                    VerboseOut(fcall, 2)
                    eval(fcall)
                    #if overviews: imgout.AddOverviews()
                    fname = glob.glob(fout+'*')[0]
                    self.products[p] = fname
                    VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname), datetime.now()-start))
                except Exception, e:
                    VerboseOut('Error creating product %s for %s: %s' % (p, bname, e), 3)
                    VerboseOut(traceback.format_exc(), 4)
            # double exception? is this necessary ?
            except Exception, e:
                VerboseOut('Error processing %s: %s' % (self.basename, e), 2)
                VerboseOut(traceback.format_exc(), 4)
        img = None
        # cleanup directory
        try:
            for bname in self.assets[''].datafiles:
                files = glob.glob(os.path.join(self.path, bname)+'*')
                print 'remove: ', files
                RemoveFiles(files)
            shutil.rmtree(os.path.join(self.path, 'modtran'))
        except:
            pass

    def filter(self, maxclouds=100):
        """ Check if tile passes filter """
        if maxclouds < 100:
            # shouldnt have to read meta again
            meta = cls.meta(tile)
            if meta['clouds'] > maxclouds:
                return False
        return True

    def meta(self):
        """ Read in Landsat MTL (metadata) file """

        # test if metadata already read in, if so, return

        datafiles = self.assets[''].datafiles()
        mtlfilename = self.assets[''].extract(([f for f in datafiles if 'MTL.txt' in f]))[0]

        VerboseOut('reading %s' % mtlfilename, 3)
        # Read MTL file
        try:
            text = open(mtlfilename, 'r').read()
        except IOError as e:
            raise Exception('({})'.format(e))

        sensor_meta = self.assets[''].sensor_meta()

        # Process MTL text - replace old metadata tags with new
        # NOTE This is not comprehensive, there may be others
        text = text.replace('ACQUISITION_DATE', 'DATE_ACQUIRED')
        text = text.replace('SCENE_CENTER_SCAN_TIME', 'SCENE_CENTER_TIME')
        for (ob, nb) in zip(sensor_meta['oldbands'], sensor_meta['bands']):
            text = re.sub(r'\WLMIN_BAND'+ob, 'RADIANCE_MINIMUM_BAND_'+nb, text)
            text = re.sub(r'\WLMAX_BAND'+ob, 'RADIANCE_MAXIMUM_BAND_'+nb, text)
            text = re.sub(r'\WQCALMIN_BAND'+ob, 'QUANTIZE_CAL_MIN_BAND_'+nb, text)
            text = re.sub(r'\WQCALMAX_BAND'+ob, 'QUANTIZE_CAL_MAX_BAND_'+nb, text)
            text = re.sub(r'\WBAND'+ob+'_FILE_NAME', 'FILE_NAME_BAND_'+nb, text)
        for l in ('LAT', 'LON', 'MAPX', 'MAPY'):
            for c in ('UL', 'UR', 'LL', 'LR'):
                text = text.replace('PRODUCT_'+c+'_CORNER_'+l, 'CORNER_'+c+'_'+l+'_PRODUCT')
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
        lat = (min(lats) + max(lats))/2.0
        lon = (min(lons) + max(lons))/2.0
        dt = datetime.strptime(mtl['DATE_ACQUIRED'] + ' ' + mtl['SCENE_CENTER_TIME'][:-2], '%Y-%m-%d %H:%M:%S.%f')
        seconds = (dt.second + dt.microsecond/1000000.)/3600.
        dectime = dt.hour + dt.minute/60.0 + seconds
        try:
            clouds = float(mtl['CLOUD_COVER'])
        except:
            clouds = 0

        filenames = []
        gain = []
        offset = []
        dynrange = []
        for i, b in enumerate(sensor_meta['bands']):
            minval = int(float(mtl['QUANTIZE_CAL_MIN_BAND_'+b]))
            maxval = int(float(mtl['QUANTIZE_CAL_MAX_BAND_'+b]))
            minrad = float(mtl['RADIANCE_MINIMUM_BAND_'+b])
            maxrad = float(mtl['RADIANCE_MAXIMUM_BAND_'+b])
            gain.append((maxrad-minrad)/(maxval-minval))
            offset.append(minrad)
            dynrange.append((minval, maxval))
            filenames.append(mtl['FILE_NAME_BAND_'+b].strip('\"'))

        _geometry = {
            'solarzenith': (90.0 - float(mtl['SUN_ELEVATION'])),
            'solarazimuth': float(mtl['SUN_AZIMUTH']),
            'zenith': 0.0,
            'azimuth': 180.0,
            'lat': lat,
            'lon': lon,
        }

        _datetime = {
            'datetime': dt,
            'JulianDay': (dt - datetime(dt.year, 1, 1)).days + 1,
            'DecimalTime': dectime,
        }

        # TODO - now that metadata part of LandsatData object some of these keys not needed
        self.metadata = {
            'filenames': filenames,
            'path': self.path,
            'gain': gain,
            'offset': offset,
            'dynrange': dynrange,
            'geometry': _geometry,
            'datetime': _datetime,
            'clouds': clouds
        }
        self.metadata.update(sensor_meta)

    def _readraw(self):  # , bandnums=[]):
        """ Read in Landsat bands using original tar.gz file """
        # make sure metadata is loaded
        self.meta()

        # TODO - allow for band numbers
        #if len(bandnums) != 0:
        #    bandnums = numpy.array(bandnums)
        #else:
        #    bandnums = numpy.arange(0, len(self.metadata['bands'])) + 1

        # Extract all files
        datafiles = self.assets[''].extract(self.metadata['filenames'])
        VerboseOut('Extracted Landsat files:', 3)
        VerboseOut(datafiles, 3)

        # TODO - replace with single call (add GeoImage(vector<string> to GIP)
        image = gippy.GeoImage(datafiles[0])
        del datafiles[0]
        for f in datafiles:
            image.AddBand(gippy.GeoImage(f)[0])
        # Set metadata
        image.SetNoData(0)
        image.SetUnits('radiance')

        # TODO - set appropriate metadata
        #for key,val in meta.iteritems():
        #    image.SetMeta(key,str(val))

        # TODO - most setting here removed when GIP updated with late binding
        # Geometry used for calculating incident irradiance
        theta = numpy.pi * self.metadata['geometry']['solarzenith']/180.0
        sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (self.metadata['datetime']['JulianDay']-4.0)/180.0))
        for bi in range(0, len(self.metadata['filenames'])):
            image.SetColor(self.metadata['colors'][bi], bi+1)
            # need to do this or can we index correctly?
            band = image[bi]
            band.SetGain(self.metadata['gain'][bi])
            band.SetOffset(self.metadata['offset'][bi])
            dynrange = self.metadata['dynrange'][bi]
            band.SetDynamicRange(dynrange[0], dynrange[1])
            band.SetEsun((self.metadata['E'][bi] * numpy.cos(theta)) / (numpy.pi * sundist * sundist))
            band.SetThermal(self.metadata['K1'][bi], self.metadata['K2'][bi])
            image[bi] = band

        return image


class LandsatData(Data):
    """ Single date and temporal extent with all assets and products (ref, toaref, ind, ndvi, etc) """
    name = 'Landsat'

    _defaultresolution = [30.0, 30.0]

    Tile = LandsatTile

    _tiles_vector = 'landsat_wrs'
    _vectoratt = 'pr'

    @classmethod
    def feature2tile(cls, feature):
        tile = super(LandsatData, cls).feature2tile(feature)
        return tile.zfill(6)

    def project(self, res=[30, 30], datadir='data_landsat'):
        """ Create reprojected, mosaiced images for this site and date """
        if res is None:
            res = [30, 30]
        super(LandsatData, self).project(res=res, datadir=datadir)


def main():
    LandsatData.main()
