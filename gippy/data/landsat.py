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
from gippy.data.core import Repository, Asset, Data
from gippy.data.inventory import DataInventory
from gippy.utils import VerboseOut, RemoveFiles

import traceback
from pdb import set_trace


class LandsatRepository(Repository):
    """ Singleton (all class methods) to be overridden by child data classes """
    _rootpath = '/titan/data/landsat'
    _tiles_vector = 'landsat_wrs'
    #_tilesdir = 'tiles.dev'

    # attribute (column) in tiles vector giving tile id
    _tile_attribute = 'pr'

    @classmethod
    def feature2tile(cls, feature):
        tile = super(LandsatRepository, cls).feature2tile(feature)
        return tile.zfill(6)


class LandsatAsset(Asset):
    """ Landsat asset (original raw tar file) """
    Repository = LandsatRepository

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

    _defaultresolution = [30.0, 30.0]

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(LandsatAsset, self).__init__(filename)
        self.basename = self.basename[:-12]
        self.sensor = self.basename[0:3]
        self.tile = self.basename[3:9]
        year = self.basename[9:13]
        doy = self.basename[13:16]
        self.date = datetime.strptime(year+doy, "%Y%j")


class LandsatData(Data):
    name = 'Landsat'

    Asset = LandsatAsset

    _prodpattern = '*.tif'
    _products = {
        #'Standard': {
        # 'rgb': 'RGB image for viewing (quick processing)',
        'rad':  {'description': 'Surface-leaving radiance'}, # , 'args': '?'},
        'ref':  {'description': 'Surface reflectance'}, #, 'args': '?'},
        'temp': {'description': 'Apparent temperature', 'toa': True},
        'acca': {'description': 'Automated Cloud Cover Assesment', 'args': '*', 'toa': True},
        #'Indices': {
        'bi':   {'description': 'Brightness Index'},
        'ndvi': {'description': 'Normalized Difference Vegetation Index'},
        'evi':  {'description': 'Enhanced Vegetation Index'},
        'lswi': {'description': 'Land Surface Water Index'},
        'ndsi': {'description': 'Normalized Difference Snow Index'},
        'satvi': {'description': 'Soil-adjusted Total Vegetation Index'},
        #'Tillage Indices': {
        'ndti': {'description': 'Normalized Difference Tillage Index'},
        'crc':  {'description': 'Crop Residue Cover'},
        'sti':  {'description': 'Standard Tillage Index'},
        'isti': {'description': 'Inverse Standard Tillage Index'},
    }
    _groups = {
        'Standard': ['rad', 'ref', 'temp', 'acca'],
        'Index': ['bi', 'ndvi', 'evi', 'lswi', 'ndsi', 'satvi'],
        'Tillage': ['ndti', 'crc', 'sti', 'isti']
    }
    _defaultproduct = 'ref'

    def process(self, products):
        """ Make sure all products have been processed """
        start = datetime.now()
        bname = os.path.basename(self.assets[''].filename)
        img = self._readraw()

        # running atmosphere
        toa = True
        for val in products.values():
            toa = toa and (self._products[val[0]].get('toa', False) or 'toa' in val)
        if not toa:
            start = datetime.now()
            atmospheres = {}
            for i in range(0,img.NumBands()):
                atmospheres[img[i].Description()] = atmosphere(i+1, self.metadata)
            VerboseOut('Ran atmospheric model in %s' % str(datetime.now()-start), 3)

        # Break down by group
        groups = self.products2groups(products)
        smeta = self.assets[''].sensor_meta()

        visbands = [b for b in smeta['colors'] if b[0:4] != "LWIR"]
        lwbands = [b for b in smeta['colors'] if b[0:4] == "LWIR"]

        # create non-atmospherically corrected apparent reflectance and temperature image
        reflimg = img
        theta = numpy.pi * self.metadata['geometry']['solarzenith']/180.0
        sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (self.metadata['datetime']['JulianDay']-4.0)/180.0))
        Esuns = dict(zip(smeta['colors'],smeta['E']))
        for b in visbands:
            reflimg[b] = img[b] * (1.0/((Esuns[b] * numpy.cos(theta)) / (numpy.pi * sundist * sundist)))
        for b in lwbands:
            reflimg[b] = (((img[b]^(-1))*smeta['K1'][5]+1).log()^(-1))*smeta['K2'][5] - 273.15;

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
                    erosion = int(val[1]) if len(val) > 1 else 5
                    dilation = int(val[2]) if len(val) > 2 else 10
                    cloudheight = int(val[3]) if len(val) > 3 else 4000
                    imgout = gippy.ACCA(reflimg, fname, s_azim, s_elev, erosion, dilation, cloudheight)
                elif val[0] == 'rad':
                    imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(visbands))
                    for i in range(0,imgout.NumBands()):
                        imgout.SetColor(visbands[i],i+1)
                    imgout.SetNoData(-32768)
                    imgout.SetGain(0.1)
                    if toa:
                        for b in visbands:
                            imgout[b].Process(img[b])
                    else:
                        for b in visbands:
                            imgout[b].Process( ((img[b]-atmospheres[b][1])/atmospheres[b][0]) )
                elif val[0] == 'ref':
                    imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(visbands))
                    for i in range(0,imgout.NumBands()):
                        imgout.SetColor(visbands[i],i+1)
                    imgout.SetNoData(-32768)
                    imgout.SetGain(0.0001)
                    if toa:
                        for b in visbands:
                            imgout[b].Process(reflimg[b])
                    else:
                        for b in visbands:
                            imgout[b].Process(((img[b]-atmospheres[b][1])/atmospheres[b][0]) * (1.0/atmospheres[b][2]))
                elif val[0] == 'temp':
                    lwbands = [b for b in smeta['colors'] if b[0:4] == "LWIR"]
                    imgout = gippy.GeoImage(fname, img, gippy.GDT_Int16, len(lwbands))
                    imgout.SetNoData(-32768)
                    imgout.SetGain(0.1)
                    e = 0.95
                    for b in lwbands:
                        if toa:
                            band = img[b]
                        else:
                            band = (img[b] - (atmospheres[b][1] + (1-e) * atmospheres[b][2])) / (atmospheres[b][0] * e)
                        band = (((band^(-1))*smeta['K1'][5]+1).log()^(-1))*smeta['K2'][5] - 273.15;
                        imgout[b].Process(band)
                fname = imgout.Filename()
                imgout = None
                self.products[key] = fname
                VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname), datetime.now()-start))
            except Exception, e:
                VerboseOut('Error creating product %s for %s: %s' % (key, bname, e), 3)
                VerboseOut(traceback.format_exc(), 4)

        # Process Indices
        indices = dict(groups['Index'], **groups['Tillage'])
        if len(indices) > 0:
            start = datetime.now()
            # TODO - this assumes atmospheric correct - what if mix of atmos and non-atmos?
            for b in visbands:
                img[b] = ((img[b]-atmospheres[b][1])/atmospheres[b][0]) * (1.0/atmospheres[b][2])
            fnames = [os.path.join(self.path, self.basename + '_' + key) for key in indices]
            prodarr = dict(zip([indices[p][0] for p in indices.keys()], fnames))
            prodout = gippy.Indices(img, prodarr)
            self.products.update(prodout)
            VerboseOut(' -> %s: processed %s in %s' % (self.basename, indices.keys(), datetime.now()-start))

        img = None
        # cleanup directory
        try:
            for bname in self.assets[''].datafiles():
                files = glob.glob(os.path.join(self.path, bname)+'*')
                RemoveFiles(files)
            shutil.rmtree(os.path.join(self.path, 'modtran'))
        except:
            #VerboseOut(traceback.format_exc(), 4)
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
        mtlfilename = [f for f in datafiles if 'MTL.txt' in f][0]
        if not os.path.exists(mtlfilename):
            mtlfilename = self.assets[''].extract([mtlfilename])[0]
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
        start = datetime.now()
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
        for bi in range(0, len(self.metadata['filenames'])):
            image.SetColor(self.metadata['colors'][bi], bi+1)
            # need to do this or can we index correctly?
            band = image[bi]
            band.SetGain(self.metadata['gain'][bi])
            band.SetOffset(self.metadata['offset'][bi])
            dynrange = self.metadata['dynrange'][bi]
            band.SetDynamicRange(dynrange[0], dynrange[1])
            image[bi] = band

        VerboseOut('%s: read in %s' % (image.Basename(), datetime.now() - start), 3)
        return image

    @classmethod
    def extra_arguments(cls):
        return {}
        return {'--noatmos': {
                'help': 'No atmospheric correction for any products',
                'default': False, 'action': 'store_true'
                }}


def main():
    DataInventory.main(LandsatData)
