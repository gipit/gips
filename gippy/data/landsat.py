#!/usr/bin/env python

import os, sys, errno
import argparse
import glob
import re

from datetime import datetime
import shutil
import numpy

from copy import deepcopy
import traceback
from collections import OrderedDict

import gippy
from gippy.atmosphere import atmosphere

from gippy.data.core import Data
from gippy.utils import VerboseOut, File2List



class LandsatData(Data):
    """ Represents a single date and temporal extent along with (existing) product variations (raw, ref, toaref, ind, ndvi, etc) """
    name = 'Landsat'
    sensors = {'LT4': 'Landsat 4', 'LT5': 'Landsat 5', 'LE7': 'Landsat 7', 'LC8': 'Landsat 8'}
    _rootdir = '/titan/data/landsat/tiles'
    _tiles_vector = 'landsat_wrs'
    _tiles_attribute = 'pr'
    _assetpattern = 'L*.tar.gz'
    #_pattern = r'^L[TEC][4578].*\.tar\.gz$'
    _prodpattern = '*.tif'
    _metapattern = 'MTL.txt'
    _defaultresolution = [30.0, 30.0]

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
        ('toarad',{
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

    @classmethod
    def inspect(cls, filename):
        """ Inspect a single file and get some metadata """
        path,basename = os.path.split(filename)
        tile = basename[3:9]
        year = basename[9:13]
        doy = basename[13:16]
        return {
            'filename': filename,
            'datafiles': [],    # not needed for inspect, populate in meta()
            'tile': tile,
            'date': datetime.strptime(year+doy,"%Y%j"),
            'basename': basename[:-12],
            #'path':os.path.join(cls._rootdir,tile,year+doy),
            'sensor': basename[0:3],
            'products': {}
        }

    def _readmeta(self, tile):
        """ Read in Landsat MTL (metadata) file """
        tdata = self.tiles[tile]
        mtlfilename = self.extracthdr(tdata['assets'][0])

        VerboseOut('reading %s' % mtlfilename, 3)
        # Read MTL file
        try:
            text = open(mtlfilename,'r').read()
        except IOError as e:
            raise Exception('({})'.format(e))

        # Find out what Landsat this is
        try:
            id = int(os.path.basename(mtlfilename)[1:2])
        except:
            id = int(os.path.basename(mtlfilename)[2:3])

        # Set sensor specific constants
        if id == 5:
            bands = ['1','2','3','4','5','6','7']
            colors = ["BLUE","GREEN","RED","NIR","SWIR1","LWIR","SWIR2"]
            # TODO - update bands with actual L5 values (these are L7)
            bandlocs = [0.4825, 0.565, 0.66, 0.825, 1.65, 11.45, 2.22]
            bandwidths = [0.065, 0.08, 0.06, 0.15, 0.2, 2.1, 0.26]
            E = [1983, 1796, 1536, 1031, 220.0, 0, 83.44]
            K1 = [0, 0, 0, 0, 0, 607.76, 0]
            K2 = [0, 0, 0, 0, 0, 1260.56, 0]
            oldbands = bands
        elif id == 7:
            #bands = ['1','2','3','4','5','6_VCID_1','6_VCID_2','7','8']
            bands = ['1','2','3','4','5','6_VCID_1','7']
            colors = ["BLUE","GREEN","RED","NIR","SWIR1","LWIR","SWIR2"]
            bandlocs = [0.4825, 0.565, 0.66, 0.825, 1.65, 11.45, 2.22]
            bandwidths = [0.065, 0.08, 0.06, 0.15, 0.2, 2.1, 0.26]
            E = [1997, 1812, 1533, 1039, 230.8, 0, 84.90]
            K1 = [0, 0, 0, 0, 0, 666.09, 0]
            K2 = [0, 0, 0, 0, 0, 1282.71, 0]
            oldbands = deepcopy(bands)
            oldbands[5] = '61'
            oldbands[6] = '62'
        elif id == 8:
            bands = ['1','2','3','4','5','6','7','9'] #,'10','11']
            colors = ["COASTAL","BLUE","GREEN","RED","NIR","SWIR1","SWIR2","CIRRUS"] #,"LWIR1","LWIR2"]
            bandlocs = [0.443, 0.4825, 0.5625, 0.655, 0.865, 1.610, 2.2, 1.375] #, 10.8, 12.0]
            bandwidths = [0.01, 0.0325, 0.0375, 0.025, 0.02, 0.05, 0.1, 0.015] #, 0.5, 0.5]
            E = [2638.35, 2031.08, 1821.09, 2075.48, 1272.96, 246.94, 90.61, 369.36] #, 0, 0]
            K1 = [0, 0, 0, 0, 0, 0, 0, 0] #774.89, 480.89]
            K2 = [0, 0, 0, 0, 0, 0, 0, 0] #1321.08, 1201.14]
            oldbands = bands
        else:
            raise Exception('Landsat%s? not recognized' % id)

        # Process MTL text - replace old metadata tags with new
        # NOTE This is not comprehensive, there may be others
        text = text.replace('ACQUISITION_DATE','DATE_ACQUIRED')
        text = text.replace('SCENE_CENTER_SCAN_TIME','SCENE_CENTER_TIME')
        for (ob,nb) in zip(oldbands,bands):
            text = re.sub(r'\WLMIN_BAND'+ob,'RADIANCE_MINIMUM_BAND_'+nb,text)
            text = re.sub(r'\WLMAX_BAND'+ob,'RADIANCE_MAXIMUM_BAND_'+nb,text)
            text = re.sub(r'\WQCALMIN_BAND'+ob,'QUANTIZE_CAL_MIN_BAND_'+nb,text)
            text = re.sub(r'\WQCALMAX_BAND'+ob,'QUANTIZE_CAL_MAX_BAND_'+nb,text)
            text = re.sub(r'\WBAND'+ob+'_FILE_NAME','FILE_NAME_BAND_'+nb,text)
        for l in ('LAT','LON','MAPX','MAPY'):
            for c in ('UL','UR','LL','LR'):
                text = text.replace('PRODUCT_'+c+'_CORNER_'+l, 'CORNER_'+c+'_'+l+'_PRODUCT')
        text = text.replace('\x00','')
        # Remove junk
        lines = text.split('\n')
        mtl = dict()
        for l in lines:
            meta = l.replace('\"',"").strip().split('=')
            if len(meta) > 1:
                key = meta[0].strip()
                item = meta[1].strip()
                if key != "GROUP" and key !="END_GROUP": mtl[key] = item

        # Extract useful metadata
        lats = ( float(mtl['CORNER_UL_LAT_PRODUCT']), float(mtl['CORNER_UR_LAT_PRODUCT']),
                float(mtl['CORNER_LL_LAT_PRODUCT']), float(mtl['CORNER_LR_LAT_PRODUCT']))
        lons = ( float(mtl['CORNER_UL_LON_PRODUCT']), float(mtl['CORNER_UR_LON_PRODUCT']),
            float(mtl['CORNER_LL_LON_PRODUCT']), float(mtl['CORNER_LR_LON_PRODUCT']))
        lat = (min(lats) + max(lats))/2.0
        lon = (min(lons) + max(lons))/2.0
        dt = datetime.strptime(mtl['DATE_ACQUIRED'] + ' ' +
            mtl['SCENE_CENTER_TIME'][:-2],'%Y-%m-%d %H:%M:%S.%f')
        seconds = (dt.second + dt.microsecond/1000000.)/3600.
        dectime = dt.hour + dt.minute/60.0 + seconds
        try:
            clouds = float(mtl['CLOUD_COVER'])
        except: clouds = 0

        # Band metadata
        bandmeta = []
        filenames = []
        for i,b in enumerate(bands):
            band = {
                'num': i+1,
                'id': b,
                'color': colors[i],
                'E': E[i],
                'wavelength': bandlocs[i],
                'bandwidth': bandwidths[i],
                'minval': int(float(mtl['QUANTIZE_CAL_MIN_BAND_'+b])),
                'maxval': int(float(mtl['QUANTIZE_CAL_MAX_BAND_'+b])),
                'minrad': float(mtl['RADIANCE_MINIMUM_BAND_'+b]),
                'maxrad': float(mtl['RADIANCE_MAXIMUM_BAND_'+b]),
                'k1': K1[i],
                'k2': K2[i]
            }
            bandmeta.append(band)
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
            'JulianDay': (dt - datetime(dt.year,1,1)).days + 1,
            'DecimalTime': dectime,
        }

        # TODO - now that metadata part of LandsatData object some of these keys not needed
        tdata['metadata'] = {
            'sensor': 'Landsat'+str(id),
            'metafilename': mtlfilename,
            'filenames': filenames,
            'geometry': _geometry,
            'datetime': _datetime,
            'bands': bandmeta,
            'clouds': clouds
        }
        return self.tiles[tile]['metadata']

    def processtile(self,tile,products):
        """ Make sure all products have been pre-processed """
        start = datetime.now()
        tdata = self.tiles[tile]
        bname = os.path.basename(tdata['assets'][0])
        try:
            img = self._readraw(tile)
        except Exception,e:
            VerboseOut('Error reading %s %s' % (bname,e), 2)
            VerboseOut(traceback.format_exc(), 4)

        #if self._products[p]['atmcorr']:
        # running atmosphere automatically, for now
        atmospheres = [atmosphere(i,tdata['metadata']) for i in range(1,img.NumBands()+1)]
        VerboseOut('%s: read in %s' % (bname,datetime.now() - start)) 

        for p,fout in products.items():
            try: 
                start = datetime.now()
                if (self._products[p]['atmcorr']):
                    for i in range(0,img.NumBands()):
                        b = img[i]
                        b.SetAtmosphere(atmospheres[i])
                        img[i] = b
                        VerboseOut('atmospherically correcting',3)
                else: img.ClearAtmosphere()

                try:
                    fcall = 'gippy.%s(img, fout)' % self._products[p]['function']
                    VerboseOut(fcall,2)
                    eval(fcall)
                    #if overviews: imgout.AddOverviews()
                    fname = glob.glob(fout+'*')[0]
                    tdata['products'][p] = fname
                    VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname),datetime.now()-start))
                except Exception,e:
                    VerboseOut('Error creating product %s for %s: %s' % (p,bname,e),3)
                    VerboseOut(traceback.format_exc(),4)
            # double exception? is this necessary ?
            except Exception,e:
                VerboseOut('Error processing %s: %s' % (info['basename'],e),2)
                VerboseOut(traceback.format_exc(), 4)

            img = None
            # cleanup directory
            try:
            #    for bname self.tiles[tile]['datafiles']:
            #        files = glob.glob(os.path.join(data['path'],bname)+'*')
            #            for f in files: os.remove(f)
                shutil.rmtree(os.path.join(tdata['path'],'modtran'))
            except: pass

    def filter(self, tile, maxclouds=100):
        """ Check if tile passes filter """
        if maxclouds < 100:
            meta = cls._readmeta(tile)
            if meta['clouds'] > maxclouds:
                return False
        return True

    @classmethod
    def feature2tile(cls,feature):
        tile = super(LandsatData, cls).feature2tile(feature)
        return tile.zfill(6)

    def project(self, res=[30,30], datadir='data_landsat'):
        """ Create reprojected, mosaiced images for this site and date """
        if res is None: res=[30,30]
        super(LandsatData, self).project(res=res, datadir=datadir)

    def _readraw(self,tile,bandnums=[]):
        """ Read in Landsat bands using original tar.gz file """
        # Read in collection metadata
        self._readmeta(tile)

        tiledata = self.tiles[tile]

        if len(bandnums) != 0:
            bandnums = numpy.array(bandnums)
        else:
            bandnums = numpy.arange(0,len(tiledata['metadata']['bands'])) + 1

        # Extract desired files from tarfile
        index = self.extractdata(tiledata['assets'][0])

        filenames = []
        for b in bandnums:
            fname = tiledata['metadata']['filenames'][b-1]
            filenames.append( os.path.join(self.path(tile,self.date),fname) )

        if gippy.Options.Verbose() > 2:
            print 'Landsat files: '
            for f in filenames: print ' %s' % f

        # Read in files
        image = gippy.GeoImage(filenames[0])
        del filenames[0]
        for f in filenames: image.AddBand(gippy.GeoImage(f)[0])

        # Geometry used for calculating incident irradiance
        theta = numpy.pi * tiledata['metadata']['geometry']['solarzenith']/180.0
        sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (tiledata['metadata']['datetime']['JulianDay']-4.0)/180.0))

        # Set metadata
        image.SetNoData(0)
        image.SetUnits('radiance')
        # TODO - set appropriate metadata
        #for key,val in meta.iteritems():
        #    image.SetMeta(key,str(val))

        for b in bandnums:
            bandmeta = tiledata['metadata']['bands'][b-1]
            gain = (bandmeta['maxrad']-bandmeta['minrad']) / (bandmeta['maxval']-bandmeta['minval'])
            offset = bandmeta['minrad']
            i = bandmeta['num']-1
            band = image[i]
            band.SetGain(gain)
            band.SetOffset(offset)
            band.SetDynamicRange(bandmeta['minval'],bandmeta['maxval'])
            band.SetEsun( (bandmeta['E'] * numpy.cos(theta)) / (numpy.pi * sundist * sundist) )
            band.SetThermal(bandmeta['k1'],bandmeta['k2'])
            image[i] = band
            image.SetColor(bandmeta['color'],i+1)
        return image

    #@classmethod
    #def add_subparsers(cls, parser):
    #    p = parser.add_parser('cloudmask', help='Create cloud masks', parents=[cls.args_inventory()],
    #        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

def main(): LandsatData.main()
