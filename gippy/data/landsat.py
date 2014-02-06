#!/usr/bin/env python

import os, sys, errno
import argparse
import glob
import re

import datetime
import shutil
import numpy
import tarfile
from copy import deepcopy
import traceback
from collections import OrderedDict

import gippy
from gippy.atmosphere import atmosphere

from gippy.data.core import Data, DataInventory, VerboseOut
from gippy.data.core import main as datamain

from pdb import set_trace

class LandsatData(Data):
    """ Represents a single date and temporal extent along with (existing) product variations (raw, ref, toaref, ind, ndvi, etc) """
    name = 'Landsat'
    sensors = {'LT4': 'Landsat 4', 'LT5': 'Landsat 5', 'LE7': 'Landsat 7', 'LC8': 'Landsat 8'}
    rootdir = '/titan/data/landsat/tiles'
    _tiles_vector = 'landsat_wrs'
    _tiles_attribute = 'pr'
    _pattern = 'L*.tar.gz'

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

        # Other
        ('acca', {
            'description': 'ACCA Cloud Algorithm',
            'function': 'ACCA',
            'atmcorr': False,
        }),
    ])

    def find_products(self, tile):
        """ Get filename info for specified tile """
        filename = glob.glob(os.path.join(cls.rootdir, tile, self.date.strftime('%Y%j'), '*.tar.gz'))
        if len(filename) == 0:
            raise Exception('No files found')
        elif len(filename) > 1:
            raise Exception('More than 1 file found for same tile/date')
        filename = filename[0]
        path,basename = os.path.split(filename)
        basename = basename[:-12]
        sensor = basename[0:3]

        self.tiles[tile] { 'path': path, 'basename': basename, 'sensor': sensor }
        prod = {'raw': filename}
        for p in self.products:
            fnames = glob.glob( os.path.join(path,'*_%s*.tif' % p) )
            for f in fnames:
                prodname = os.path.splitext(os.path.split(f)[1][len(basename)+1:])[0]
                prod[ prodname ] = f
        self.tiles[tile]['products'] = prod

    @classmethod
    def filter(cls, tile, filename, maxclouds=100):
        """ Check if tile passes filter """
        # TODO - remove filename parameter
        if maxclouds < 100:
            meta = cls._readmeta(tile,filename=filename)
            if meta['clouds'] > maxclouds:
                return False
        return True
        
    @classmethod
    def archive(cls, path=''):
        """ Move landsat files from current directory to archive location """
        fnames = glob.glob(os.path.join(path,cls.pattern))

        numadded = 0
        for f in fnames:
            pathrow = f[3:9]
            year = f[9:13]
            doy = f[13:16]
            path = os.path.join(cls.rootdir,pathrow,year+doy)

            # Make directory
            try:
                #pass
                os.makedirs(path)
            except OSError as exc: # Python >2.5
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    #print 'Directory already exists'
                    pass
                else:
                    raise Exception('Unable to make product directory %s' % path)

            # Move file
            try:
                newf = os.path.join(path,f)
                if not os.path.exists(newf):
                    # Check for older versions
                    existing_files = glob.glob(os.path.join(path,'*tar.gz'))
                    if len(existing_files) > 0:
                        print 'Other version of %s already exist:' % f
                        for ef in existing_files: print '\t%s' % ef
                    shutil.move(f,newf)
                    #print f, ' -> ',path
                    numadded = numadded + 1
            except shutil.Error as err:
                print err
                print f, ' -> problem archiving file'
                #raise Exception('shutils error %s' % err)
                #if exc.errno == errno.EEXIST:
                #    print f, ' removed, already in archive'
        print '%s files added to archive' % numadded
        if numadded != len(fnames):
            print '%s files not added to archive' % (len(fnames)-numadded)

    def process(self, overwrite=False, suffix=''): # , overviews=False):
        """ Make sure all products exist for all tiles, process if necessary """
        for tile, data in self.tiles.items():
            if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
            fout_base = os.path.join(data['path'], data['basename'] + '_')
            fouts = {}
            runatm = False
            # generate all product names and check if already existing
            for product in self.products:
                if self._products[product]['atmcorr']: runatm = True
                if product not in self._products.keys():
                    raise Exception('Product %s not recognized' % product)
                fout = fout_base + product + suffix
                if len(glob.glob(fout+'.*')) == 0 or overwrite: fouts[product] = fout

            if len(fouts) > 0:
                start = datetime.datetime.now()
                try:
                    # TODO - If only doing temp then don't waste time with other bands
                    img = self._readraw(tile)
                except Exception,e:
                    print 'Error reading data %s' % data['products']['raw']
                    VerboseOut('%s %s' % (data['basename'],e), 2)
                    VerboseOut(traceback.format_exc(), 3)
                    return
                if runatm: atmospheres = [atmosphere(i,data['metadata']) for i in range(1,img.NumBands()+1)]
                VerboseOut('%s: read in %s' % (data['basename'],datetime.datetime.now() - start))

                for product,fname in fouts.iteritems():
                    try:
                        start = datetime.datetime.now()

                        if (self._products[product]['atmcorr']):
                            for i in range(0,img.NumBands()):
                                b = img[i]
                                b.SetAtmosphere(atmospheres[i])
                                img[i] = b
                            VerboseOut('atmospherically correcting',2)
                        else:
                            img.ClearAtmosphere()

                        try:
                            fcall = 'gippy.%s(img, fname)' % self._products[product]['function']
                            VerboseOut(fcall,2)
                            imgout = eval(fcall)
                        except Exception,e:
                            print 'Error creating product %s for %s: %s' % (product,fname,e)
                            VerboseOut(traceback.format_exc(),2)

                        #if overviews: imgout.AddOverviews()
                        #fname = imgout.Filename()
                        imgout = None
                        dur = datetime.datetime.now() - start
                        data['products'][product] = fname
                        VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname),dur))
                    except Exception,e:
                        print 'Error processing %s' % fname
                        VerboseOut('%s %s' % (data['basename'],e),2)
                        VerboseOut(traceback.format_exc(), 3)
                img = None
                # cleanup directory
                try:
                    for bname in data['metadata']['filenames']: 
                        files = glob.glob(os.path.join(data['path'],bname)+'*')
                        for f in files: os.remove(f)
                    shutil.rmtree(os.path.join(dirname,'modtran'))
                except: pass

    def project(self, res=[30,30], datadir='data_landsat'):
        """ Create reprojected, mosaiced images for this site and date """
        if res is None: res=[30,30]
        super(LandsatData, self).project(res=res, datadir=datadir)

    def _readmeta(self, tile, filename=None):
        """ Read in Landsat MTL (metadata) file """
        if filename is None:
            filename = self.tiles[tile]['products']['raw']
        mtlfilename = glob.glob(os.path.join(os.path.dirname(filename),'*MTL.txt'))

        if len(mtlfilename) == 0:
            # Extract MTL file
            if tarfile.is_tarfile(filename):
                tfile = tarfile.open(filename)
            else:
                raise Exception('Not a valid landsat tar file')
            try:
                mtl = ([f for f in tfile.getnames() if "MTL.txt" in f])[0]
            except:
                raise Exception(': possibly an (unsupported) NLAPS processed file')
            tfile.extract(mtl,os.path.dirname(filename))
            mtlfilename = os.path.join(os.path.dirname(filename),mtl)
        else:
            mtlfilename = mtlfilename[0]

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
        dt = datetime.datetime.strptime(mtl['DATE_ACQUIRED'] + ' ' + 
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
            'JulianDay': (dt - datetime.datetime(dt.year,1,1)).days + 1, 
            'DecimalTime': dectime,
        }

        # TODO - now that metadata part of LandsatData object some of these keys not needed
        self.tiles[tile]['metadata'] = {
            'sensor': 'Landsat'+str(id),
            'metafilename': mtlfilename, 
            'filenames': filenames,
            'geometry': _geometry,
            'datetime': _datetime,
            'bands': bandmeta,
            'clouds': clouds
        }
        return self.tiles[tile]['metadata']

    def opentile(self, tile, product=''):
        """ Open and return tile/product GeoImage """
        if product != '':
            return gippy.GeoImage()
          

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
        filename = tiledata['products']['raw']
        if tarfile.is_tarfile(filename):
            tfile = tarfile.open(filename)
        else:
            raise Exception('%s not a valid landsat tar file' % os.path.basename(filename))
            return

        filenames = []
        for b in bandnums:
            fname = tiledata['metadata']['filenames'][b-1]
            if not os.path.exists(fname):
                tfile.extract(fname,tiledata['path'])
            filenames.append(os.path.join(tiledata['path'],fname))

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

def main(): datamain(LandsatData)
