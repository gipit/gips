#!/usr/bin/env python

import os, sys, errno
import argparse
import glob
import re

import datetime
import shutil
import numpy
import tarfile
from copy import deepcopy as deepcopy
import traceback
from collections import OrderedDict

import gippy
from gippy.atmosphere import atmosphere

from gippy.data.core import Data, DataInventory

from pdb import set_trace

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
    ('ind', {
        'description': 'Atmospherically corrected common indices (NDVI,EVI,LSWI,NDSI,BI)',
        'function': 'Indices',
        'atmcorr': True,
    }),
    ('toaind', {
        'description': 'Top of atmosphere common indices (NDVI,EVI,LSWI,NDSI,BI)',
        'function': 'Indices',
        'atmcorr': False,
    }),
])

global _verbose

class LandsatInventory(DataInventory):
    sensors = {'LT4': 'Landsat 4', 'LT5': 'Landsat 5', 'LE7': 'Landsat 7', 'LC8': 'Landsat 8'}
    _colors = {'Landsat 4':'bright yellow', 'Landsat 5':'bright red', 'Landsat 7':'bright green', 'Landsat 8':'bright blue'}
    _rootdir = '/titan/data/landsat/unprocessed'
    _tile_attribute = 'pr'

    @staticmethod
    def get_tile_vector():
        """ Get GeoVector of the sensor tile grid """
        dbstr  = "PG:dbname=geodata host=congo port=5432 user=ags"
        return gippy.GeoVector(dbstr, layer='landsat_wrs')

    def __init__(self, site=None, tiles=None, dates=None, days=None, products=None):
        super(LandsatInventory, self).__init__(site, tiles, dates, days)

        # get all potential matching dates for tiles
        dates = []
        for t in self.tile_names:
            #set_trace()
            for d in os.listdir(self.path(t)):
                date = datetime.datetime.strptime(os.path.basename(d),'%Y%j').date()
                doy = int(date.strftime('%j'))
                if (self.start_date <= date <= self.end_date) and (self.start_day <= doy <= self.end_day): 
                    if date not in dates: dates.append(date)

        # for each date, find all files (define self.data)
        self.numfiles = 0
        self.data = {}
        for date in sorted(dates):
            datedir = date.strftime('%Y%j')
            fnames = []
            for t in self.tile_names: fnames.extend( glob.glob(os.path.join(self.path(t),datedir,'*.tar.gz')) )
            files = []
            for fname in fnames:
                # this is a LandsatData object
                data = LandsatData(fname, products=products)
                files.append(data)
            self.numfiles = self.numfiles + len(fnames)
            # assuming first file is same sensor for all (and for landsat it is except for initial validation flyover L8-L7)
            self.data[date] = { files[0].sensor: files }

    def printcalendar(self,md=False,products=False):
        print 'Landsat Inventory:'
        super(LandsatInventory, self).printcalendar(md,products)
        if self.site is not None:
            print 'Tile Coverage:'
            for t in sorted(self.tiles): print ' P/R %s: %2.0f%%' % (t,self.tiles[t]*100)

class LandsatData(Data):
    """ Represents a single tile and all (existing) product variations (raw, ref, toaref, ind, ndvi, etc) """
    def __init__(self,filename, products=None):
        self.filename = filename
        self.path,basename = os.path.split(filename)
        self.basename = basename[:-12]
        self.tile = self.basename[3:9]
        self.sensor = self.basename[0:3]

        # find all products for this tile
        product_filenames = glob.glob( os.path.join(self.path,'*.tif'))
        self.products = {}
        for f in product_filenames:
            prod = os.path.splitext(os.path.split(f)[1][len(self.basename)+1:])[0]
            self.products[ prod ] = f
        # TODO - prune products based on filter

    def read(self,product='raw'):
        if product == 'raw':
            return self.readraw()
        else:
            try:
                return GeoImage(self.products[product])
            except Exception,e:
                print 'Unknown product %s' % product

    def _readmeta(self):
        """ Read in Landsat MTL (metadata) file """
        mtlfilename = glob.glob(os.path.join(os.path.dirname(self.filename),'*MTL.txt'))

        if len(mtlfilename) == 0:
            # Extract MTL file
            if tarfile.is_tarfile(self.filename):
                tfile = tarfile.open(self.filename)
            else:
                raise Exception('Not a valid landsat tar file')
            try:
                mtl = ([f for f in tfile.getnames() if "MTL.txt" in f])[0]
            except:
                raise Exception(': possibly an (unsupported) NLAPS processed file')
            tfile.extract(mtl,os.path.dirname(self.filename))
            mtlfilename = os.path.join(os.path.dirname(self.filename),mtl)
        else:
            mtlfilename = mtlfilename[0]

        if _verbose > 1: print 'reading %s' % mtlfilename
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
            colors = ["Blue","Green","Red","NIR","SWIR1","LWIR","SWIR2"]
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
            colors = ["Blue","Green","Red","NIR","SWIR1","LWIR","SWIR2"]
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
            colors = ["Coastal","Blue","Green","Red","NIR","SWIR1","SWIR2","Cirrus"] #,"LWIR1","LWIR2"]
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
        self.metadata = {
            'sensor': 'Landsat'+str(id),
            'metafilename': mtlfilename, 
            'filenames': filenames,
            'geometry': _geometry,
            'datetime': _datetime,
            'bands': bandmeta
        }

    def readraw(self,bandnums=[],verbose=False):
        """ Read in Landsat bands using original tar.gz file """
        # Read in collection metadata
        self._readmeta()

        if len(bandnums) != 0: 
            bandnums = numpy.array(bandnums)
        else:
            bandnums = numpy.arange(0,len(self.metadata['bands'])) + 1
        
        # Extract desired files from tarfile
        if tarfile.is_tarfile(self.filename):
            tfile = tarfile.open(self.filename)
        else:
            raise Exception('%s not a valid landsat tar file' % os.path.basename(self.filename))
            return

        filenames = []
        for b in bandnums:
            fname = self.metadata['filenames'][b-1]
            if not os.path.exists(fname):
                tfile.extract(fname,self.path)
            filenames.append(os.path.join(self.path,fname))

        if _verbose > 1: 
            print 'Landsat files: '
            for f in filenames: print ' %s' % f

        # Read in files
        image = gippy.GeoImage(filenames[0])
        del filenames[0]
        for f in filenames: image.AddBand(gippy.GeoImage(f)[0])

        # Geometry used for calculating incident irradiance
        theta = numpy.pi * self.metadata['geometry']['solarzenith']/180.0
        sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (self.metadata['datetime']['JulianDay']-4.0)/180.0))

        # Set metadata
        image.SetNoData(0)
        image.SetUnits('radiance')
        # TODO - set appropriate metadata
        #for key,val in meta.iteritems():
        #    image.SetMeta(key,str(val))

        for b in bandnums:
            bandmeta = self.metadata['bands'][b-1]
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

def process(inventory, products=['toarad'], overwrite=False, suffix='', overviews=False):
    
    if inventory == None: raise Exception('No data inventory provided')
    if products == None: raise Exception('No products specified')

    if _verbose > 0:
        print 'Requested %s products for %s files' % (len(products), inventory.numfiles)
    
    if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
    for date in inventory.dates:
        for sensor in inventory.data[date]:
            for dat in inventory.data[date][sensor]:
                fout_base = os.path.join(dat.path, dat.basename + '_')
                fouts = {}
                runatm = False
                # generate all product names and check if already existing
                for product in products:
                    if _products[product]['atmcorr']: runatm = True
                    if product not in _products.keys():
                        raise Exception('Product %s not recognized' % product)
                    fout = fout_base + product + suffix
                    if len(glob.glob(fout+'.*')) == 0 or overwrite: fouts[product] = fout

                if len(fouts) > 0:
                    start = datetime.datetime.now()
                    try:
                        # TODO - If only doing temp then don't waste time with other bands
                        img = dat.read()
                    except Exception,e:
                        print 'Error reading data %s' % dat.filename
                        if _verbose > 1:
                            print '%s %s' % (dat.basename,e)
                            if _verbose > 2: print traceback.format_exc()
                        continue
                    if runatm: atmospheres = [atmosphere(i,dat.meta) for i in range(1,img.NumBands()+1)]
                    if _verbose > 0:
                        print '%s: read in %s' % (dat.basename,datetime.datetime.now() - start)

                    for product,fname in fouts.iteritems():
                        try:
                            start = datetime.datetime.now()

                            if (_products[product]['atmcorr']):
                                for i in range(0,img.NumBands()):
                                    b = img[i]
                                    b.SetAtmosphere(atmospheres[i])
                                    img[i] = b
                                if _verbose > 1: print 'atmospherically correcting'
                            else:
                                img.ClearAtmosphere()

                            try:
                                fcall = 'gippy.%s(img, fname)' % _products[product]['function']
                                if _verbose > 1: print fcall
                                imgout = eval(fcall)
                            except Exception,e:
                                print 'Error creating product %s for %s: %s' % (product,fname,e)
                                if _verbose > 1: print traceback.format_exc()                              

                            if overviews: imgout.AddOverviews()
                            fname = imgout.Filename()
                            imgout = None
                            if _verbose > 0:
                                dur = datetime.datetime.now() - start
                                print ' -> %s: processed in %s' % (os.path.basename(fname),dur)
                        except Exception,e:
                            print 'Error processing %s' % fname
                            if _verbose > 1:
                                print '%s %s' % (dat.basename,e)
                                if _verbose > 2: print traceback.format_exc()
                    img = None
                    # cleanup directory
                    try:
                        for bname in dat.metadata['filenames']: 
                            files = glob.glob(os.path.join(dat.path,bname)+'*')
                            for f in files: os.remove(f)
                        shutil.rmtree(os.path.join(dirname,'modtran'))
                    except: pass
    if _verbose > 0: print 'Completed processing'
    # TODO - multithreaded ?
    #if args.multi:
    #    for f in fnames:
    #        pool = multiprocessing.Pool(cpus)
    #        pool.apply_async(landsat.process, 
    #            [f,args.units, args.atmcorr, args.product, args.datatype, args.verbose, args.overwrite, args.suffix])
    #    pool.close()
    #    pool.join()
        #multiprocessing.Process(target=Process, 
        #    args=(f, args.units,args.atmcorr, args.product, args.datatype, args.verbose))
        #jobs.append(p)
        #p.start()
    #else:

def project(inventory, site, res):
    for date in inventory.dates:
        for sensor in inventory.data[date]:
            products = inventory.data[date][sensor][0].products.keys()
            for data in inventory.data[date][sensor]:
                if products != data.products.keys():
                    raise Exception('SOMETHING IS AFOOT')
            for p in products:
                start = datetime.datetime.now()
                fnamesin = []
                for data in inventory.data[date][sensor]:
                    fnamesin.append( data.products[p] )
                fnameout = date.strftime('%Y%j') + '_%s_%s' % (p,sensor)
                imgout = gippy.CookieCutter(fnamesin, fnameout, site, res[0], res[1])
                print 'Merged %s files -> %s in %s' % (len(fnamesin),imgout.Basename(),datetime.datetime.now() - start)

# This is broken
def archive(dir=''):
    if (dir == ''):
        fnames = glob.glob('L*.tar.gz')
    else:
        fnames = glob.glob(os.path.join('L*.tar.gz'))
    numadded = 0
    for f in fnames:
        pathrow = f[3:9]
        year = f[9:13]
        doy = f[13:16]
        path = os.path.join(origpath,pathrow,year+doy)
        try:
            os.makedirs(path)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise Exception('Unable to make product directory %s' % path)
        try:
            newf = os.path.join(path,f)
            if os.path.exists(newf):
                print '%s: already in archive' % f
                os.remove(f)
            else:
                shutil.move(f,newf)
                try:
                    #mtl = _readmeta(newf)
                    print f, ' -> ',path
                    numadded = numadded + 1
                except Exception,e:
                    print f, ' -> problem with file'
                    unarchive(os.path.dirname(newf))
        except shutil.Error as err:
            print err
            #raise Exception('shutils error %s' % err)
            #if exc.errno == errno.EEXIST:
            #    print f, ' removed, already in archive'
            #else:
            #    pass
    print '%s files added to archive' % numadded
    if numadded != len(fnames):
        print '%s files not added to archive' % (len(fnames)-numadded)

def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='Landsat Processing Utilities', formatter_class=argparse.RawTextHelpFormatter)
    subparser = parser0.add_subparsers(dest='command')

    # Global options
    gparser = argparse.ArgumentParser(add_help=False, formatter_class=dhf)
    group = gparser.add_argument_group('Data Inventory Options')
    group.add_argument('-s','--site',help='Vector file for region of interest', default=None)
    group.add_argument('-d','--dates',help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    group.add_argument('--days',help='Include only those that fall within these days of year (doy1,doy2)',default=None)
    group.add_argument('--tiles', nargs='*', help='the Landsat pathrow(s)', default=[])
    group.add_argument('-p','--products', nargs='*', help='Process/filter these products') #default=False)
    group.add_argument('-v','--verbose', help='Verbosity level', default=0, type=int)

    # Help
    parser = subparser.add_parser('help',help='Print extended help', parents=[gparser], formatter_class=dhf)

    # Inventory
    parser = subparser.add_parser('inventory',help='Get Landsat Inventory', parents=[gparser], formatter_class=dhf)
    parser.add_argument('--md',help='Show dates using MM-DD',action='store_true',default=False)

    # Processing
    parser = subparser.add_parser('process',help='Process Landsat scenes', parents=[gparser],formatter_class=dhf)
    group = parser.add_argument_group('Processing Options')
    group.add_argument('--overwrite', help='Overwrite output files if they exist', default=False, action='store_true')
    group.add_argument('--suffix', help='Append string to end of filename (before extension)',default='')
    #group.add_argument('--nooverviews', help='Do not add overviews to output', default=False, action='store_true')
    #pparser.add_argument('--link', help='Create links in current directory to output', default=False, action='store_true')
    #pparser.add_argument('--multi', help='Use multiple processors', default=False, action='store_true')

    # Project
    parser = subparser.add_parser('project',help='Create project', parents=[gparser], formatter_class=dhf)
    group = parser.add_argument_group('Project options')
    group.add_argument('--res',nargs=2,help='Resolution of output rasters', default=[30,30], type=float)

    # Links
    parser = subparser.add_parser('link',help='Link to Landsat Products', parents=[gparser], formatter_class=dhf)
    parser.add_argument('--hard',help='Create hard links instead of symbolic', default=False,action='store_true')

    # Misc
    parser_archive = subparser.add_parser('archive',help='Move files from this directory to Landsat archive')

    args = parser0.parse_args()

    if args.command == 'help':
        parser0.print_help()
        print '\navailable products:'
        for key,val in _products.items(): 
            print '    {:<20}{:<100}'.format(key, val['description'])
        exit(1)

    global _verbose 
    _verbose = args.verbose
    gippy.Options.SetVerbose(_verbose)
    gippy.Options.SetChunkSize(128.0)   # replace with option

    try:
        inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)
    except Exception,e:
        print 'Error getting inventory: %s' % (e)
        exit(1)

    if args.command == 'inventory':
        if args.products is None:
            inv.printcalendar(args.md)
        else: inv.printcalendar(args.md,True)
        
    elif args.command == 'link':
        inv.createlinks(args.hard)

    elif args.command == 'process':
        try:
            #merrafname = fetchmerra(meta['datetime'])
            process(inv,products=args.products,overwrite=args.overwrite,suffix=args.suffix) #, nooverviews=args.nooverviews)
            #inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)
        except Exception,e:
            print 'Error processing: %s' % e
            if _verbose > 2: print traceback.format_exc()

    elif args.command == 'project':
        process(inv,products=args.products) 
        inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)
        project(inv, args.site, args.res)

    elif args.command == 'archive': archive()
    #elif args.command == 'clean': clean()
    else:
        print 'Command %s not recognized' % cmd
