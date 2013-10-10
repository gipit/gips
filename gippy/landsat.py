#!/usr/bin/env python

import os, sys, errno
import argparse
import glob
import re
import agspy.utils.dateparse as dateparse
import datetime
import os, errno
import shutil
import numpy
import tarfile
from copy import deepcopy as deepcopy
import traceback

import gippy
from gippy.atmosphere import atmosphere
from gippy.GeoVector import GeoVector
import ogr, osr
from shapely.wkb import loads

from pdb import set_trace

# VARIABLES TO EDIT ################################################



####################################################################



# This should only be called from read function
def _readmeta(filename):
    """ Read in Landsat MTL (metadata) file """
    if filename[-7:] != 'MTL.txt':
        mtl = glob.glob(os.path.join(os.path.dirname(filename),'*MTL.txt'))
        if len(mtl) == 0:
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
            filename = os.path.join(os.path.dirname(filename),mtl)
        else:
            filename = mtl[0]

    # Read MTL file
    try:
        text = open(filename,'r').read()
    except IOError as e:
        raise Exception('({})'.format(e))

    # Find out what Landsat this is
    try:
        id = int(os.path.basename(filename)[1:2])
    except:
        id = int(os.path.basename(filename)[2:3])

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

    return {
        'sensor': 'Landsat'+str(id),
        'metafilename': filename, 
        'filenames': filenames,
        'geometry': _geometry,
        'datetime': _datetime,
        'bands': bandmeta
    }

def read(filename, bandnums=[],verbose=False):
    """ Read in Landsat bands using original tar.gz file """
    
    # Read in collection metadata
    meta = _readmeta(filename)

    if len(bandnums) != 0: 
        bandnums = numpy.array(bandnums)
    else:
        bandnums = numpy.arange(0,len(meta['bands'])) + 1

    # Extract desired files from tarfile
    if tarfile.is_tarfile(filename):
        tfile = tarfile.open(filename)
    else:
        raise Exception('%s not a valid landsat tar file' % os.path.basename(filename))
        return
    filenames = []
    for b in bandnums:
        dirname = os.path.dirname(filename)
        fname = meta['filenames'][b-1]
        if not os.path.exists(fname):
            tfile.extract(fname,dirname)
        filenames.append(os.path.join(dirname,fname))

    # Read in files
    image = gippy.GeoImage(filenames[0])
    del filenames[0]
    for f in filenames: image.AddBand(gippy.GeoImage(f)[0])

    # Geometry used for calculating incident irradiance
    theta = numpy.pi * meta['geometry']['solarzenith']/180.0
    sundist = (1.0 - 0.016728 * numpy.cos(numpy.pi * 0.9856 * (meta['datetime']['JulianDay']-4.0)/180.0))

    # Set metadata
    image.SetNoData(0)
    image.SetUnits('radiance')
    # TODO - set appropriate metadata
    #for key,val in meta.iteritems():
    #    image.SetMeta(key,str(val))

    for b in bandnums:
        bandmeta = meta['bands'][b-1]
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

def link(f,hard=False):
    """ Create link to file in current directory """
    faux = f + '.aux.xml'
    if hard:
        try:
            os.link(f,os.path.basename(f))
            os.link(faux,os.path.basename(faux))
        except:
            pass
    else: 
        try:
            os.symlink(f,os.path.basename(f))
            if os.path.isfile(faux):
                os.symlink(faux,os.path.basename(faux))
        except:
            pass

class DataInventory(object):
    _rootdir = ''
    _origdir = ''
    _proddir = ''

    _colorcodes = {
        'black':    '0;30',     'bright gray':  '0;37',
        'blue':     '0;34',     'white':        '1;37',
        'green':    '0;32',     'bright blue':  '1;34',
        'cyan':     '0;36',     'bright green': '1;32',
        'red':      '0;31',     'bright cyan':  '1;36',
        'purple':   '0;35',     'bright red':   '1;31',
        'yellow':   '0;33',     'bright purple':'1;35',
        'dark gray':'1;30',     'bright yellow':'1;33',
        'normal':   '0'
    }

    def __getitem__(self,date):
        return self.files[date]

    def origpath(self):
        return os.path.join(self._rootdir,self._origdir)

    def prodpath(self):
        return os.path.join(self._rootdir,self._proddir)

    def tilepath(self,tile): 
        return os.path.join(self.origpath(),tile)

    def _colorize(self,txt,color): 
        return "\033["+self._colorcodes[color]+'m' + txt + "\033[0m"

    def get_tiles(self):
        return [k for k,v in self.tiles.items()]

    def get_dates(self):
        return [k for k,v in self.files.items()]

    def size(self): 
        return len(self.files)

    def site_to_tiles(self,vector):
        """ Identify sensor tiles that fall within vector """
        geom = vector.union()
        ogrgeom = ogr.CreateGeometryFromWkb(geom.wkb)
        tvector = self.get_tile_vector()
        tlayer = tvector.layer
        tlayer.SetSpatialFilter(ogrgeom)
        tiles = {}
        tlayer.ResetReading()
        feat = tlayer.GetNextFeature()
        fldindex = feat.GetFieldIndex(self._tile_attribute)
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkb())
            area = geom.intersection(tgeom).area
            if area != 0: 
                tile = str(feat.GetField(fldindex))
                if len(tile) == 5: tile = '0' + tile
                tiles[tile] = area/geom.area
            feat = tlayer.GetNextFeature()  
        self.tiles = tiles
        return tiles

    def printcalendar(self,inventories,md=False, products=False):
        #import calendar
        #cal = calendar.TextCalendar()
        oldyear = ''
        for d in sorted(inventories):        
            if args.md:
                daystr = str(d.month) + '-' + str(d.day)
            else:
                daystr = str(d.timetuple().tm_yday)
                if len(daystr) == 1:
                    daystr = '00' + daystr
                elif len(daystr) == 2:
                    daystr = '0' + daystr
            if d.year != oldyear:
                sys.stdout.write('\n{:>5}:'.format(d.year))
                if products: sys.stdout.write('\n')
            sys.stdout.write(_colorize('{:>7}'.format(daystr),inventories[d]['color']))
            if products:
                sys.stdout.write('        ')
                for v in inventories[d]['products']:
                    sys.stdout.write(_colorize('{:<9}'.format(v),inventories[d]['color']))
                sys.stdout.write('\n')
            oldyear = d.year
        sys.stdout.write('\n')

    #def __str__(self):
    #    s = self.date.strftime('%Y-%j') + ':'
    #    for f in self.files: s = s + ' ' + f['basename']
    #    return s

class LandsatInventory(DataInventory):
    _colors = {4:'bright yellow', 5:'bright red', 7:'bright green', 8:'bright blue'}
    _rootdir = '/titan/data/landsat'
    _origdir = 'unprocessed'
    _proddir = 'products/gippy-beta'
    _tile_attribute = 'pr'

    @staticmethod
    def get_tile_vector():
        """ Get GeoVector of the sensor tile grid """
        dbstr  = "PG:dbname=geodata host=congo port=5432 user=ags"
        return GeoVector(dbstr, layer='landsat_wrs')

    def __init__(self, site=None, dates=None, days=None, tiles=None):
        if tiles is None: tiles = os.listdir(self.origpath)
        if dates is None: dates='1984,2050'
        if args.site is not None: 
            self.site_to_tiles(GeoVector(args.site))
        else: self.tiles = dict((t,0) for t in tiles)
        self.start_date, self.end_date = dateparse.range(dates)
        if days: 
            days = days.split(',')
            self.days = ( int(days[0]), int(days[1]) )
        else: self.start_day, self.end_day = (1,366)
        # get all potential matching dates
        dates = []
        for t in tiles: 
            for d in os.listdir(self.tilepath(t)):
                date = datetime.datetime.strptime(os.path.basename(d),'%Y%j').date()
                doy = int(date.strftime('%j'))
                if (self.start_date <= date <= self.end_date) and (self.start_day <= doy <= self.end_day): 
                    if date not in dates: dates.append(date)
        # for each date, find all files
        self.files = {}
        for date in sorted(dates):
            datedir = date.strftime('%Y%j')
            fnames = []
            for t in tiles: fnames.extend( glob.glob(os.path.join(self.tilepath(t),datedir,'*.tar.gz')) )
            files = []
            for fname in fnames:
                basename = os.path.basename(fname)
                tile = basename[3:9]
                product_path = os.path.split(fname)[0].replace(self.origpath(),self.prodpath())
                product_fnames = glob.glob(os.path.join(product_path, '*.tif'))
                products = [os.path.basename(f)[17:-4] for f in product_fnames]
                e = {'filename': fname, 'basename':basename, 'products': sorted(products), 
                    'tile': tile, 'color':self._colors[int(basename[2:3])]}
                files.append(e)
            self.files[date] = files

def outfile(filename, product="rad", atmcorr=False, overwrite=False, suffix=""):
    """ Generate output filename and create needed directories """
    # Generate output filename
    basename = filename[:-12]
    if atmcorr:
        atmtag = 'A'
    else: atmtag = 'N'
    if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
    ofile = basename.replace(origpath,prodpath) \
        + '_' + atmtag + '_' + product.lower() + suffix

    # Check for existing matching file
    existing_ofiles = glob.glob(ofile+'.*')
    if len(existing_ofiles) > 0 and overwrite == False:
        raise Exception('-> %s: output already exists' % os.path.basename(ofile))

    # Make output directory
    try:
        path = os.path.dirname(ofile)
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise Exception('Unable to make product directory')
    return ofile

def process(img, fname_out, product='rad', verbose=1, overviews=False):
    """ Process an image into various products """
    gippy.Options.SetVerbose(verbose)
    gippy.Options.SetChunkSize(64.0)

    product = product.lower()
    #if product == 'cind': product = 'ind'
    #if product == 'radi': product = 'rad'
    #if product == 'refl': product = 'ref'

    if product == 'cind':
        imgout = gippy.Indices(img,fname_out)
    elif product == 'rgb':
        img.PruneToRGB()
        imgout = gippy.RGB(img,fname_out)
    elif product == 'radi':
        imgout = gippy.Rad(img,fname_out)
    elif product == 'refl':
        imgout = gippy.Ref(img,fname_out)
    elif product == 'fmask':
        imgout = gippy.algorithms.fmask.process(img,fname_out)
    else:
        raise Exception('product %s not recognized' % product)

    if overviews: imgout.AddOverviews()
    fname = imgout.Filename()
    imgout = None
    return fname

def batchprocess(fnames, products=['rad'], atmcorr=False, 
    verbose=1, overwrite=False, suffix='', overviews=False):
    
    for fin in fnames:
        # Generate all fout names
        fouts = []
        for product in products:
            try:
                fname = outfile(fin, product=product, atmcorr=atmcorr,overwrite=overwrite, suffix=suffix)
                fouts.append((product,fname))
            except Exception,e:
                print '%s %s' % (os.path.basename(fin)[:16],e)

        if len(fouts) > 0:
            # Copy MTL file
            #mtlout = meta['metafilename'].replace(origpath,prodpath)
            #mtl = glob.glob(os.path.join(os.path.dirname(fin),'*MTL.txt'))[0]
            #mtlout = mtl.replace(origpath,prodpath)
            #if not os.path.exists(os.path.dirname(mtlout)):
            #    os.makedirs(os.path.dirname(mtlout))
            #try:
            #    shutil.copy(mtl,mtlout)
            #except: pass
            start = datetime.datetime.now()
            try:
                # TODO - If only doing temp then don't waste time with other bands
                img = read(fin, verbose=verbose)
                # TODO - shouldn't have to _readmeta again
                meta = _readmeta(fin)
            except Exception,e:
                print '%s %s' % (os.path.basename(fin)[:16],e)
                if verbose > 1: print traceback.format_exc()
                continue

            if atmcorr:
                for i in range(0,img.NumBands()):
                    band = img[i]
                    band.SetAtmosphere(atmosphere(i+1,meta))
                    img[i] = band

            print '%s: read in %s' % (os.path.basename(fin)[:-12],datetime.datetime.now() - start)
            for fout in fouts:
                try:
                    start = datetime.datetime.now()
                    process(img, fout[1], fout[0], overviews=overviews)
                    dur = datetime.datetime.now() - start
                    print '  -> %s: processed in %s' % (os.path.basename(fout[1]),dur)
                except Exception,e:
                    print '%s %s' % (os.path.basename(fin)[:16],e)
                    if verbose > 1: print traceback.format_exc()
            img = None
        _cleandir(os.path.dirname(fin))
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

#def unarchive(path):
#    try:
#        datafile = glob.glob(os.path.join(path,'*tar.gz'))[0]
#        bname = os.path.basename(datafile)
#        newf = os.path.join(qdir,bname)
#        shutil.move(datafile,newf)
#        print '%s: removed from archive' % os.path.basename(bname)
#    except Exception,e: pass
#    shutil.rmtree(path)

def _cleandir(dir):
    try:
        shutil.rmtree(os.path.join(dir,'modtran'))
    except: pass
    # Remove extraneous files
    fnames = glob.glob(os.path.join(dir,'*'))
    for f in fnames:
        extension = os.path.splitext(f)[1][1:].strip()
        if extension != 'txt' and extension != 'gz' and extension != 'tar':
            os.remove(f)

def clean():
    """ Clean out archive directories, keeping tar.gz and mtl files """
    prdirs = os.listdir(origpath)
    print 'Cleaning Landsat archive: %s pathrows' % len(prdirs)
    for pr in prdirs:
        print 'Cleaning pathrow %s' % pr
        datedirs = os.listdir(os.path.join(origpath,pr))
        for dd in datedirs:
            _cleandir(os.path.join(origpath,pr,dd))



if __name__ == "__main__":
    prog = os.path.split(__file__)[1]
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(prog=prog, description='Landsat Archive Utilities', formatter_class=dhf)
    subparser = parser0.add_subparsers(dest='command')

    # Global options
    gparser = argparse.ArgumentParser(add_help=False, formatter_class=dhf)
    gparser.add_argument('-s','--site',help='Vector file for region of interest', default=None)
    gparser.add_argument('-d','--dates',help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    gparser.add_argument('--days',help='Include only those that fall within these days of year (doy1,doy2)',default=None)
    gparser.add_argument('--tiles', nargs='*', help='the Landsat pathrow(s)', default=[])

    # Inventory
    parser = subparser.add_parser('inventory',help='Get Landsat Inventory', parents=[gparser], formatter_class=dhf)
    parser.add_argument('--md',help='Show dates using MM-DD',action='store_true',default=False)
    parser.add_argument('-p','--products', help='Also list all products available',action='store_true',default=False)

    # Links
    parser = subparser.add_parser('link',help='Link to Landsat Products', parents=[gparser], formatter_class=dhf)
    parser.add_argument('-p','--product',nargs='*',help='Limit links to product (e.g., refl)',default="")
    #parser.add_argument('--hard',help='Create hard links instead of symbolic  (titan only)', default=False,action='store_true')

    # Misc
    parser_archive = subparser.add_parser('archive',help='Move files from this directory to Landsat archive')
    parser_clean = subparser.add_parser('clean',help='Clean archive of all temporary files')

    # General processing options
    pparser = argparse.ArgumentParser(add_help=False,formatter_class=dhf)
    pparser.add_argument('--overviews', help='Add Overviews to output', default=False, action='store_true')
    pparser.add_argument('-v','--verbose', help='Verbosity level', default=0, type=int)
    pparser.add_argument('--overwrite', help='Overwrite output files if they exist', default=False, action='store_true')
    pparser.add_argument('--suffix', help='Append string to end of filename (before extension)',default='')
    #pparser.add_argument('--link', help='Create links in current directory to output', default=False, action='store_true')
    #pparser.add_argument('--multi', help='Use multiple processors', default=False, action='store_true')

    parser = subparser.add_parser('process',help='Process Landsat scenes', parents=[gparser,pparser],formatter_class=dhf)
    parser.add_argument('-a','--atmcorr',action='store_true',help='Atmospheric Correction', default=False)
    parser.add_argument('-p','--product',nargs='*',help='Output product to create (radi,refl,temp,cind,fmask)',default=['radi'])

    args = parser0.parse_args()

    # Getting list of tiles

    inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles)

    set_trace()
    if inv.size() == 0:
        print 'No data matching criteria in inventory'
        exit(1)

    if args.command == 'inventory':
        print 'Landsat Inventory'
        _printsimplecal(dates,args.md,args.products)
        for s in sorted(_colors): print _colorize(' Landsat%s' % s, _colors[s])
        print 'Total of %i matching data files on %s dates' % (len(fnames),len(dates))
        if args.site is not None:
            print 'Tile Coverage:'
            for t,cov in sorted(tile_coverage.items()): print ' P/R %s: %2.0f%%' % (t,cov*100)
        
    elif args.command == 'link':
        for d in dates:
            print d
            for f in fnames: link(f,hard=args.hard)
        #for pr in args.pathrow:
        #    fnames = inventory(pr,args.dates,products=True,days=args.days)
        #    fnames = [f for f in fnames if args.product in f]

    elif args.command == 'process':
        #merrafname = fetchmerra(meta['datetime'])
        for pr in args.pathrow:
            fnames = inventory(pr,args.dates,doy=args.doy)
            print 'Processing %s products for %s files' % (len(args.product), len(fnames))
            batchprocess(fnames, products=args.product, atmcorr=args.atmcorr, 
                verbose=args.verbose, overwrite=args.overwrite, suffix=args.suffix, overviews=args.overviews)

    elif args.command == 'archive': archive()
    elif args.command == 'clean': clean()
    else:
        print 'Command %s not recognized' % cmd