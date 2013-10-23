#!/usr/bin/env python

import os, sys, errno
import argparse
import glob
import re
import agspy.utils.dateparse as dateparse
import datetime
import shutil
import numpy
import tarfile
from copy import deepcopy as deepcopy
import traceback

import gippy
from gippy.atmosphere import atmosphere

from gippy.data.core import DataInventory

from pdb import set_trace

class LandsatInventory(DataInventory):
    sensors = {'LT4': 'Landsat 4', 'LT5': 'Landsat 5', 'LE7': 'Landsat 7', 'LC8': 'Landsat 8'}
    _colors = {'Landsat 4':'bright yellow', 'Landsat 5':'bright red', 'Landsat 7':'bright green', 'Landsat 8':'bright blue'}
    _rootdir = '/titan/data/landsat'
    _origdir = 'unprocessed'
    _proddir = 'products/gippy-beta'
    _tile_attribute = 'pr'

    @staticmethod
    def get_tile_vector():
        """ Get GeoVector of the sensor tile grid """
        dbstr  = "PG:dbname=geodata host=congo port=5432 user=ags"
        return gippy.GeoVector(dbstr, layer='landsat_wrs')

    def __init__(self, site=None, tiles=None, dates=None, days=None, products=None):
        # Spatial extent (define self.tiles)
        if tiles is None: tiles = os.listdir(self.origpath())
        self.site = site
        if self.site is not None: 
            self.site_to_tiles(gippy.GeoVector(self.site))
        else: self.tiles = dict((t,1) for t in tiles)

        # Temporal extent (define self.dates and self.days)
        if dates is None: dates='1984,2050'
        self.start_date,self.end_date = dateparse.range(dates)
        if days: 
            days = days.split(',')
        else: days = (1,366)
        self.start_day,self.end_day = ( int(days[0]), int(days[1]) )

        # Products of interest (define self.products)
        if products is None:
            products = ['']
        elif len(products) == 0:
            products = ['']
        self.products = products

        # get all potential matching dates for tiles
        dates = []
        for t in self.tile_names:
            for d in os.listdir(self.origpath(t)):
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
            for t in self.tile_names: fnames.extend( glob.glob(os.path.join(self.origpath(t),datedir,'*.tar.gz')) )
            files = []
            for fname in fnames:
                path,basename = os.path.split(fname)
                basename = basename[:-12]
                tile = basename[3:9]
                # find products for this tile
                product_path = os.path.split(fname)[0].replace(self.origpath(),self.prodpath())
                product_fnames = []
                products = {}
                for p in self.products:
                    product_fnames.extend( glob.glob( os.path.join(product_path,'*'+p+'.tif') ) )
                for f in product_fnames:
                    prod = os.path.splitext(os.path.split(f)[1][len(basename)+1:])[0]
                    products[ prod ] = f
                #f = fname.split()
                e = {'tile': tile, 'basename':basename, 'filename': fname, 'products': products, 'productpath':product_path}
                files.append(e)
            self.numfiles = self.numfiles + len(fnames)
            self.data[date] = { basename[0:3]: files }
            #color':self._colors['Landsat'+basename[2:3]]

    def printcalendar(self,md=False,products=False):
        print 'Landsat Inventory:'
        super(LandsatInventory, self).printcalendar(md,products)
        if self.site is not None:
            print 'Tile Coverage:'
            for t in sorted(self.tiles): print ' P/R %s: %2.0f%%' % (t,self.tiles[t]*100)


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




def process(inventory, products=['rad'], verbose=0, overwrite=False, suffix='', overviews=False):

    print 'Requested %s products for %s files' % (len(products), inventory.numfiles)
    
    if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix

    for date in inventory.dates:
        for sensor in inventory.data[date]:
            for dat in inventory.data[date][sensor]:
                fout_base = os.path.join(dat['productpath'], dat['basename'] + '_')

                # Make output directory
                if not os.path.isdir(dat['productpath']):
                    try:
                        os.makedirs(dat['productpath'])
                    except OSError as exc:
                        raise Exception('Unable to make product directory %s' % dat['productpath'])

                # Copy MTL file
                #mtlout = meta['metafilename'].replace(origpath,prodpath)
                #mtl = glob.glob(os.path.join(os.path.dirname(fin),'*MTL.txt'))[0]
                #mtlout = mtl.replace(origpath,prodpath)
                #if not os.path.exists(os.path.dirname(mtlout)):
                #    os.makedirs(os.path.dirname(mtlout))
                #try:
                #    shutil.copy(mtl,mtlout)
                #except: pass

                fouts = {}
                runatm = False
                # generate all product names
                for product in products:
                    if product[0:2] == 'A_':
                        runatm = True
                    fout = fout_base + product + suffix
                    if len(glob.glob(fout+'.*')) == 0 or overwrite: fouts[product] = fout

                if len(fouts) > 0:
                    # Read in data
                    start = datetime.datetime.now()
                    try:
                        # TODO - If only doing temp then don't waste time with other bands
                        img = read(dat['filename'], verbose=verbose)
                        # TODO - shouldn't have to _readmeta again
                        meta = _readmeta(dat['filename'])
                    except Exception,e:
                        print '%s %s' % (dat['basename'],e)
                        if verbose: print traceback.format_exc()
                        continue
                    if runatm:
                        atmospheres = [atmosphere(i,meta) for i in range(1,img.NumBands()+1)]
                    print '%s: read in %s' % (dat['basename'],datetime.datetime.now() - start)

                    for product,fname in fouts.iteritems():
                        try:
                            start = datetime.datetime.now()

                            if (product[0:1] == 'A'):
                                for i in range(0,img.NumBands()): img[i].SetAtmosphere(atmospheres[i])
                            else:
                                img.ClearAtmosphere()
                            product = product[2:]
                            if product == 'cind':
                                imgout = gippy.Indices(img,fname)
                            elif product == 'rgb':
                                img.PruneToRGB()
                                imgout = gippy.RGB(img,fname)
                            elif product == 'radi':
                                imgout = gippy.Rad(img,fname)
                            elif product == 'refl':
                                imgout = gippy.Ref(img,fname)
                            elif product == 'raw':
                                imgout = gippy.Copy(img,fname)
                            elif product == 'fmask':
                                imgout = gippy.algorithms.fmask.process(img,fname)
                            else:
                                raise Exception('product %s not recognized' % product)

                            if overviews: imgout.AddOverviews()
                            fname = imgout.Filename()
                            imgout = None
                            
                            dur = datetime.datetime.now() - start
                            print ' -> %s: processed in %s' % (os.path.basename(fname),dur)
                        except Exception,e:
                            print '%s %s' % (dat['basename'],e)
                            if verbose > 1: print traceback.format_exc()
                    img = None
                    _cleandir(os.path.dirname(dat['filename']))
    print 'Completed processing'
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

def project(inventory, site, res):
    for date in inventory.dates:
        for sensor in inventory.data[date]:
            products = inventory.data[date][sensor][0]['products'].keys()
            for data in inventory.data[date][sensor]:
                if products != data['products'].keys():
                    raise Exceptio
                    n('SOMETHING IS AFOOT')
            for p in products:
                start = datetime.datetime.now()
                fnamesin = []
                for data in inventory.data[date][sensor]:
                    fnamesin.append( data['products'][p] )
                fnameout = date.strftime('%Y%j') + '_%s_%s' % (p,sensor)
                imgout = gippy.CookieCutter(fnamesin, fnameout, site, res[0], res[1])
                print 'Merged %s files -> %s in %s' % (len(fnamesin),imgout.Basename(),datetime.datetime.now() - start)

#def unarchive(path):
#    try:
#        datafile = glob.glob(os.path.join(path,'*tar.gz'))[0]
#        bname = os.path.basename(datafile)
#        newf = os.path.join(qdir,bname)
#        shutil.move(datafile,newf)
#        print '%s: removed from archive' % os.path.basename(bname)
#    except Exception,e: pass
#    shutil.rmtree(path)

#def clean():
#    """ Clean out archive directories, keeping tar.gz and mtl files """
#    prdirs = os.listdir(origpath)
#    print 'Cleaning Landsat archive: %s pathrows' % len(prdirs)
#    for pr in prdirs:
#        print 'Cleaning pathrow %s' % pr
#        datedirs = os.listdir(os.path.join(origpath,pr))
#        for dd in datedirs:
#            _cleandir(os.path.join(origpath,pr,dd))

if __name__ == "__main__":
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='Landsat Archive Utilities') #, formatter_class=dhf)
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
    #parser_clean = subparser.add_parser('clean',help='Clean archive of all temporary files')

    args = parser0.parse_args()

    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetChunkSize(128.0)   # replace with option

    inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)
    
    if args.command == 'inventory':
        if args.products is None:
            inv.printcalendar(args.md)
        else: inv.printprodcal(args.md)
        
    elif args.command == 'link':
        inv.createlinks(args.hard)

    elif args.command == 'process':
        #merrafname = fetchmerra(meta['datetime'])
        process(inv,products=args.products,verbose=args.verbose,overwrite=args.overwrite,suffix=args.suffix) #, nooverviews=args.nooverviews)
        #inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)

    elif args.command == 'project':
        process(inv,products=args.products,verbose=args.verbose) 
        inv = LandsatInventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)
        project(inv, args.site, args.res)

    elif args.command == 'archive': archive()
    #elif args.command == 'clean': clean()
    else:
        print 'Command %s not recognized' % cmd