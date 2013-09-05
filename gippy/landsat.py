#!/usr/bin/env python

import os, sys, errno
import glob
import re
import agspy.utils.utils as utils
import datetime
import os, errno
import shutil
import numpy
import tarfile
from copy import deepcopy as deepcopy

import gippy
from gippy.atmosphere import atmosphere

import traceback

#import agspy.config.datadir as datadir

topdir = '/titan/data/landsat'
import socket
host = socket.gethostname()
if host != 'nile' and host !='volga' and host != 'congo': topdir = '/mnt' + topdir
origdir = 'unprocessed'
proddir = 'products/gippy-beta'

rawdir = os.path.join(topdir,origdir)
outdir = os.path.join(topdir,proddir)
qdir = os.path.join(topdir,'quarantine')

def readmeta(filename):
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

    dirname = os.path.dirname(filename)
    
    # Read in collection metadata
    meta = readmeta(filename)

    if len(bandnums) != 0: 
        bandnums = numpy.array(bandnums)
    else:
        bandnums = numpy.arange(0,len(meta['bands'])) + 1

    # Extract desired files from tarfile
    # TODO - Check to see if already extracted
    if tarfile.is_tarfile(filename):
        tfile = tarfile.open(filename)
    else:
        raise Exception('%s not a valid landsat tar file' % os.path.basename(filename))
        return
    filenames = []
    for b in bandnums:
        fname = meta['filenames'][b-1]
        if not os.path.exists(fname):
            #if verbose > 1: print '\tExtracting ',fname
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

def link(pathrow,dates=None,hard=False,filt=''):
    """ Create links to matched products in current directory """
    fnames = inventory(pathrow,dates,products=True)
    fnames = [f for f in fnames if filt in f]
    for f in fnames:
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

def inventory(pathrow,dates=None,products=False):
    """ Get listing of filenames that match pathrow and date range """
    # Raw data
    if products == True:
        fnames = glob.glob(os.path.join(outdir, pathrow, '*', '*.tif'))
    else: # Raw data
        fnames = glob.glob(os.path.join(rawdir, pathrow, '*', '*.tar.gz'))
    if dates == None:
        return sorted(fnames)
    dates = utils.daterange(dates)
    # Find matching dates
    matched = []
    for f in fnames:
        basename = os.path.basename(f)
        try:
            year = int(basename[9:13])
            doy = int(basename[13:16])
            date = datetime.date(year,1,1) + datetime.timedelta(doy-1)
            if dates[0] <= date <= dates[1]:
                matched.append(f)
        except:
            pass
    return sorted(matched)

def outfile(filename, product="radi", atmcorr=False, overwrite=False, suffix=""):
    """ Generate output filename and create needed directories """
    # Generate output filename
    basename = filename[:-12]
    if atmcorr:
        atmtag = 'A'
    else: atmtag = 'N'
    if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
    ofile = basename.replace(rawdir,outdir) \
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

def process(img, fname_out, product='radi', datatype='Int16', verbose=1, overviews=False):

    # Atmospheric correct
    #if atmcorr:
        #meta = readmeta(img.Filename())
        #atm = atmosphere(meta=meta)
        #for b in range(0,len(bandnums)):
        #band.SetAtmosphere( atm[bandi+1] )
            #if 'Xa' in atm[bandi+1]: a.Xa = atm[bandi+1]['Xa']
            #if 'Xb' in atm[bandi+1]: a.Xb = atm[bandi+1]['Xb']
            #if 'Xc' in atm[bandi+1]: a.Xc = atm[bandi+1]['Xc']

    product = product.lower()

    # Set gippy verbosity to 1 less (so a v of 1 is verbose only for python code)
    gippy.Options.SetVerbose(verbose)
    gippy.Options.SetChunkSize(64.0)

    # Copy input into output
    if product == 'radi' or product == 'refl' or product == 'temp':
        # Datatype of output
        if datatype.lower() == 'int16':
            dtype = gippy.GDT_Int16
        elif datatype.lower() == 'float32':
            dtype = gippy.GDT_Float32

        # Create output file and set default parameters
        imgout = gippy.GeoImage(fname_out, img, dtype)
        imgout.SetNoData(-32768)

        # Units
        if product == 'radi':
            units = gippy.RADIANCE
        else: units = gippy.REFLECTIVITY

        if dtype == gippy.GDT_Int16:
            if units == gippy.RADIANCE:
                imgout.SetGain(0.1)
            elif units == gippy.REFLECTIVITY:
                imgout.SetGain(0.0001)
                thermalbands = ["LWIR"]
                for tb in thermalbands:
                    try:
                        band = imgout[tb]
                        band.SetGain(0.01)
                        imgout[tb] = band
                    except:
                        pass
        else: imgout.SetGain(1.0)
        imgout = gippy.Copy(img, imgout, units)
    # TODO - this should be an algorithm ? (landsat alg instead of landsat process)
    elif product == 'cind':
        imgout = gippy.Indices(img,fname_out)
    elif product == 'temptest':
        meta = readmeta(img.Filename())
        imgout = gippy.GeoImage(fname_out,img,gippy.GDT_Int16,3)
        imgout.SetNoData(-32768)
        imgout.SetGain(0.01)
        band = img["LWIR"]
        band.ClearAtmosphere()
        gippy.Copy(band,imgout[0],gippy.REFLECTIVITY)
        band.SetAtmosphere(atmosphere(6,meta=meta))
        gippy.Copy(band,imgout[1],gippy.REFLECTIVITY)
        band.SetAtmosphere(atmosphere(6,meta=meta,merraprofile=True))
        gippy.Copy(band,imgout[2],gippy.REFLECTIVITY)
    else:
        raise Exception("Unknown product name %s" % product)

    if overviews: imgout.AddOverviews()
    # Clean up extracted tif files in original dir
    # this should be part of the destructor - add python class destructor for GeoImage
    fname = imgout.Filename()
    imgout = None
    return fname

def batchprocess(fnames, products=['radi'], atmcorr=False, 
    datatype='Int16', verbose=1, overwrite=False, suffix='', overviews=False):
    
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
            #mtlout = meta['metafilename'].replace(rawdir,outdir)
            #mtl = glob.glob(os.path.join(os.path.dirname(fin),'*MTL.txt'))[0]
            #mtlout = mtl.replace(rawdir,outdir)
            #if not os.path.exists(os.path.dirname(mtlout)):
            #    os.makedirs(os.path.dirname(mtlout))
            #try:
            #    shutil.copy(mtl,mtlout)
            #except: pass

            start = datetime.datetime.now()
            try:
                # TODO - If only doing temp then don't waste time with other bands
                img = read(fin, verbose=verbose)
                # TODO - shouldn't have to readmeta again
                meta = readmeta(fin)
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
                    process(img, fout[1], fout[0], datatype=datatype, overviews=overviews)
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
        path = os.path.join(topdir,origdir,pathrow,year+doy)
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
                    #mtl = readmeta(newf)
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

def unarchive(path):
    try:
        datafile = glob.glob(os.path.join(path,'*tar.gz'))[0]
        bname = os.path.basename(datafile)
        newf = os.path.join(qdir,bname)
        shutil.move(datafile,newf)
        print '%s: removed from archive' % os.path.basename(bname)
    except Exception,e: pass
    shutil.rmtree(path)

def _cleandir(dir):
    try:
        shutil.rmtree(os.path.join(dir,'modtran'))
    except: pass

    # Check validity
    #try:
    #    mtl = readmeta(os.path.join(dir,'dummy'))
    #except Exception,e:
    #    print 'bad file? - should unarchive it'
        #unarchive(dir)

    # Remove extraneous files
    fnames = glob.glob(os.path.join(dir,'*'))
    for f in fnames:
        extension = os.path.splitext(f)[1][1:].strip()
        if extension != 'txt' and extension != 'gz' and extension != 'tar':
            os.remove(f)

def clean():
    """ Clean out archive directories, keeping tar.gz and mtl files """
    prdirs = os.listdir(rawdir)
    print 'Cleaning Landsat archive: %s pathrows' % len(prdirs)
    for pr in prdirs:
        print 'Cleaning pathrow %s' % pr
        datedirs = os.listdir(os.path.join(rawdir,pr))
        for dd in datedirs:
            _cleandir(os.path.join(rawdir,pr,dd))
