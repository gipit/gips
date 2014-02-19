#!/usr/bin/env python

import os
import datetime
import glob
import tarfile
import copy
import numpy

# move VerboseOut, File2List to utils or main gippy
from gippy.data.core import Data, VerboseOut, File2List, List2File, RemoveFiles
import gippy
from pdb import set_trace

from collections import OrderedDict

class SARData(Data):
    """ Represents a single date and temporal extent along with (existing) product variations """
    name = 'SAR'
    sensors = {
        'AFBS':'PALSAR FineBeam Single Polarization',
        'AFBD':'PALSAR FineBeam Dual Polarization',
        'AWB1':'PALSAR WideBeam (ScanSAR Short Mode)',
        'JFBS':'JERS-1 FineBeam Single Polarization'
    }
    _defaultresolution = [0.000834028356964,0.000834028356964]
    _rootdir = '/titan/data/SAR/tiles.dev'
    #_datedir = '%Y%j'
    _tiles_vector = '/titan/data/SAR/tiles.shp'
    _pattern = 'KC_*.tar.gz'
    _prodpattern = '*'
    _metapattern = '.hdr'
    _products = OrderedDict([
        ('sign', {
            'description': 'Sigma nought (radar backscatter coefficient)',
        }),
        ('date', {
            'description': 'Pixel observation dates (days since launch)',
        }),
    ])

    # SAR specific constants
    # launch dates for PALSAR (A) and JERS-1 (J)
    _launchdate = {'A': datetime.date(2006,1,24), 'J': datetime.date(1992,2,11)}
    _databands = ["sl_HH","sl_HV"]

    _cycledates = {
        7:  '20-Oct-06',
        8:  '05-Dec-06',
        9:  '20-Jan-07',
        10: '07-Mar-07',
        11: '22-Apr-07',
        12: '07-Jun-07',
        13: '23-Jul-07',
        14: '07-Sep-07',
        15: '23-Oct-07',
        16: '08-Dec-07',
        17: '23-Jan-08',
        18: '09-Mar-08',
        19: '24-Apr-08',
        20: '09-Jun-08',
        21: '25-Jul-08',
        22: '09-Sep-08',
        23: '25-Oct-08',
        24: '10-Dec-08',
        25: '25-Jan-09',
        26: '12-Mar-09',
        27: '27-Apr-09',
        28: '12-Jun-09',
        29: '28-Jul-09',
        30: '12-Sep-09',
        31: '28-Oct-09',
        32: '13-Dec-09',
        33: '28-Jan-10',
        34: '15-Mar-10',
        35: '30-Apr-10',
        36: '15-Jun-10',
        37: '31-Jul-10',
        38: '15-Sep-10',
        39: '31-Oct-10',
        40: '16-Dec-10',
        41: '31-Jan-11',
        42: '18-Mar-11',
        43: '03-May-11'
    }

    @classmethod
    def inspect(cls, filename):
        """ Inspect a single file and get some basic info """
        path, fname = os.path.split(filename)
        tile = fname[10:17]
        #start = datetime.datetime.now()
        # read date
        indexfile = os.path.join(path,fname+'.index')
        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:            
            tfile = tarfile.open(filename)
            datafiles = tfile.getnames()
            List2File(datafiles,indexfile)
        for f in datafiles: 
            if f[-3:] == 'hdr': 
                hdrfile = f
            if f[-4:] == 'date': 
                datefile = f
                bname = f[:-5]

        # Check if inspecting a file in the repository
        if cls._rootdir in path:
            date = datetime.datetime.strptime(os.path.basename(path),cls._datedir).date()
            dates = date
            #VerboseOut('Date from repository = '+str(dates),4)
        else:
            # extract header and date image
            tfile = tarfile.open(filename)
            tfile.extract(hdrfile,path)
            hdrfile = os.path.join(path,hdrfile)
            os.chmod(hdrfile,0664)
            meta = cls._meta(hdrfile)
            tfile.extract(datefile,path)
            os.chmod(os.path.join(path,datefile),0664)
            # Write ENVI header for date image
            List2File(meta['envihdr'],datefile+'.hdr')
            dateimg = gippy.GeoImage(datefile)
            dateimg.SetNoData(0)
            dates = [cls._launchdate[fname[-9]] + datetime.timedelta(days=int(d)) for d in numpy.unique(dateimg.Read()) if d != 0]
            if not dates: raise Exception('%s: no valid dates' % fname)
            date = min(dates)
            #stats = dateimg[0].ComputeStats()[0]
            #date = cls._launchdate[fname[-9]] + datetime.timedelta(days=int(stats[0]))
            RemoveFiles([hdrfile,datefile,datefile+'.hdr'])
            #VerboseOut('Date from image: %s' % str(date),3) 
            # If year provided check
            if fname[7] == 'Y' and fname[8:10] != '00':
                ydate = datetime.datetime.strptime(fname[8:10],'%y')
                if date.year != ydate.year:
                    raise Exception('%s: Date %s outside of expected year (%s)' % (fname, str(date),str(ydate)))
            # If widebeam check cycle dates
            if fname[7] == 'C':
                cdate = datetime.datetime.strptime(cls._cycledates[int(fname[8:10])],'%d-%b-%y').date()
                if not (cdate <= date <= (cdate + datetime.timedelta(days=45))):
                    raise Exception('%s: Date %s outside of cycle range (%s)' % (fname, str(date),str(cdate)))
            #VerboseOut('%s: inspect %s' % (fname,datetime.datetime.now()-start), 4)

        return {
            'filename': filename,
            'datafiles': datafiles,
            'tile': tile, 
            'date': dates,
            'basename': bname,
            'sensor': fname[-9:-8] + fname[-15:-12],
            # unique to SARData
            'hdrfile': hdrfile
        }

    def meta(self, tile):
        """ Get metadata for this tile """
        #meta = self.inspect(filename)
        meta = self.tiles[tile]
        # add info from headerfile
        tfile = tarfile.open(meta['filename'])
        path = self.path(meta['tile'],meta['date'])
        hdrfile = os.path.join(path,meta['hdrfile'])
        if not os.path.exists(hdrfile):
            tfile.extract(meta['hdrfile'],path)
        meta.update( self._meta(hdrfile) ) 
        return meta

    @classmethod
    def _meta(cls,hdrfile):
        """ Get some metadata from header file """
        hdr = File2List( hdrfile )
        meta = {}
        meta['proj'] = (
            'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984", SPHEROID["WGS_1984",6378137.0,298.257223563]],' + 
            'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
        meta['size'] = [int(hdr[23]), int(hdr[24])]
        lat = [ float(hdr[12]), float(hdr[14]) ]
        lat = [ min(lat), max(lat) ]
        meta['lat'] = lat
        lon = [ float(hdr[13]), float(hdr[15]) ]
        lon = [ min(lon), max(lon)]
        meta['lon'] = lon
        meta['res'] = [ (lon[1]-lon[0])/(meta['size'][0]-1), (lat[1]-lat[0])/(meta['size'][1]-1) ]
        meta['envihdr'] = ['ENVI','samples = %s' % meta['size'][0], 'lines = %s' % meta['size'][1],
            'bands = 1','header offset = 0','file type = ENVI Standard','data type = 12',
            'interleave = bsq', 'sensor type = Unknown', 'byte order = 0',
            'coordinate system string = ' + meta['proj'],
            'data ignore value = 0',
            'map info = {Geographic Lat/Lon, 1, 1, %s, %s, %s, %s}' % (lon[0], lat[1], meta['res'][0], meta['res'][1]) ]
        meta['CF'] = float(hdr[21])
        return meta

    def _extract(self,tile):
        meta = self.meta(tile)
        path,bname = os.path.split(meta['filename'])
        if tarfile.is_tarfile(meta['filename']):
            tfile = tarfile.open(meta['filename'])
        else: raise Exceptin('%s is not a valid tar file' % bname)
        tfile.extractall(path)
        prods = {}
        for f in self.tiles[tile]['datafiles']:
            filename = os.path.join(path,f)
            try:
                os.chmod(filename,0664)
            except: pass
            #if not os.path.exists(fname+'.hdr'):
            bandname = f[len(meta['basename'])+1:]
            envihdr = copy.deepcopy(meta['envihdr'])
            if bandname in ['mask','linci']: envihdr[6] = 'data type = 1'
            envihdr.append('band names={%s}' % bandname)
            List2File(envihdr,filename+'.hdr') 
            prods[bandname] = filename
        # update with available products
        self.tiles[tile]['products'].update(prods)

    def processtile(self, tile, products):
        """ Make sure all products have been pre-processed """
        # is this needed? or is processtile not called if products empty
        if len(products) == 0: return

        # extract all data from archive
        self._extract(tile)
        meta = self.meta(tile)

        if 'sign' in products.keys():
            avail = self.tiles[tile]['products']
            bands = [b for b in self._databands if b in avail]
            print avail
            img = gippy.GeoImage(avail[bands[0]]) 
            del bands[0]
            for b in bands: img.AddBand(gippy.GeoImage(avail[b])[0])
            img.SetNoData(0)
            mask = gippy.GeoImage(avail['mask'],False)
            img.AddMask(mask[0] == 255)
            # apply date mask
            dateimg = gippy.GeoImage(avail['date'],False)
            dateday = (meta['date'] - self._launchdate[meta['sensor'][0]]).days
            img.AddMask(dateimg[0] == dateday)
            imgout = gippy.SigmaNought(img, products['sign'], meta['CF'])
            avail['sign'] = imgout.Filename()
            print self.tiles[tile]

        # Remove unused stuff - always leave date product
        for k in ['linci','mask'] + self._databands:
            if k in self.tiles[tile]['products']:
                #RemoveFiles([self.tiles[tile]['products'][k],self.tiles[tile]['products'][k]+'.hdr',self.tiles[tile]['products'][k]+'.aux.xml'])
                del self.tiles[tile]['products'][k]

    @classmethod
    def feature2tile(cls,feature):
        """ Get tile designation from a geospatial feature (i.e. a row) """
        fldindex_lat = feature.GetFieldIndex("lat")
        fldindex_lon = feature.GetFieldIndex("lon")
        lat = int(feature.GetField(fldindex_lat)+0.5)
        lon = int(feature.GetField(fldindex_lon)-0.5)
        if lat < 0:
            lat_h = 'S'
        else: lat_h = 'N'
        if lon < 0:
            lon_h = 'W'
        else: lon_h = 'E'
        tile = lat_h + str(abs(lat)).zfill(2) + lon_h + str(abs(lon)).zfill(3)
        return tile

def main(): SARData.main()
