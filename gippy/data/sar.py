#!/usr/bin/env python

import os
import datetime
import glob
import tarfile
import copy
import numpy
from collections import OrderedDict
from pdb import set_trace

import gippy
from gippy.data.core import Data
from gippy.utils import VerboseOut, File2List, List2File, RemoveFiles

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
    _rootdir = '/titan/data/SAR'
    _tiledir = os.path.join(_rootdir,'tiles')
    _stagedir = os.path.join(_rootdir,'.stage')
    #_datedir = '%Y%j'

    _tiles_vector = os.path.join(_rootdir,'vectors','tiles.shp')

    _prodpattern = '*'
    _metapattern = '.hdr'

    _assets = {
        '': {
            'pattern': 'KC_*.tar.gz'
        }
    }

    _products = OrderedDict([
        ('sign', {
            'description': 'Sigma nought (radar backscatter coefficient)',
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
            RemoveFiles([hdrfile,datefile],['.hdr','.aux.xml'])
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
            'asset': '',
            'filename': filename,
            'datafiles': datafiles,
            'tile': tile, 
            'date': dates,
            'basename': bname,
            'sensor': fname[-9:-8] + fname[-15:-12],
            'products': {},
            # unique to SARData
            'hdrfile': hdrfile
        }

    def meta(self, tile):
        """ Get metadata for this tile """
        tdata = self.tiles[tile]
        # add info from headerfile
        tfile = tarfile.open(tdata['assets'][0])
        for f in tdata['datafiles'][ tdata['assets'][0] ]:
            if f[-3:] == 'hdr': hdr_bname = f
        hdr_fname = os.path.join(tdata['path'],hdr_bname)
        if not os.path.exists(hdr_fname):
            tfile.extract(hdr_bname, tdata['path'])
        #tdata.update( self._meta(hdr_fname) ) 
        return self._meta(hdr_fname)

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
        tdata = self.tiles[tile]
        meta = self.meta(tile)
        asset_name = self.tiles[tile]['assets'][0]
        if tarfile.is_tarfile(asset_name):
            tfile = tarfile.open(asset_name)
        else: raise Exception('%s is not a valid tar file' % os.path.basename(asset_name))
        tfile.extractall(tdata['path'])
        datafiles = {}
        for f in tdata['datafiles'][asset_name]:
            filename = os.path.join(tdata['path'],f)
            try:
                os.chmod(filename,0664)
            except: pass
            #if not os.path.exists(fname+'.hdr'):
            bandname = f[len(tdata['basename'])+1:]
            envihdr = copy.deepcopy(meta['envihdr'])
            if bandname in ['mask','linci']: envihdr[6] = 'data type = 1'
            envihdr.append('band names={%s}' % bandname)
            List2File(envihdr,filename+'.hdr') 
            datafiles[bandname] = filename
        return datafiles

    def processtile(self, tile, products):
        """ Make sure all products have been pre-processed """
        if len(products) == 0: raise Exception('Tile %s: No products specified' % tile)
        # extract all data from archive
        tdata = self.tiles[tile]
        datafiles = self._extract(tile)
        meta = self.meta(tile)
        if 'sign' in products.keys():
            bands = [b for b in self._databands if b in datafiles]
            img = gippy.GeoImage(datafiles[bands[0]]) 
            del bands[0]
            for b in bands: img.AddBand(gippy.GeoImage(datafiles[b])[0])
            img.SetNoData(0)
            mask = gippy.GeoImage(datafiles['mask'],False)
            img.AddMask(mask[0] == 255)
            # apply date mask
            dateimg = gippy.GeoImage(datafiles['date'],False)
            dateday = (self.date - self._launchdate[self.sensor[0]]).days
            img.AddMask(dateimg[0] == dateday)
            imgout = gippy.SigmaNought(img, products['sign'], meta['CF'])
            self.tiles[tile]['products']['sign'] = imgout.Filename()

        # Remove unused stuff
        rmfiles = [k for k in ['linci','mask','date'] + self._databands if k in tdata['products']]
        for f in rmfiles:
            RemoveFiles([self.tiles[tile]['products'][f],self.tiles[tile]['products'][f]+'.hdr',self.tiles[tile]['products'][f]+'.aux.xml'])

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
