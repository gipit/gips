#!/usr/bin/env python

import os
import datetime
import glob
import tarfile
import copy

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

    _rootdir = '/titan/data/SAR/tiles'
    #_datedir = '%Y%j'
    _tiles_vector = '/titan/data/SAR/tiles.shp'
    _pattern = 'KC_*.tar.gz'
    _prodpattern = '*.tif'
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
            VerboseOut('Date from repository = '+str(date),2)
        else:
            # extract header and date image
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
            stats = dateimg[0].ComputeStats()[0]
            date = cls._launchdate[fname[-9]] + datetime.timedelta(days=int(stats[0]))
            RemoveFiles([hdrfile,datefile,datefile+'.hdr'])
            #VerboseOut('Date from image: %s' % str(date),3) 
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
            'date': date,
            'basename': bname,
            'path': os.path.join(cls._rootdir,tile,date.strftime(cls._datedir)),
            'sensor': fname[-9:-8] + fname[-15:-12],
            # unique to SARData
            'hdrfile': hdrfile
        }

    def meta(self, tile):
        """ Get metadata for this tile """
        filename = self.tiles[tile]['filename']
        meta = self.inspect(filename)
        # add info from headerfile
        tfile = tarfile.open(meta['filename'])
        if not os.path.exists(os.path.join(meta['path'],meta['hdrfile'])):
            tfile.extract(meta['hdrfile'],path)
        meta.update( _meta(meta['hdrfile']) ) 
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

    @classmethod
    def extractdata(cls,info):
        """ Extract and create ENVI files from original data """
        index = super(SARData, cls).extractdata(filename=info['filename'])
        meta = cls._meta(index['headerfile'])
        toc = {}
        for fname in index['datafiles']:
            bandname = os.path.basename(fname)[len(info['basename'])+1:]
            envihdr = copy.deepcopy(meta['envihdr'])
            if bandname in ['mask','linci']: envihdr[6] = 'data type = 1'
            envihdr.append('band names={%s}' % bandname)
            List2File(envihdr,fname+'.hdr') 
            toc[bandname] = fname
        return toc

    #def opentile(self, tile, product=''):

    def process(self, overwrite=False, suffix=''):
        """ Make sure all files have been pre-processed """
        if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
        for tile, data in self.tiles.items():

            toc = self.extractdata(data)
            data['products']['date'] = toc['date']

            if 'sign' in self.products:
                existing_bands = []
                for b in self._databands:
                    if b in toc:  existing_bands.append(b)

                img = gippy.GeoImage(toc[existing_bands[0]])
                del existing_bands[0]
                for f in existing_bands: img.AddBand(gippy.GeoImage(toc[f])[0])

                img.SetNoData(0)
                mask = gippy.GeoImage(toc['mask'],False)
                img.AddMask(mask[0] == 255)

                # Process products
                fout = os.path.join(data['path'],data['basename']+'_sign')
                imgout = gippy.SigmaNought(img, fout)
                data['products']['sign'] = imgout.Filename()
                img = None

            # Remove unused stuff
            for k in ['linci','mask'] + self._databands:
                if k in toc: RemoveFiles([toc[k],toc[k]+'.hdr',toc[k]+'.aux.xml'])
                    

    @classmethod
    def feature2tile(cls,feature):
        """ Get tile designaation from a geospatial feature (i.e. a row) """
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
