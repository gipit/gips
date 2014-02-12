#!/usr/bin/env python

import os
import datetime
import glob
import tarfile

from gippy.data.core import Data, VerboseOut, File2List
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
    _datedir = '%Y'
    _tiles_vector = '/titan/data/SAR/tiles.shp'
    _pattern = 'KC_*.tar.gz'
    _prodpattern = '*.tif'
    _metapattern = '.hdr'
    _products = OrderedDict([
        ('sign', {
            'description': 'Sigma nought (radar backscatter coefficient)',
        }),
    ])

    @classmethod
    def inspect(cls, filename):
        """ Inspect a single file and get some metadata """
        path, basename = os.path.split(filename)
        # extract metadata file
        meta = File2List( os.path.join(path,cls.extracthdr(filename)) )
        tile = basename[10:17]
        datestr = meta[2].zfill(4)
        if datestr == '0000': datestr = '1996'
        date = datetime.datetime.strptime(datestr, '%Y')

        tfile = tarfile.open(filename)
        filenames = tfile.getnames()
        for f in filenames: 
            if f[-4:] == 'date': bname = f[:-5]

        return {
            'tile': tile, 
            'basename': bname,
            'sensor': basename[-9:-8] + basename[-15:-12],
            'path': os.path.join(cls._rootdir,tile,date.strftime('%Y')),
            'res': float(meta[7]),
            'CF': float(meta[21])
        }

    @classmethod
    def feature2tile(cls,feature):
        """ Get tile designaation from a geospatial feature (i.e. a row) """
        fldindex_lat = feature.GetFieldIndex("lat")
        fldindex_lon = feature.GetFieldIndex("lon")
        lat = abs(int(feature.GetField(fldindex_lat)+0.5))
        lon = abs(int(feature.GetField(fldindex_lon)-0.5))
        if lat < 0:
            lat_h = 'S'
        else: lat_h = 'N'
        if lon < 0:
            lon_h = 'S'
        else: lon_h = 'N'
        tile = lat_h + str(lat).zfill(2) + lon_h + str(lon).zfill(3)
        return tile

    @classmethod
    def archive(cls, path=''):
        super(SARData, cls).archive(path=path)
        # remove leftover header files
        hdrfiles = glob.glob( os.path.join(path,'*'+cls._metapattern) )
        for f in hdrfiles: os.remove(f)

    def process(self, overwrite=False, suffix=''):
        """ Make sure all files have been pre-processed """
        if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
        for tile, data in self.tiles.items():

            # Create readable ENVI file from raw originals
            hdrfile = self.extracthdr(data['products']['raw'])
            hdr = File2List(hdrfile)
            datafiles = self.extract(data['products']['raw'])
            proj = (
                'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984", SPHEROID["WGS_1984",6378137.0,298.257223563]],' + 
                'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
            size = [int(hdr[23]), int(hdr[24])]
            lat = [ float(hdr[12]), float(hdr[14]) ]
            lat = [ min(lat), max(lat) ]
            lon = [ float(hdr[13]), float(hdr[15]) ]
            lon = [ min(lon), max(lon)]
            res = [ (lon[1]-lon[0])/(size[0]-1), (lat[1]-lat[0])/(size[1]-1) ]
            envihdr = ['ENVI','samples = %s' % size[0], 'lines = %s' % size[1],
                'bands = 1','header offset = 0','file type = ENVI Standard','data type = 2',
                'interleave = bsq', 'sensor type = Unknown', 'byte order = 0',
                'coordinate system string = ' + proj,
                'data ignore value = 0',
                'map info = {Geographic Lat/Lon, 1, 1, %s, %s, %s, %s}' % (lon[0], lat[1], res[0], res[1]) ]
            bandfiles = []
            bands = ["HH","HV"]
            for fname in datafiles:
                f = open(fname+'.hdr','w')
                f.write('\n'.join(envihdr)+'\n')
                f.write('band names = {%s}' % fname[len(os.path.join(data['path'],data['basename']))+1:])
                f.close()
                if fname[-2:] in bands: bandfiles.append(os.path.join(data['path'],fname))

            # Convert to product - sigma naut TODO - check for existence
            img = gippy.GeoImage(bandfiles[0],False)
            del bandfiles[0]
            for f in bandfiles: img.AddBand(gippy.GeoImage(f,False)[0])
            imgout = gippy.SigmaNought(img, os.path.join(data['path'],data['basename']+'_sign') )
            data['products']['sign'] = imgout.Filename()

            return
            # Clean up
            for df in datafiles:
                files = glob.glob(df+'*')
                for f in files:
                    try:
                        os.remove(f)
                    except: pass

def main(): SARData.main()