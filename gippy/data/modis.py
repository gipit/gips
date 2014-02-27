import os
import re
import time
import datetime
import urllib
from osgeo import gdal
from collections import OrderedDict

from gippy.data.core import Data
from gippy.utils import File2List, List2File, VerboseOut

from pdb import set_trace


class ModisData(Data):
    """ Represent a single day and temporal extent of MODIS data along with product variations """

    name = 'Modis'

    sensors = {
        'MOD': 'Terra',
        'MYD': 'Aqua',
        'MCD': 'Combined'
    }
    _rootdir = '/titan/data/modis/tiles'
    _tiles_vector = '/titan/data/vector/MODIS/modis_sinusoidal/modis_sinusoidal_grid_world.shp'

    #_assetpattern = 'M?D????.????????.h??v??.???.hdf'

    #_pattern = 'M?D*.hdf'

    _assets = {
        'MOD11A1': {
            'pattern': 'MOD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005'
        },
        'MYD11A1': {
            'pattern': 'MYD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.005'
        }
    }

    _stage = os.path.join(_rootdir, 'stage')

    _prodpattern = '*.tif'

    _defaultresolution = [926.625433138333392, -926.625433139166944]

    _products = OrderedDict([
        ('temp', {
            'description': 'Surface temperature observations',
            # the list of asset types associated with this product
            'assets': ['MOD11A1', 'MYD11A1'],
        }),
        ('sds1', {
            'description': 'First SDS in the file',
            # the list of asset types associated with this product
            'assets': ['MOD11A1'],
        }),
    ])

    _launchdate = {
        'MOD': datetime.date(1999, 12, 1),
        'MYD': datetime.date(2002, 5, 1),
        'MCD': None
    }

    # MOD11A1.A2005241.h12v04.005.2008059145629.hdf

    @classmethod
    def inspect(cls, filename):
        """ Inspect a single file and get some metadata """
        path, basename = os.path.split(filename)

        tile = basename[17:23]
        year = basename[9:13]
        doy = basename[13:16]
        sensor = basename[:3]

        indexfile = filename + '.index'
        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:            
            gdalfile = gdal.Open(filename)
            subdatasets = gdalfile.GetSubDatasets()
            datafiles = [s[0] for s in subdatasets]
            List2File(datafiles, indexfile)

        return {
            # required
            'asset': basename[0:7],
            'filename': filename,
            'datafiles': datafiles,
            'tile': tile,
            'date': datetime.datetime.strptime(year+doy, "%Y%j").date(),
            'basename': 'MODIS'+tile+'_'+year+doy,
            'sensor': sensor,
            'products': {'sds1': datafiles[0]}
            # optional
        }

    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")
        h = str(int(feature.GetField(fldindex_h))).zfill(2)
        v = str(int(feature.GetField(fldindex_v))).zfill(2)
        tile = "h%sv%s" % (h, v)
        return tile


    @classmethod
    def fetch_asset(cls, asset, tile, date):

        print
        print "trying to fetch:"
        print asset, tile, date

        VerboseOut('about to fetch',4)

        httploc = cls._assets[asset]['url']

        year, month, day = date.timetuple()[:3]
        doy = date.timetuple()[7]

        pattern = ''.join(['(', asset, '.A', str(year), str(doy).zfill(3), '.', tile, '.005.\d{13}.hdf)'])
        mainurl = ''.join([httploc, '/', str(year), '.', '%02d'%month, '.', '%02d'%day])

        print "mainurl"
        print mainurl

        print "pattern"
        print pattern

        try:
            VerboseOut('opening listing', 1)
            listing = urllib.urlopen(mainurl).readlines()
        except Exception, e:
            listing = None
            print 'unable to access %s' % mainurl
            return 2
        
        cpattern = re.compile(pattern)
        name = None

        success = False

        for item in listing:
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                print 'found in', item.strip()
                name = cpattern.findall(item)[0]
                print 'it is', name
                url = ''.join([mainurl, '/', name])
                print 'the url is', url
                try:
                    urllib.urlretrieve(url, os.path.join(cls._stage, name))
                    print "retrieved %s" % name
                    success =  True
                except Exception, e:
                    print e
                    print 'unable to retrieve %s from %s' % (name, url)


        if not success:
            print "did not find a match for %s in listing of %s" % (pattern, mainurl)
            return 1
        else:
            return 0


def main(): ModisData.main()
