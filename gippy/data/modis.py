import os
import datetime
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
    _assetpattern = 'M?D????.????????.h??v??.???.hdf'

    _prodpattern = '*.tif'
    #_metapattern = 'MTL.txt'
    _defaultresolution = [926.625433138333392, -926.625433139166944]

    _products = OrderedDict([
        ('temp', {
            'description': 'Surface temperature observations',
            # the list of asset types associated with this product
            'depends': ['MOD11A1', 'MYD11A1'],
        }),
        ('sds1', {
            'description': 'First SDS in the file',
            # the list of asset types associated with this product
            'depends': ['MOD11A1'],
        }),
    ])

    _launchdate = {
        'MOD': datetime.date(1999, 12, 1),
        'MYD': datetime.date(2002, 5, 1),
        'MCD': None
    }

    _stage = '/titan/data/modis/stage'

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
    def fetch(cls, tile, date, products):

        VerboseOut('about to fetch',2)

        #VerboseOut(dir(self), 1)

        VerboseOut(self.products,3)
        VerboseOut(self._products,3)

        set_trace()
        VerboseOut(self.path( self.tiles.keys()[0] , self.date ),3)

        assets = set()
        for product in products:
            assets.update(self._products[product]['depends'])

        httploc = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005' # only some datasets available here

        year, month, day = date.timetuple()[:3]
        doy = date.timetuple[7]

        for asset in assets:

            pattern = ''.join(['(', asset, '.A', year, doy, '.', tile, '.005.\d{13}.hdf)'])
            mainurl = ''.join([httploc, '/', year, '.', '%02d'%month, '.', '%02d'%day])

            try:
                VerboseOut('opening listing', 1)
                # listing = urllib.urlopen(mainurl).readlines()

            except Exception, e:
                listing = None
                print 'unable to access %s' % mainurl
                continue

            cpattern = re.compile(pattern)
            name = None
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
                        urllib.urlretrieve(url, os.path.join(self._stage, name))
                        retrieved.append(name)
                    except Exception, e:
                        print 'unable to retrieve %s from %s' % (name, url)

            time.sleep(2)

def main(): ModisData.main()
