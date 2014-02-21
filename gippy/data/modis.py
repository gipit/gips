import os
import datetime
from osgeo import gdal
from collections import OrderedDict

from gippy.data.core import Data, DataInventory
from gippy utils import File2List, List2File


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
    _pattern = 'M?D????.????????.h??v??.???.hdf'

    _prodpattern = '*.tif'
    #_metapattern = 'MTL.txt'
    _defaultresolution = [926.625433138333392, -926.625433139166944]

    _products = OrderedDict([
        ('temp', {
            'description': 'Surface temperature observations',
            'depends': ['MOD11A1', 'MYD11A1'],
        }),
        ('sds1', {
            'description': 'First SDS in the file',
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
            'filename': filename,
            'datafiles': datafiles,
            'tile': tile,
            'date': datetime.datetime.strptime(year+doy, "%Y%j").date(),
            'basename': 'MODIS_',
            'sensor': sensor,

            # optional

            'products': {'sds1': datafiles[0]}

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


    
    def fetch(self):

        datasets = set()
        for product in self.products:
            datasets.add(self._products[product]['depends'])

        for tile in self.tiles:
            for dataset in datasets:
                http_fetch(self, tile, date, dataset)



def main(): ModisData.main()
