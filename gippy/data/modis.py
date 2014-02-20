
import datetime

from collections import OrderedDict

from gippy.data.core import Data, DataInventory

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

        gdalfile = gdal.Open(filename)
        subdatasets = gdalfile.GetSubDatasets()
        sds_names = [s[0] for s in subdatasets]

        return {

            'filename': filename,

            # these are just HDF sds names, should they be paths?
            'datafiles': sds_names,

            'tile': tile,

            # should this be a datetime date?
            'date': datetime.datetime.strptime(year+doy, "%Y%j"),

            # ? the file doesn't know what other files it might be used with to make products
            'basename': 'M',

            'sensor': sensor

            # what is this?
            'path': os.path.join(cls._rootdir,tile,year+doy),

        }


    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")

        h = str(feature.GetField(fldindex_h)).zfill(2)
        v = str(feature.GetField(fldindex_v)).zfill(2)

        tile = "h%sv%s" % (h, v)

        return tile


def main(): ModisData.main()
