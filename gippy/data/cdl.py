#!/usr/bin/env python

import os
import glob
from datetime import datetime

from gippy.data.core import Asset, Tile, Data
from agspy.utils.table import Table

from pdb import set_trace


class CDLAsset(Asset):
    _rootpath = '/titan/data/CDL'
    _datedir = '%Y'
    _sensors = {
        'cdl': {'description': 'Crop Data Layer'}
    }
    _assets = {
        '': {
            'pattern': 'CDL*.tif'
        }
    }

    def __init__(self, filename):
        """ Inspect a CDL file """
        super(CDLAsset, self).__init__(filename)
        # TODO - get tile (state) so we can archive
        self.sensor = 'cdl'
        self.basename = self.basename[0:9]
        self.date = datetime.strptime(self.basename[4:8], self._datedir)
        self.products['cdl'] = filename

    @classmethod
    def path(cls, tile, date=''):
        return os.path.join(cls._rootpath, cls._tilesdir, tile)

    # _datedir used inappropriately elsewhere. if it wasn't this could be removed
    @classmethod
    def find_assets(cls, tile, date, asset=''):
        files = glob.glob(os.path.join(cls.path(tile), 'CDL_%s_*.tif' % date.strftime(cls._datedir)))
        found = []
        for f in files:
            found.append(cls(f))
        return found

    @classmethod
    def archive(cls, path=''):
        raise Exception('Archive not supported')


class CDLTile(Tile):
    """ A tile (CONUS State) of CDL """

    Asset = CDLAsset

    _products = {
        'cdl': {'description': 'Crop Data Layer'}
    }


class CDLData(Data):
    """ A CDL (Crop Data Layer) object """
    name = 'CDL'

    _defaultresolution = [30.0, 30.0]

    Tile = CDLTile

    _tiles_vector = 'usa_states'
    _vectoratt = 'state_name'

    _legend_file = os.path.join(CDLAsset._rootpath, 'CDL_Legend.csv')
    _legend = map(lambda x: x.lower(), Table(csvfile=_legend_file)['ClassName'])

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates for tile """
        files = glob.glob(os.path.join(CDLAsset._rootpath, CDLAsset._tilesdir, tile, 'CDL*.tif'))
        return [datetime.strptime(os.path.basename(f)[4:8], '%Y').date() for f in files]

    #def find_data(self, tile):
    #    """ Find all data/products for given tile, save in self.tiles dictionary """
    #    filename = self.find_raw(tile)
    #    if len(filename) == 0:
    #        return {}
    #    if len(filename) > 1:
    #        raise Exception('More than 1 file found for same tile/date')
    #    return self.inspect(filename[0])

    def process(self, **kwargs):
        pass

    @classmethod
    def get_code(cls, cropname):
        ''' Retrieve CDL code for the given crop name (lower case) '''
        return cls._legend.index(cropname)

    @classmethod
    def get_cropname(cls, code):
        '''Retrieve name associated with given crop code'''
        return cls._legend[code]


def main():
    CDLData.main()
