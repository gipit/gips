#!/usr/bin/env python

import os
import datetime
import glob
from gippy.data.core import Data

from agspy.utils.table import Table

# Debugging
from pdb import set_trace

class CDLData(Data):
    """ A CDL (Crop Data Layer) object """
    name = "CDL"
    sensors = {'cdl': 'CDL'}
    _rootdir = '/titan/data/CDL/tiles'
    _datedir = '%Y'
    _pattern = 'CDL*.tif'
    _tiles_vector = 'usa_states'
    _tiles_attribute = 'state_name'
    _products = {
        'cdl': {'description': 'Crop Data Layer'}
    }
    _legend_file = _rootdir + '/../CDL_Legend.csv'
    _legend = map(lambda x: x.lower(), Table(csvfile=_legend_file)['ClassName'])

    @classmethod
    def inspect(cls, filename):
        bname = os.path.basename(filename)
        path = os.path.dirname(filename)
        tile = os.path.basename(path)
        return {
            'tile': '', 
            'date': datetime.datetime.strptime(bname[4:8],cls._datedir),
            'basename': bname[0:9],
            'path': '',
            'sensor': 'cdl',
        }

    def find_data(self, tile):
        """ Find all data for given tile, save in self.tiles dictionary """
        filename = self.find_original(tile)
        meta = self.inspect(filename)
        products = {'raw': filename}
        meta['products'] = products
        self.tiles[tile] = meta

    def find_original(self, tile):
        filename = glob.glob(os.path.join(self._rootdir, tile, 'CDL_%s_*.tif' % self.date.strftime('%Y')))
        if len(filename) == 0:
            raise Exception('No data for this tile/date')
        elif len(filename) > 1:
            raise Exception('More than 1 file found for same tile/date')
        return filename[0]

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates for tile """
        files = glob.glob(os.path.join(cls._rootdir,tile,'CDL*.tif'))
        return [datetime.datetime.strptime(os.path.basename(f)[4:8],'%Y').date() for f in files]

    @classmethod
    def archive(cls, path=''):
        raise Exception('Archive not supported')

    @classmethod
    def get_code(cls, cropname):
        ''' Retrieve CDL code for the given crop name (lower case) '''
        return cls._legend.index(cropname)

    @classmethod
    def get_cropname(cls, code):
        '''Retrieve name associated with given crop code'''
        return cls._legend[code]

def main(): CDLData.main()
