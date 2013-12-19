#!/usr/bin/env python

import os
import datetime
import glob
from gippy.data.core import Data, main as datamain

from agspy.utils.table import Table

# Debugging
from pdb import set_trace

class CDLData(Data):
    """ A CDL (Crop Data Layer) object """
    sensors = {'cdl': 'CDL'}
    sensor = 'cdl'
    rootdir = '/titan/data/CDL'
    _products = {
        'cdl': {'description': 'Crop Data Layer'}
    }
    _tiles_vector = 'usa_states'
    _tiles_attribute = 'state_abbr'
    _legend_file = rootdir + '/CDL_Legend.csv'
    _legend = map(lambda x: x.lower(), Table(csvfile=_legend_file)['ClassName'])

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates for tile """
        files = glob.glob(os.path.join(cls.path(tile),'CDL*.tif'))
        return [datetime.datetime.strptime(os.path.basename(f)[4:8],'%Y').date() for f in files]

    @classmethod
    def find_products(cls, tile, date, products):
        """ Get filename for specified tile. Return (path,basename,filename,sensor) """
        filename = glob.glob(os.path.join(cls.rootdir, tile, 'CDL_%s_*.tif' % date.strftime('%Y')))
        path, basename = os.path.split(filename[0])
        return {'path': path, 'basename': basename, 'sensor': 'cdl', 'products': {'raw':filename[0],'cdl':filename[0]}} 

    @classmethod
    def path(cls, tile='', date=''):
        """ Path to tile directory. Date ignored since all in same directory """
        if tile == '':
            return cls.rootdir
        else:
            return os.path.join(cls.rootdir, tile)

    @classmethod
    def get_code(cls, cropname):
        ''' Retrieve CDL code for the given crop name (lower case) '''
        return cls._legend.index(cropname)

    @classmethod
    def get_cropname(cls, code):
        '''Retrieve name associated with given crop code'''
        return cls._legend[code]

def main(): datamain(CDLData)
