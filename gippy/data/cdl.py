#!/usr/bin/env python

from gippy.data.core import Data, main as datamain

from agspy.utils.table import Table

class CDLData(Data):
    """ A CDL (Crop Data Layer) object """
    sensors = {'cdl': 'CDL'}
    rootdir = '/titan/data/CDL'
    _tiles_vector = 'states'
    _tiles_attribute = 'abbr'
    _legend_file = rootdir + '/CDL_Legend.csv'
    _legend = map(lambda x: x.lower(), Table(csvfile=_legend_file)['ClassName'])

    @classmethod
    def get_code(cls, cropname):
        ''' Retrieve CDL code for the given crop name (lower case) '''
        return cls._legend.index(cropname)

    @classmethod
    def get_cropname(cls, code):
        '''Retrieve name associated with given crop code'''
        return cls._legend[code]



def main(): datamain(CDLData)
