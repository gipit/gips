#!/usr/bin/env python
################################################################################
#    GIPPY: Geospatial Image Processing library for Python
#
#    Copyright (C) 2014 Matthew A Hanson
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>
################################################################################

import os
import glob
from datetime import datetime

from gips.core import Repository, Asset, Data
from gips.inventory import DataInventory
from agspy.utils.table import Table
import gips.settings as settings


class CDLRepository(Repository):
    repo = settings.REPOS['CDL']
    _rootpath = repo.get('rootpath', Repository._rootpath)
    _tiles_vector = repo.get('tiles_vector', Repository._tiles_vector)
    _tile_attribute = repo.get('tile_attribute', Repository._tile_attribute)
    _datedir = '%Y'

    _defaultresolution = [30.0, 30.0]

    @classmethod
    def path(cls, tile, date=''):
        return os.path.join(cls._rootpath, cls._tilesdir, tile)

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates for tile """
        files = glob.glob(os.path.join(CDLRepository.path(tile), 'CDL*.tif'))
        return [datetime.strptime(os.path.basename(f)[4:8], '%Y').date() for f in files]


class CDLAsset(Asset):
    Repository = CDLRepository
    _sensors = {
        '': {'description': 'Crop Data Layer'}
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
        bname = os.path.basename(filename)
        self.date = datetime.strptime(bname[4:8], self.Repository._datedir)
        self.products['cdl'] = filename

    # _datedir used inappropriately elsewhere. if it wasn't this could be removed
    @classmethod
    def discover(cls, tile, date, asset=''):
        pattern = os.path.join(cls.Repository.path(tile), 'CDL_%s_*.tif' % date.strftime(cls.Repository._datedir))
        files = glob.glob(pattern)
        found = []
        for f in files:
            found.append(cls(f))
        return found

    @classmethod
    def archive(cls, path=''):
        raise Exception('Archive not supported')


class CDLData(Data):
    """ A tile (CONUS State) of CDL """
    name = 'CDL'
    Asset = CDLAsset
    _products = {
        'cdl': {'description': 'Crop Data Layer'}
    }

    _legend_file = os.path.join(CDLRepository._rootpath, 'CDL_Legend.csv')
    _legend = map(lambda x: x.lower(), Table(csvfile=_legend_file)['ClassName'])

    @classmethod
    def get_code(cls, cropname):
        ''' Retrieve CDL code for the given crop name (lower case) '''
        return cls._legend.index(cropname)

    @classmethod
    def get_cropname(cls, code):
        '''Retrieve name associated with given crop code'''
        return cls._legend[code]


def main():
    DataInventory.main(CDLData)
