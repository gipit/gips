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

from gippy.data.core import Repository, Asset, Data
from pdb import set_trace


class ModisAtmosRepository(Repository):
    _rootpath = '/titan/data/atmos'
    _tilesdir = ''

    @classmethod
    def path(cls, tile='', date=''):
        path = os.path.join(cls._rootpath)
        if date != '':
            path = os.path.join(path, str(date.strftime(cls._datedir)))
        return path

    @classmethod
    def find_tiles(cls):
        return []


class ModisAtmosAsset(Asset):
    Repository = ModisAtmosRepository

    _assets = {
        'MOD08': {
            'pattern': 'MOD08_D3*hdf',
            'url': ''
        },
        'MYD08': {
            'pattern': 'MYD08_D3*hdf',
            'url': ''
        }
    }


class ModisAtmosData(Data):
    name = 'Globally Gridded Atmospheric Data'
    Asset = ModisAtmosAsset

    _products = {

    }
