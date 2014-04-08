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
import datetime
import gdal

from gippy.data.inventory import DataInventory
from gippy.data.core import Repository, Asset, Data
from gippy.utils import File2List, List2File, VerboseOut

from pdb import set_trace


class ModisAtmosRepository(Repository):
    _rootpath = '/titan/data/modisatmos'
    _datedir = '%Y%j'

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

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(ModisAtmosAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        self.asset = bname[0:5]
        self.tile = ''
        year = bname[10:14]
        doy = bname[14:17]
        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()
        self.sensor = bname[:3]
        datafiles = self.datafiles()
        self.products = {}

    dev

    def datafiles(self):
        indexfile = self.filename + '.index'

        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:
            gdalfile = gdal.Open(self.filename)
            subdatasets = gdalfile.GetSubDatasets()
            datafiles = [s[0] for s in subdatasets]
            List2File(datafiles, indexfile)
        return datafiles


class ModisAtmosData(Data):
    name = 'Globally Gridded Atmospheric Data'
    Asset = ModisAtmosAsset

    _products = {

    }


def main():
    DataInventory.main(ModisAtmosData)
