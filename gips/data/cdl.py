#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  matt.a.hanson@gmail.com
#
#    Copyright (C) 2014 Applied Geosolutions
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
from datetime import datetime
from csv import DictReader

from gips.data.core import Repository, Asset, Data


class CDLRepository(Repository):
    name = 'CDL'
    description = 'Crop Data Layer'
    _datedir = '%Y'
    _defaultresolution = [30.0, 30.0]


class CDLAsset(Asset):
    Repository = CDLRepository
    _sensors = {
        'cdl': {'description': 'Crop Data Layer'}
    }
    _assets = {
        '': {
            'pattern': '*.tif'
        }
    }

    def __init__(self, filename):
        """ Inspect a CDL file """
        super(CDLAsset, self).__init__(filename)
        # TODO - get tile (state) so we can archive
        bname = os.path.basename(filename)
        try:
            self.date = datetime.strptime(bname[4:8], self.Repository._datedir)
        except:
            self.date = datetime.strptime(bname[13:17], self.Repository._datedir)
        self.products['cdl'] = filename
        self.sensor = 'cdl'

    @classmethod
    def archive(cls, path=''):
        raise Exception('Archive not supported')


class CDLData(Data):
    """ A tile (CONUS State) of CDL """
    name = 'CDL'
    version = '0.9.0'
    Asset = CDLAsset
    _products = {
        'cdl': {'description': 'Crop Data Layer'}
    }

    _legend_file = os.path.join(CDLRepository.rootpath(), 'CDL_Legend.csv')
    _legend = [row['ClassName'].lower() for row in DictReader(open(_legend_file))]

    @classmethod
    def get_code(cls, cropname):
        ''' Retrieve CDL code for the given crop name (lower case) '''
        return cls._legend.index(cropname)

    @classmethod
    def get_cropname(cls, code):
        '''Retrieve name associated with given crop code'''
        return cls._legend[code]
