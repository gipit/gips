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
import re
import datetime
from csv import DictReader
from zipfile import ZipFile

from dbfread import DBF

from gips.data.core import Repository, Asset, Data
from gippy import GeoImage
from osgeo import gdal

import imghdr

# make the compiler spell-check the one sensor, product, and asset type in the driver
_cdl = 'cdl'
_cdlmkii = 'cdlmkii'

class cdlRepository(Repository):
    name = 'CDL'
    description = 'Crop Data Layer'
    _datedir = '%Y'
    _defaultresolution = [30.0, 30.0]
    _tile_attribute = 'STATE_ABBR'


class cdlAsset(Asset):
    Repository = cdlRepository

    _sensors = {_cdl: {'description': 'Crop Data Layer'}}
    _assets = {
        _cdl: {
            # CDL assets are named just like products: tile_date_sensor_asset-product.tif
            'pattern': r'^(?P<tile>[A-Z]{2})_(?P<date>\d{4})_' + _cdl + '_' + _cdl + '\.tif$'
        },
        _cdlmkii: {
            'pattern': r'^(?P<tile>[A-Z]{2})_(?P<date>\d{4})_' + _cdl + '_' + _cdlmkii + '\.zip$',
            'description': '',
        },
    }

    def __init__(self, filename):
        """Use the given filename to set metadata."""
        super(cdlAsset, self).__init__(filename)
        self.tile, date_str = self.parse_asset_fp().group('tile', 'date')
        self.date = datetime.datetime.strptime(date_str, self.Repository._datedir).date()
        if re.match(self._assets[_cdl]['pattern'], os.path.basename(filename)):
            self.asset = _cdl
        else:
            self.asset = _cdlmkii
        self.sensor = _cdl
        self.products[self.asset] = filename # magically it is also a product


class cdlData(Data):
    """ A tile (CONUS State) of CDL """
    name = 'CDL'
    version = '0.9.0'
    Asset = cdlAsset
    _products = {
        _cdl: {
            'description': 'Crop Data Layer',
            'assets': [_cdl, _cdlmkii],
            'bands': [{'name': _cdl, 'units': 'none'}],
            # presently 'startdate' & 'latency' are permitted to be unspecified
            # by DH gips.utils.get_data_variables
        }
    }

    @Data.proc_temp_dir_manager
    def process(self, products, overwrite=False, **kwargs):
        for asset_type, asset in self.assets.iteritems():
            if asset_type != _cdlmkii:  # with older cdl products, the asset is the product
                continue

            fname = self.temp_product_filename(_cdl, _cdlmkii)
            fname_without_ext, _ = os.path.splitext(fname)

            with ZipFile(asset.filename, 'r') as zipfile:
                for member in zipfile.infolist():
                    member_ext = member.filename.split('.', 1)[1]
                    extracted = zipfile.extract(member, fname_without_ext)
                    os.rename(extracted, fname_without_ext + '.' + member_ext)


            image = GeoImage(fname, True)
            image[0].SetNoData(0)
            image = None

            image = gdal.Open(fname, gdal.GA_Update)
            dbf = DBF(fname + '.vat.dbf')
            for i, record in enumerate(dbf):
                image.SetMetadataItem(str("CLASS_NAME_%s" % record['CLASS_NAME']), str(i))
            image = None

            archive_fp = self.archive_temp_path(fname)
            self.AddFile(_cdl, _cdl, archive_fp)

    def legend(self):
        """Open the legend file, keeping it memoized for future calls."""
        if getattr(self, "_legend", None) is None:
            if self.assets.keys()[0] == _cdlmkii:
                self._legend = [''] * 256
                im = gdal.Open(os.path.splitext(self.assets[_cdlmkii].filename)[0] + '.tif')
                for key, val in im.GetMetadata().iteritems():
                    if key[0:10] == 'CLASS_NAME':
                        self._legend[int(val)] = key[11:]
            else:
                legend_fp = os.path.join(cdlRepository.get_setting('repository'), 'CDL_Legend.csv')
                self._legend = [row['ClassName'].lower() for row in DictReader(open(legend_fp))]
        return self._legend

    def get_code(self, cropname):
        """Retrieve CDL code for the given crop name (lower case)."""
        return self.legend().index(cropname)

    def get_cropname(self, code):
        """Retrieve name associated with given crop code."""
        return self.legend()[code]
