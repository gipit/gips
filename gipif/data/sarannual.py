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
import glob
import tarfile
import copy
import numpy
from collections import OrderedDict

import gippy
from gipif.core import Repository, Asset, Data
from gipif.inventory import DataInventory
from gipif.utils import VerboseOut, File2List, List2File, RemoveFiles, Colors
import gipif.settings as settings


class SARAnnualRepository(Repository):
    repo = settings.REPOS['SARannual']
    _rootpath = repo.get('rootpath', Repository._rootpath)
    _tiles_vector = repo.get('tiles_vector', Repository._tiles_vector)
    _tile_attribute = repo.get('tile_attribute', Repository._tile_attribute)
    _datedir = '%Y'

    @classmethod
    def feature2tile(cls, feature):
        """ Get tile designation from a geospatial feature (i.e. a row) """
        fldindex_lat = feature.GetFieldIndex("lat")
        fldindex_lon = feature.GetFieldIndex("lon")
        lat = int(feature.GetField(fldindex_lat)+0.5)
        lon = int(feature.GetField(fldindex_lon)-0.5)
        if lat < 0:
            lat_h = 'S'
        else:
            lat_h = 'N'
        if lon < 0:
            lon_h = 'W'
        else:
            lon_h = 'E'
        tile = lat_h + str(abs(lat)).zfill(2) + lon_h + str(abs(lon)).zfill(3)
        return tile


class SARAnnualAsset(Asset):
    Repository = SARAnnualRepository
    _sensors = {
        #'AFBS': 'PALSAR FineBeam Single Polarization',
        'PALSAR': {'description': 'PALSAR Mosaic (FineBeam Dual Polarization)', 'colors': Colors.BOLD + Colors.RED},
        #'AWB1': 'PALSAR WideBeam (ScanSAR Short Mode)',
        #'JFBS': 'JERS-1 FineBeam Single Polarization'
    }
    _assets = {
        'MOS': {
            'pattern': '???????_??_MOS.tar.gz'
        },
        'FNF': {
            'pattern': '???????_??_FNF.tar.gz'
        },
    }

    _defaultresolution = [0.00044444444, 0.00044444444]

    def __init__(self, filename):
        """ Inspect a single file and get some basic info """
        super(SARAnnualAsset, self).__init__(filename)
        bname = os.path.basename(filename)
        self.asset = bname[11:14]
        self.tile = bname[0:7]
        self.sensor = 'PALSAR'
        self.date = datetime.datetime.strptime(bname[8:10], '%y')
        self.rootname = bname[0:10]

    def extract(self, filenames=[]):
        """ Extract filesnames from asset """
        files = super(SARAnnualAsset, self).extract(filenames)
        datafiles = {}
        for f in files:
            bname = os.path.basename(f)
            if f[-3:] != 'hdr':
                bandname = bname[len(self.rootname)+1:]
                datafiles[bandname] = f
        return datafiles


class SARAnnualData(Data):
    """ Tile of data """
    name = 'SARAnnual'
    Asset = SARAnnualAsset

    _pattern = '*'
    _products = OrderedDict([
        ('sign', {
            'description': 'Sigma nought (radar backscatter coefficient)',
            'assets': 'MOS',
        }),
        ('fnf', {
            'description': 'Forest/NonForest Mask',
            'assets': 'FNF',
        })
    ])
    _groups = {
        'Standard': _products.keys(),
    }

    def meta(self, tile):
        """ Get metadata for this tile """
        return {'CF': -83.0}

    def process(self, products):
        """ Process all requested products for this tile """
        if len(products) == 0:
            raise Exception('Tile %s: No products specified' % tile)

        for key, val in products.items():
            fname = os.path.join(self.path, self.basename + '_' + key)
            if val[0] == 'sign':
                datafiles = self.assets['MOS'].extract()
                bands = [datafiles[b] for b in ["sl_HH", "sl_HV"] if b in datafiles]
                if len(bands) > 0:
                    img = gippy.GeoImage(bands)
                    img.SetNoData(0)
                    mask = gippy.GeoImage(datafiles['mask'], False)
                    img.AddMask(mask[0] == 255)
                    imgout = gippy.GeoImage(fname, img, gippy.GDT_Float32)
                    imgout.SetNoData(-32768)
                    for b in range(0, imgout.NumBands()):
                        imgout.SetColor(img[b].Description(), b+1)
                        (img[b].pow(2).log10() * 10 - 83.0).Process(imgout[b])
                    fname = imgout.Filename()
                    img = None
                    imgout = None
                    for key, f in datafiles.items():
                        if key != 'hdr':
                            RemoveFiles([f], ['.hdr', '.aux.xml'])
            if val[0] == 'fnf':
                datafiles = self.assets['FNF'].extract()
                if 'C' in datafiles:
                    # rename both files to product name
                    os.rename(datafiles['C'], fname)
                    os.rename(datafiles['C']+'.hdr', fname+'.hdr')
                    img = gippy.GeoImage(fname)
                    img.SetNoData(0)
                    img = None
            self.products[key] = fname


def main():
    DataInventory.main(SARAnnualData)
