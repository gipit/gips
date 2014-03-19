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
from pdb import set_trace

import gippy
from gippy.data.core import Repository, Asset, Tile, Data
from gippy.utils import VerboseOut, File2List, List2File, RemoveFiles


class SARAnnualRepository(Repository):
    _rootpath = '/titan/data/SARannual'
    _datedir = '%Y'
    #_tilesdir = 'tiles.dev'

    _tiles_vector = 'sar_tiles'

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
        'PALSAR': {'description': 'PALSAR Mosaic (FineBeam Dual Polarization)'},
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
        self.asset = self.basename[11:14]
        self.tile = self.basename[0:7]
        self.sensor = 'PALSAR'
        self.date = datetime.datetime.strptime(self.basename[8:10], '%y')
        self.basename = self.basename[0:10]

    def extract(self, filenames=[]):
        """ Extract filesnames from asset """
        files = super(SARAnnualAsset, self).extract(filenames)
        datafiles = {}
        for f in files:
            bname = os.path.basename(f)
            if f[-3:] != 'hdr':
                bandname = bname[len(self.basename)+1:]
                datafiles[bandname] = f
        return datafiles


class SARAnnualTile(Tile):
    """ Tile of data """

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

        if 'sign' in products.keys():
            datafiles = self.assets['MOS'].extract()
            bands = [b for b in ["sl_HH", "sl_HV"] if b in datafiles]
            print datafiles
            print bands
            if len(bands) > 0:
                img = gippy.GeoImage(datafiles[bands[0]])
                del bands[0]
                for b in bands:
                    img.AddBand(gippy.GeoImage(datafiles[b])[0])
                img.SetNoData(0)
                mask = gippy.GeoImage(datafiles['mask'], False)
                img.AddMask(mask[0] == 255)
                imgout = gippy.SigmaNought(img, products['sign'][0], -83.0)
                self.products['sign'] = imgout.Filename()
                img = None
                imgout = None
                for key, f in datafiles.items():
                    if key != 'hdr':
                        RemoveFiles([f], ['.hdr', '.aux.xml'])

        if 'fnf' in products.keys():
            datafiles = self.assets['FNF'].extract()
            if 'C' in datafiles:
                # rename both files to product name
                newfilename = datafiles['C'][:-1]+'fnf'
                os.rename(datafiles['C'], newfilename)
                os.rename(datafiles['C']+'.hdr', newfilename+'.hdr')
                img = gippy.GeoImage(newfilename)
                img.SetNoData(0)
                self.products['fnf'] = newfilename
                img = None


class SARAnnualData(Data):
    """ Represents a single date and temporal extent along with (existing) product variations """

    Tile = SARAnnualTile


def main():
    SARAnnualData.main()
