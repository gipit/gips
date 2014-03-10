#!/usr/bin/env python

import os
import datetime
import glob
import tarfile
import copy
import numpy
from collections import OrderedDict
from pdb import set_trace

import gippy
from gippy.data.core import Asset, Tile, Data
from gippy.utils import VerboseOut, File2List, List2File, RemoveFiles


class SARAnnualAsset(Asset):
    _rootpath = '/titan/data/SARannual'
    _datedir = '%Y'
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
    _prodpattern = '*'
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

    Asset = SARAnnualAsset

    def meta(self, tile):
        """ Get metadata for this tile """
        return {'CF': -83.0}

    def process(self, products):
        """ Make sure all products have been pre-processed """
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
                imgout = gippy.SigmaNought(img, products['sign'], -83.0)
                self.products['sign'] = imgout.Filename()
                img = None
                imgout = None
                for key, f in datafiles.items():
                    if key != 'hdr':
                        RemoveFiles([f], ['.hdr', '.aux.xml'])

        if 'fnf' in products.keys():
            datafiles = self.assets['FNF'].extract()
            set_trace()
            if 'C' in datafiles:
                self.products['fnf'] = datafiles['C']


class SARAnnualData(Data):
    """ Represents a single date and temporal extent along with (existing) product variations """
    name = 'SARannual'

    _defaultresolution = [0.00044444444, 0.00044444444]

    Tile = SARAnnualTile

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


def main():
    SARAnnualData.main()
