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

import sys
import os
import re
import time
import datetime
import urllib
from osgeo import gdal
from collections import OrderedDict

import math
import numpy as np

import gippy
from gipif.core import Repository, Asset, Data
from gipif.inventory import DataInventory
from gipif.utils import File2List, List2File, VerboseOut
import gipif.settings as settings


def binmask(arr, bit):
    """ Return boolean array indicating which elements as binary have a 1 in
        a specified bit position. Input is Numpy array.
    """
    return arr & (1 << (bit-1)) == (1 << (bit-1))


class MerraRepository(Repository):
    repo = settings.REPOS['merra']
    _rootpath = repo.get('rootpath', Repository._rootpath)
    _tiles_vector = repo.get('tiles_vector', Repository._tiles_vector)
    _tile_attribute = repo.get('tile_attribute', Repository._tile_attribute)

    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")
        h = str(int(feature.GetField(fldindex_h))).zfill(2)
        v = str(int(feature.GetField(fldindex_v))).zfill(2)
        tile = "h%sv%s" % (h, v)
        return tile


class MerraAsset(Asset):
    Repository = MerraRepository

    _sensors = {
        'MAT': {'description': 'MERRA IAU'},
        'MST': {'description': 'MERRA Land'},
        'MAI': {'description': 'MERRA DAS'}
    }

    _asset_layers = {
        'MAT1NXSLV': {}
    }

    _assets = {
        'MAT1NXSLV': {

            # MERRA300.prod.assim.tavg1_2d_slv_Nx.20100101.hdf
            'pattern': 'MERRA300.prod.assim.tavg1_2d_slv_Nx*hdf',

            'url': 'ftp://goldsmr2.sci.gsfc.nasa.gov/data/s4pa/MERRA/MAT1NXSLV.5.2.0/',
            'startdate': datetime.date(1979, 1, 11),
            'latency': -45
        },
    }

    # Not used anywhere
    _launchdate = {
    }

    # Should this be specified on a per asset basis?
    _defaultresolution = (0.666666666666667, -0.498614958448753)


    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(MerraAsset, self).__init__(filename)


        bname = os.path.basename(filename)

        self.asset = bname[0:7]
        self.tile = bname[17:23]
        year = bname[9:13]
        doy = bname[13:16]

        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()
        self.sensor = bname[:3]

        datafiles = self.datafiles()


    def datafiles(self):

        indexfile = self.filename + '.index'

        # TODO: get File2List to handle missing file and empty file in the same way
        try:
            datafiles = File2List(indexfile)
        except:
            datafiles = []

        if not datafiles:
            gdalfile = gdal.Open(self.filename)
            subdatasets = gdalfile.GetSubDatasets()
            datafiles = [s[0] for s in subdatasets]
            List2File(datafiles, indexfile)

        return datafiles

    @classmethod
    def _filepattern(cls, asset, tile, date):
        """ used by fetch only """

        year, month, day = date.timetuple()[:3]

        # MERRA300.prod.assim.tavg1_2d_slv_Nx.20100101.hdf

        # pattern = ''.join(['(', asset, '.A', str(year), str(doy).zfill(3), '.', tile, '.005.\d{13}.hdf)'])

        datestr = str(year) + str(month).zfill(2) + str(day).zfill(2)

        # pattern = ''.join(['(', 'MERRA300.prod.assim.tavg1_2d_slv_Nx.', datestr, '.hdf)'])

        pattern = "(MERRA300.prod.assim.tavg1_2d_slv_Nx.%s.hdf)" % datestr

        print "pattern", pattern

        return pattern


    @classmethod
    def _remote_subdirs(cls, asset, tile, date):
        """ used by fetch only """
        year, month, day = date.timetuple()[:3]
        httploc = cls._assets[asset]['url']

        # mainurl = ''.join([httploc, '/', str(year), '.', '%02d' % month, '.', '%02d' % day])

        mainurl = "%s/%04d/%02d" % (httploc, year, month)

        # mainurl = ''.join([httploc, '/', str(year), '/', '%02d' % month])

        print "mainurl", mainurl

        return mainurl


    @classmethod
    def fetch(cls, asset, tile, date):

        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)

        print dir(cls.Repository)

        tilesvector = cls.Repository.tiles_vector()


        # if date.date() < cls._assets[asset]['startdate']:
        #     print "date is too early"
        #     return 3

        # if date > datetime.datetime.now() - datetime.timedelta(cls._assets[asset]['latency']):
        #     print "date is too recent"
        #     return 3

        pattern = cls._filepattern(asset, tile, date)
        mainurl = cls._remote_subdirs(asset, tile, date)
        outdir = cls.Repository.spath()

        VerboseOut('%s: mainurl %s, pattern %s' % (asset, mainurl, pattern), 4)

        print "outdir", outdir
        print "pattern", pattern
        print "mainurl", mainurl


        # assetname = ""
        # tileid = ""
        
        # h = int(tileid[1:3])
        # v = int(tileid[4:6])

        # print h, v

        # x0 = ORIG[0] + 12.0*(h - 1)
        # y0 = ORIG[1] + 10.0*(v - 1)
        # dx = 12.0
        # dy = 10.0

        # nx = int(round(dx/RES[0]))
        # ny = int(round(-dy/RES[1]))
        # ix0 = int(round((x0 - ORIG[0])/RES[0]))
        # iy0 = int(round(-(y0 - ORIG[1])/RES[1]))
        # ix1 = ix0 + nx
        # iy1 = iy0 + ny

        # print "x0, y0", x0, y0
        # print "nx, ny", nx, ny
        # print "ix0, iy0", ix0, iy0
        # print "ix1, iy1", ix1, iy1

        # dataset = open_url(LOC)

        # keys = dataset.keys()

        # names = []
        # for i, key in enumerate(keys):
        #     x = dataset[key]
        #     print i, key, x.shape
        #     if x.shape == (24, 361, 540):
        #         names.append(key)        

        # print names
        # print len(names)

        # assert assetname in names, "asset name is not in SDS list"

        # data = dataset[assetname][:, iy0:iy1, ix0:ix1].astype('float32')
        # print data.shape

        # proj = raster.create_proj(4326)
        # geo = raster.create_geo((x0, y0), RES[0], RES[1])
        # meta = {}

        # print type(data)

        # print data.dtype

        # print dataset['Time'][:]
        # print dataset['TIME'][:]

        # outfilename = assetname + ".tif"
        # names = ['%02d30GMT' % i for i in range(24)]
        # raster.write_raster(outfilename, data, proj, geo, meta, bandnames=names)



 

class MerraData(Data):
    """ A tile of data (all assets and products) """
    name = 'Merra'
    Asset = MerraAsset
    _pattern = '*.tif' 
    _products = {
        'temp': {
            'description': 'Air temperature data',
            'assets': ['MAT1NXSLV']
        },
    }
    
    # def process(self, products):



def main():
    DataInventory.main(MerraData)
