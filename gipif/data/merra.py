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
        # there is not a MERRA sensor. These are product groups
        'MAT1NXSLV': {'description': 'MERRA IAU 2d atmospheric single-level diagnostics'},
        'MST1NXMLD': {'description': 'MERRA-Land 2d land surface diagnostics'},
        'MAI1NXINT': {'description': 'MERRA IAU 2d Vertical integrals'}
    }

    # _asset_layers = {
    #     'MAT1NXSLV': {}
    # }

    _assets = {
        'TS': {
            'description': 'Surface skin temperature',
            'pattern': 'MERRA_MAT1NXSLV_TS_*.tif',
            'url': 'http://goldsmr2.sci.gsfc.nasa.gov:80/opendap/MERRA/MAT1NXSLV.5.2.0',
            'source': 'MERRA%s.prod.assim.tavg1_2d_slv_Nx.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 11),
            'latency': -45
        },
        'T2M': {
            'description': 'Temperature at 2 m above the displacement height',
            'pattern': 'MERRA_MAT1NXSLV_T2M_*.tif',
            'url': 'http://goldsmr2.sci.gsfc.nasa.gov:80/opendap/MERRA/MAT1NXSLV.5.2.0',
            'source': 'MERRA%s.prod.assim.tavg1_2d_slv_Nx.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 11),
            'latency': -45
        },
        'T10M': {
            'description': 'Temperature at 10 m above the displacement height',
            'pattern': 'MERRA_MAT1NXSLV_T10M_*.tif',
            'url': 'http://goldsmr2.sci.gsfc.nasa.gov:80/opendap/MERRA/MAT1NXSLV.5.2.0',
            'source': 'MERRA%s.prod.assim.tavg1_2d_slv_Nx.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 11),
            'latency': -45
        },
    }

    # Not used anywhere
    _launchdate = {
    }

    # Should this be specified on a per asset basis?
    _defaultresolution = (0.666666666666667, -0.50)

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(MerraAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        bname = os.path.splitext(bname)[0]

        print
        print "init"

        print "bname", bname
        parts = bname.split('_')
        print "parts", parts

        # MERRA_MAT1NXSLV_%s_%s_%s.tif

        self.sensor = parts[1]
        self.asset = parts[2]
        self.tile = parts[3]
        year = parts[4]
        doy = parts[5]

        print "self.sensor, self.asset, self.tile, year, doy:", self.sensor, self.asset, self.tile, year, doy
        print

        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()

        datafiles = self.datafiles()


    def datafiles(self):
        """ Get list of datafiles from asset (if archive file) """

        print "self.filename", self.filename

        # try:
        #     datafiles = File2List(indexfile)
        # except:
        #     datafiles = []


    # def datafiles(self):

    #     indexfile = self.filename + '.index'

    #     # TODO: get File2List to handle missing file and empty file in the same way
    #     try:
    #         datafiles = File2List(indexfile)
    #     except:
    #         datafiles = []

    #     if not datafiles:
    #         gdalfile = gdal.Open(self.filename)
    #         subdatasets = gdalfile.GetSubDatasets()
    #         datafiles = [s[0] for s in subdatasets]
    #         List2File(datafiles, indexfile)

    #     return datafiles

    # @classmethod
    # def _filepattern(cls, asset, tile, date):
    #     """ used by fetch only """
    #     year, month, day = date.timetuple()[:3]
    #     # MERRA300.prod.assim.tavg1_2d_slv_Nx.20100101.hdf
    #     # pattern = ''.join(['(', asset, '.A', str(year), str(doy).zfill(3), '.', tile, '.005.\d{13}.hdf)'])
    #     datestr = str(year) + str(month).zfill(2) + str(day).zfill(2)
    #     # pattern = ''.join(['(', 'MERRA300.prod.assim.tavg1_2d_slv_Nx.', datestr, '.hdf)'])
    #     pattern = "(MERRA300.prod.assim.tavg1_2d_slv_Nx.%s.hdf)" % datestr
    #     print "pattern", pattern
    #     return pattern

    # @classmethod
    # def _remote_subdirs(cls, asset, tile, date):
    #     """ used by fetch only """
    #     year, month, day = date.timetuple()[:3]
    #     httploc = cls._assets[asset]['url']
    #     # mainurl = ''.join([httploc, '/', str(year), '.', '%02d' % month, '.', '%02d' % day])
    #     mainurl = "%s/%04d/%02d" % (httploc, year, month)
    #     # mainurl = ''.join([httploc, '/', str(year), '/', '%02d' % month])
    #     print "mainurl", mainurl
    #     return mainurl


    @classmethod
    def fetch(cls, asset, tile, date):

        from pydap.client import open_url
        from agspy.utils import raster

        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)

        outdir = cls.Repository.spath()

        tilesvector = cls.Repository.tiles_vector()

        print "asset", asset
        print "tile", tile
        print "date", date
        print "outdir", outdir

        for fid in tilesvector.get_fids():
            feature = tilesvector.get_feature(fid)
            if feature['tileid'] == tile:
                break

        bounds = eval(feature['bounds'])
        print "bounds", bounds
        x0 = bounds[0]
        y0 = bounds[1]

        y1 = bounds[3]

        RES = (0.666666666666667, 0.50)
        ORIG = (-180., -90.)

        dx = 12.0
        dy = 10.0

        nx = int(round(dx/RES[0]))
        ny = int(round(dy/RES[1]))

        # verify bounds
        h = int(tile[1:3])
        v = int(tile[4:6])

        print h, v
        x0_alt = ORIG[0] + dx*(h - 1)
        y0_alt = ORIG[1] + dy*(18 - v)

        print "x0, y0", x0, y0
        print "x0_alt, y0_alt", x0_alt, y0_alt
        assert x0 == x0_alt
        assert y0 == y0_alt


        ix0 = int(round((x0 - ORIG[0])/RES[0]))
        iy0 = int(round((y0 - ORIG[1])/RES[1]))
        ix1 = ix0 + nx
        iy1 = iy0 + ny

        print "x0, y0", x0, y0
        print "nx, ny", nx, ny
        print "ix0, iy0", ix0, iy0
        print "ix1, iy1", ix1, iy1

        url = cls._assets[asset]['url']

        if date.year > 2000:
            ver = '300'
        elif date.year > 1992:
            ver = '200'
        else:
            ver = '100'

        source = cls._assets[asset]['source'] % (ver, date.year, date.month, date.day)

        loc = "%s/%04d/%02d/%s" % (url, date.year, date.month, source)

        print "loc", loc

        dataset = open_url(loc)

        print dataset

        keys = dataset.keys()
        names = []
        for i, key in enumerate(keys):
            x = dataset[key]
            # print i, key, x.shape
            if x.shape == (24, 361, 540):
                names.append(key)        

        # print names
        # print len(names)

        assert asset in names, "asset name is not in SDS list"

        data = dataset[asset][:, iy0:iy1, ix0:ix1].astype('float32')
        data = data[:,::-1,:]

        print data.shape
        print type(data)
        print data.dtype

        proj = raster.create_proj(4326)
        geo = (x0, RES[0], 0.0, y1, 0.0, -RES[1])

        description = cls._assets[asset]['description']

        meta = {'ASSET':asset, 'TILE':tile, 'DATE':str(date.date()),
            'DESCRIPTION':description, 'SOURCE':loc}

        names = ['%02d30GMT' % i for i in range(24)]

        doy = date.strftime('%j')

        print "asset, tile, date.year, doy", asset, tile, date.year, doy
        outfilename = "MERRA_MAT1NXSLV_%s_%s_%4d_%s.tif" % (asset, tile, date.year, doy)
        outpath = os.path.join(outdir, outfilename)

        print "Writing to", outpath
        raster.write_raster(outpath, data, proj, geo, meta, bandnames=names)



class MerraData(Data):
    """ A tile of data (all assets and products) """
    name = 'Merra'
    Asset = MerraAsset
    _pattern = '*.tif' 
    _products = {
        'temp_modis': {
            'description': 'Air temperature data',
            'assets': ['TS', 'T2M', 'T10M']
        },
    }
    

    def process(self, products, **kwargs):

        print kwargs
        print products

        for key, val in products.items():

            outfname = os.path.join(self.path, self.basename + '_' + key)        
            VerboseOut("outfname: %s" % outfname, 4)

            if val[0] == "temp_modis":
                VERSION = "1.0"
                assets = self._products['temp_modis']['assets']

                allsds = []
                missingassets = []
                availassets = []
                assetids = []

                for asset in assets:
                    try:
                        sds = self.assets[asset].datafiles()
                    except Exception,e:
                        missingassets.append(asset)
                    else:
                        assetids.append(assets.index(asset))
                        availassets.append(asset)
                        allsds.extend(sds)

                if missingassets:
                    VerboseOut('There are missing assets: %s,%s,%s' % (str(self.date), str(self.id), str(missingassets)), 4)
                    continue

                print "assetids", assetids
                print "availassets", availassets
                print "missingassets", missingassets
                for i,sds in enumerate(allsds):
                    print "i, sds", i,sds




def main():
    DataInventory.main(MerraData)
