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
import re
import time
import datetime
import urllib
from osgeo import gdal
from collections import OrderedDict

from gippy.data.core import Asset, Tile, Data
from gippy.utils import File2List, List2File, VerboseOut

from pdb import set_trace


class ModisAsset(Asset):
    _rootpath = '/titan/data/modis'

    _stagedir = os.path.join(_rootpath, 'stage')

    _sensors = {
        'MOD': {'description': 'Terra'},
        'MYD': {'description': 'Aqua'},
        'MCD': {'description': 'Combined'}
    }
    _assets = {
        'MOD11A1': {
            'pattern': 'MOD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005'
        },
        'MYD11A1': {
            'pattern': 'MYD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.005'
        }
    }

    _launchdate = {
        'MOD': datetime.date(1999, 12, 1),
        'MYD': datetime.date(2002, 5, 1),
        'MCD': None
    }

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(ModisAsset, self).__init__(filename)

        print "asset. filename"
        print filename

        self.asset = self.basename[0:7]
        self.tile = self.basename[17:23]
        year = self.basename[9:13]
        doy = self.basename[13:16]
        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()
        self.sensor = self.basename[:3]
        self.basename = 'MODIS'+self.tile+'_'+year+doy

        datafiles = self.datafiles()

        # I don't understand what this is used for
        # because the asset could be many different things
        self.products = {'sds1': datafiles[0]}


    def datafiles(self):
        indexfile = self.filename + '.index'

        print indexfile

        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:
            gdalfile = gdal.Open(self.filename)
            subdatasets = gdalfile.GetSubDatasets()
            datafiles = [s[0] for s in subdatasets]
            List2File(datafiles, indexfile)
        return datafiles

    @classmethod
    def fetch(cls, asset, tile, date):
        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)

        httploc = cls._assets[asset]['url']

        year, month, day = date.timetuple()[:3]
        doy = date.timetuple()[7]

        pattern = ''.join(['(', asset, '.A', str(year), str(doy).zfill(3), '.', tile, '.005.\d{13}.hdf)'])
        mainurl = ''.join([httploc, '/', str(year), '.', '%02d' % month, '.', '%02d' % day])

        print "mainurl"
        print mainurl

        print "pattern"
        print pattern

        try:
            VerboseOut('opening listing', 1)
            listing = urllib.urlopen(mainurl).readlines()
        except Exception, e:
            listing = None
            print 'unable to access %s' % mainurl
            return 2

        cpattern = re.compile(pattern)
        name = None

        success = False

        for item in listing:
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                print 'found in', item.strip()
                name = cpattern.findall(item)[0]
                print 'it is', name
                url = ''.join([mainurl, '/', name])
                print 'the url is', url
                try:
                    urllib.urlretrieve(url, os.path.join(cls._stagedir, name))
                    print "retrieved %s" % name
                    success = True
                except Exception, e:
                    print e
                    print 'unable to retrieve %s from %s' % (name, url)

        if not success:
            print "did not find a match for %s in listing of %s" % (pattern, mainurl)
            return 1
        else:
            return 0


class ModisTile(Tile):
    """ A tile of data (all assets and products) """
    Asset = ModisAsset
    _prodpattern = '*.tif'
    _products = OrderedDict([
        ('temp', {
            'description': 'Surface temperature observations',
            # the list of asset types associated with this product
            'assets': ['MOD11A1', 'MYD11A1'],
        }),
        ('sds1', {
            'description': 'First SDS in the file',
            # the list of asset types associated with this product
            'assets': ['MOD11A1'],
        }),
    ])


class ModisData(Data):
    """ Represent a single day and temporal extent of MODIS data along with product variations """

    name = 'Modis'

    _tiles_vector = 'modis_sinusoidal_grid_world.shp'

    _defaultresolution = [926.625433138333392, -926.625433139166944]

    Tile = ModisTile

    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")
        h = str(int(feature.GetField(fldindex_h))).zfill(2)
        v = str(int(feature.GetField(fldindex_v))).zfill(2)
        tile = "h%sv%s" % (h, v)
        return tile


def main():
    ModisData.main()
