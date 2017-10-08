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

import sys
import os
from datetime import datetime
import traceback

import gippy
from gippy.algorithms import CookieCutter
from gips.utils import VerboseOut, Colors, mosaic, mkdir
from gips import utils


class Tiles(object):
    """ Collection of files for single date and one or more regions (tiles) """

    def __init__(self, dataclass, spatial, date=None, products=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.coverage      dict of tile id: %coverage with site
        self.tiles         mapping of tile IDs to Data instances, eg for modis:
                           {'h12v04': <modis.modisData object at 0x7fb34ab9e550>}
        """
        self.dataclass = dataclass
        self.spatial = spatial
        self.products = products if products is not None else dataclass.RequestedProducts()
        self.date = date
        # For each tile locate files/products
        self.tiles = {}

    def __len__(self):
        return len(self.tiles)

    def __getitem__(self, key):
        return self.tiles[key]

    @property
    def sensor_set(self):
        """ Return list of sensors used in all tiles """
        s = set()
        for t in self.tiles:
            s.update(self.tiles[t].sensor_set)
        return list(s)

    def which_sensor(self, key):
        """ Get sensor code used for provided asset or product key """
        for t in self.tiles:
            if key in self.tiles[t].sensors:
                return self.tiles[t].sensors[key]

    def process(self, *args, **kwargs):
        """ Calls process for each tile """
        [t.process(*args, products=self.products.products, **kwargs) for t in self.tiles.values()]

    def mosaic(self, datadir, res=None, interpolation=0, crop=False,
               overwrite=False, alltouch=False):
        """For each product, combine its tiles into a single mosaic.

        Warp if res provided."""
        if self.spatial.site is None:
            raise Exception('Site required for creating mosaics')
        start = datetime.now()
        bname = self.date.strftime('%Y%j')
        for product in self.products.products:
            sensor = self.which_sensor(product)
            if sensor is None:
                continue
            # create data directory when it is needed
            mkdir(datadir)
            # TODO - this is assuming a tif file.  Use gippy FileExtension function when it is exposed
            fout = os.path.join(datadir, '%s_%s_%s' % (bname, sensor, product)) + '.tif'
            if not os.path.exists(fout) or overwrite:
                with utils.error_handler("Error mosaicking " + fout, continuable=True):
                    filenames = [self.tiles[t].filenames[(sensor, product)] for t in self.tiles]
                    images = gippy.GeoImages(filenames)
                    if self.spatial.site is not None and res is not None:
                        CookieCutter(
                            images, self.spatial.site, fout, res[0], res[1],
                            crop, interpolation, {}, alltouch,
                        )
                    else:
                        mosaic(images, fout, self.spatial.site)
        t = datetime.now() - start
        VerboseOut('%s: created project files for %s tiles in %s' % (self.date, len(self.tiles), t), 2)

    def asset_coverage(self):
        """ Calculates % coverage of site for each asset """
        asset_coverage = {}
        for a in self.dataclass.Asset._assets:
            cov = 0.0
            norm = float(len(self.spatial.coverage)) if self.spatial.site is None else 1.0
            for t in self.tiles:
                if a in self.tiles[t].assets:
                    cov = cov + (self.spatial.coverage[t][0] / norm)
            asset_coverage[a] = cov * 100
        return asset_coverage

    def product_coverage(self):
        """ Calculated % coverage of site for each product """
        coverage = {}
        for p in self.dataclass._products.keys():
            assets = self.dataclass.products2assets([p])
            norm = float(len(self.spatial.coverage)) if self.spatial.site is None else 1.0
            cov = 0.0
            for t in self.tiles:
                allassets = True
                for a in assets:
                    if a not in self.tiles[t].assets:
                        allassets = False
                if allassets:
                    cov = cov + (self.spatial.coverage[t][0] / norm)
            coverage[p] = cov * 100
        return coverage

    # TODO - remove this nonsense
    def pprint_asset_header(self):
        """ Print asset header """
        self.dataclass.pprint_asset_header()

    def _colorize_product(self, prod, colors=None):
        color = ['', '']
        if colors is not None:
            s = self.which_sensor(prod)
            if s is not None:
                color = [colors[s], Colors.OFF]
        return color[0] + prod + color[1]


    def pprint(self, dformat='%j', colors=None):
        """ Print coverage for each and every asset """
        #assets = [a for a in self.dataclass.Asset._assets]
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        asset_coverage = self.asset_coverage()
        for a in sorted(asset_coverage):
            cov = asset_coverage[a]
            if cov > 0:
                text = self._colorize_product(
                    '  {:>4.1f}%   '.format(cov), colors
                )
            else:
                text = '          '
            sys.stdout.write(text)

        products = [p for t in self.tiles for p in self.tiles[t].products]
        # Check product is available for all tiles before reporting as processed
        prods = []
        for p in set(products):
            if products.count(p) == len(self.tiles):
                prods.append(p)
        for p in sorted(set(prods)):
            text = self._colorize_product(p, colors)
            sys.stdout.write('  ' + text)
        sys.stdout.write('\n')
