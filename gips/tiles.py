#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  mhanson@ags.io
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
from glob import glob
from itertools import groupby

import gippy
from gips.core import SpatialExtent, TemporalExtent, Products, DataNotFoundException
from gips.utils import VerboseOut, Colors, basename, mosaic


class Files(object):
    """ Collection of geospatial imagery files for single date and spatial region """

    def __init__(self, filenames=None):
        """ filenames naming convention: 'date_sensor_product' or 'tileid_date_sensor_product' """
        self.filenames = {}
        self.sensors = {}
        self.date = None
        self.AddFiles(filenames)

    def __getitem__(self, key):
        """ Get filename for product key """
        if type(key) == tuple:
            return self.filenames[key]
        else:
            return self.filenames[(self.sensor_set[0], key)]

    def __str__(self):
        """ Text representation """
        return '%s: %s' % (self.date, ' '.join(self.product_set))

    @property
    def doy(self):
        return self.date.strftime('%j')

    @property
    def numfiles(self):
        return len(self.filenames)

    @property
    def path(self):
        """ Get path to where files are stored (assuming same path for all) """
        # TODO - is same path for all a safe assumption?
        return os.path.abspath(self.filenames.values()[0])

    @property
    def sensor_set(self):
        """ Return list of sensors used """
        return list(set(sorted(self.sensors.values())))

    @property
    def products(self):
        """ Get list of products """
        return sorted([k[1] for k in self.filenames.keys()])

    @property
    def product_set(self):
        """ Return list of products available """
        return list(set(self.products))

    def AddFiles(self, filenames):
        """ Add filenames to existing filenames """
        for f in filenames:
            bname = basename(f)
            parts = bname.split('_')
            if len(parts) < 3 or len(parts) > 4:
                raise Exception('%s does not follow naming convention' % f)
            offset = 1 if len(parts) == 4 else 0
            if self.date is None:
                self.date = datetime.strptime(parts[0 + offset], '%Y%j').date()
            else:
                date = datetime.strptime(parts[0 + offset], '%Y%j').date()
                if date != self.date:
                    raise Exception('Mismatched dates: %s' % ' '.join(filenames))
            sensor = parts[1 + offset]
            product = parts[2 + offset]
            self.filenames[(sensor, product)] = f
            # TODO - currently assumes single sensor for each product
            self.sensors[product] = sensor

    def update(self, files):
        """ Merge files dictionary into this instance """
        if files.numfiles > 0:
            self.AddFiles(files.filenames)

    def open(self, product, sensor=None, update=False):
        """ Open and return a GeoImage """
        if sensor is None:
            sensor = self.sensors[product]
        fname = self.filenames[(sensor, product)]
        if os.path.exists(fname):
            return gippy.GeoImage(fname)
        else:
            raise Exception('%s_%s product does not exist' % (sensor, product))

    # TODO - make general product_filter function
    def masks(self, patterns=None):
        """ List all products that are masks """
        if not patterns:
            patterns = ['acca', 'fmask', 'mask']
        m = []
        for p in self.products:
            if any(pattern in p for pattern in patterns):
                m.append(p)
        return m

    def pprint_header(self):
        """ Print product inventory header """
        return Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE') + '{:^10}'.format('Products') + Colors.OFF

    def pprint(self, dformat='%j', colors=None):
        """ Print product inventory for this date """
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        if colors is None:
            sys.stdout.write('  '.join(sorted(self.products)))
        else:
            for p in sorted(self.products):
                sys.stdout.write(colors[self.sensors[p]] + p + Colors.OFF + '  ')
        sys.stdout.write('\n')

    """
    @classmethod
    def discover(cls, _files):
        # Factory function returns instance for every date in 'files' (filenames or directory)
        if not isinstance(_files, list) and os.path.isdir(os.path.abspath(_files)):
            _files = glob(os.path.join(_files, '*.tif'))

        # Check files for 3 or 4 parts and a valid date
        files = []
        for f in _files:
            parts = basename(f).split('_')
            if len(parts) == 3 or len(parts) == 4:
                try:
                    datetime.strptime(parts[len(parts) - 3], '%Y%j')
                except:
                    continue
                files.append(f)

        if len(files) == 0:
            raise DataNotFoundException('Files: No valid files given')

        # files will have 3 or 4 parts, so sind is 0 or 1
        sind = len(basename(files[0]).split('_')) - 3
        instances = []
        func = lambda x: datetime.strptime(basename(x).split('_')[sind], '%Y%j').date()
        for date, fnames in groupby(sorted(files), func):
            instances.append(cls(list(fnames)))
        return instances
    """


class Tiles(object):
    """ Collection of files for single date and multiple regions (tiles) """

    def __init__(self, dataclass, spatial=None, date=None, products=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.coverage      dict of tile id: %coverage with site
        self.tiles              dict of tile id: tile instance
        """
        self.dataclass = dataclass
        self.spatial = spatial if spatial is not None else SpatialExtent()
        self.date = date
        self.products = products if products is not None else Products(dataclass)

        # For each tile locate files/products
        VerboseOut('%s: searching %s tiles for products and assets' % (self.date, len(self.spatial)), 4)
        self.tiles = {}
        for t in self.spatial.tiles:
            try:
                tile = dataclass(t, self.date)
                good = tile.filter(**kwargs)
                # TODO - fix custom filter based on dataclass...sensor should pass to 'dataclass'
                if good:  # and tile.sensor in sensors:
                    self.tiles[t] = tile
            except DataNotFoundException:   # don't have data for this tile/date
                pass
            except Exception, e:            # re-raise all other exceptions
                raise Exception(e)

        if len(self.tiles) == 0:
            raise DataNotFoundException('No tiles found for %s ' % date)

    @classmethod
    def discover(cls, dataclass, spatial=None, temporal=None, products=None, **kwargs):
        """ Factory function, create tile class for each date for this dataclass """

        temporal = temporal if temporal is not None else TemporalExtent()
        instances = {}
        for date in temporal.prune_dates(spatial.available_dates):
            try:
                instances[date] = cls(dataclass, spatial, date, products, **kwargs)
            except:
                # No tiles for this date!
                VerboseOut('No tiles for %s' % date)
                raise Exception('No tiles for %s' % date)
                #VerboseOut(traceback.format_exc(), 5)
                continue

        return instances

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

    def process(self, overwrite=False, **kwargs):
        """ Determines what products need to be processed for each tile and calls Data.process """
        for tileid, tile in self.tiles.items():
            toprocess = {}
            for pname, args in self.requested_products.items():
                if pname not in tile.products or overwrite:
                    toprocess[pname] = args
            if len(toprocess) != 0:
                VerboseOut('Processing products for tile %s: %s' % (tileid, ' '.join(toprocess.keys())), 2)
                self.tiles[tileid].process(toprocess, **kwargs)

    def project(self, datadir, res, crop=False, nowarp=False, nomosaic=False, **kwargs):
        """ Create image of final product (reprojected/mosaiced) """
        start = datetime.now()
        bname = self.date.strftime('%Y%j')

        if self.site is None:
            nomosaic = True
            nowarp = True
        if not hasattr(res, "__len__"):
            res = [res, res]

        if nomosaic:
            # Tile project
            # TODO - allow hard and soft link options
            for t in self.tiles:
                tiledir = datadir.replace('TILEID', t)
                if not os.path.exists(tiledir):
                    os.makedirs(tiledir)
                for p in self.requested_products:
                    sensor = self.which_sensor(p)
                    filename = self.tiles[t].products[p]
                    fout = os.path.join(tiledir, t + '_' + bname + ('_%s_%s.tif' % (sensor, p)))
                    if not os.path.exists(fout):
                        try:
                            VerboseOut("Creating %s" % os.path.basename(fout))
                            if nowarp:
                                gippy.GeoImage(filename).Process(fout)
                            else:
                                # Warp each tile
                                gippy.CookieCutter([filename], fout, self.site, res[0], res[1], crop)
                        except Exception:
                            VerboseOut("Problem creating %s" % fout, 2)
                            VerboseOut(traceback.format_exc(), 3)
        else:
            # Shapefile project
            if not os.path.exists(datadir):
                os.makedirs(datadir)
            for product in self.requested_products:
                sensor = self.which_sensor(product)
                fout = os.path.join(datadir, bname + ('_%s_%s.tif' % (sensor, product)))
                if not os.path.exists(fout):
                    try:
                        filenames = [self.tiles[t].products[product] for t in self.tiles]
                        # TODO - cookiecutter should validate pixels in image.  Throw exception if not
                        if nowarp:
                            mosaic(filenames, fout, self.site)
                        else:
                            gippy.CookieCutter(filenames, fout, self.site, res[0], res[1], crop)
                    except:
                        VerboseOut("Problem projecting %s" % fout, 2)
                        VerboseOut(traceback.format_exc(), 3)
        t = datetime.now() - start
        VerboseOut('%s: created project files for %s tiles in %s' % (self.date, len(self.tiles), t), 2)

    def pprint_header(self):
        """ Print header info for coverage """
        header = Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE')
        for a in sorted(self.dataclass.Asset._assets.keys()):
            header = header + ('{:^10}'.format(a if a != '' else 'Coverage'))
        return header + '{:^10}'.format('Products') + Colors.OFF

    def asset_coverage(self):
        """ Calculates % coverage of site for each asset """
        asset_coverage = {}
        for a in self.dataclass.Asset._assets:
            cov = 0.0
            norm = float(len(self.coverage)) if self.spatial.site is None else 1.0
            for t in self.tiles:
                if a in self.tiles[t].assets:
                    cov = cov + (self.spatial.coverage[t][0] / norm)
            asset_coverage[a] = cov * 100
        return asset_coverage

    def product_coverage(self):
        """ Calculated % coverage of site for each product """
        pass
        #asset_coverage = self.asset_coverage()
        # loop through products and calculate overlap % for assets

    def pprint(self, dformat='%j', colors=None):
        """ Print coverage for each and every asset """
        #assets = [a for a in self.dataclass.Asset._assets]
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        asset_coverage = self.asset_coverage()
        for a in sorted(asset_coverage):
            color = ['', '']
            if colors is not None:
                s = self.which_sensor(a)
                if s is not None:
                    color = [colors[s], Colors.OFF]
            cov = asset_coverage[a]
            if cov > 0:
                sys.stdout.write(color[0] + '  {:>4.1f}%   '.format(cov) + color[1])
            else:
                sys.stdout.write('          ')
        products = [p for t in self.tiles for p in self.tiles[t].products]
        prods = []
        for p in set(products):
            if products.count(p) == len(self.tiles):
                prods.append(p)
        #prods = []
        #for t in self.tiles:
        #    for p in self.tiles[t].products:
                #prods.append(p)
        for p in sorted(set(prods)):
            color = colors[self.which_sensor(p)]
            sys.stdout.write('  ' + color + p + Colors.OFF)
        sys.stdout.write('\n')
