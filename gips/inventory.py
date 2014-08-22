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
import glob
import argparse
from datetime import datetime as dt
import traceback
import textwrap
from pprint import pprint
from itertools import groupby

import gippy
from gips.tiles import Tiles
from gips.utils import VerboseOut, parse_dates, Colors, basename
from gips.GeoVector import GeoVector
from gips.version import __version__
from pdb import set_trace


class Products(object):
    """ Collection of products (files) for a single date """

    def __init__(self, filenames):
        """ Find all data for this date and sensor """
        if not isinstance(filenames, list):
            raise TypeError('requires list (of filenames)')
        self.filenames = {}
        self.sensors = set()
        self.products = set()
        ind = len(os.path.basename(filenames[0]).split('_')) - 3
        for f in filenames:
            parts = basename(f).split('_')
            date = dt.strptime(parts[ind], '%Y%j').date()
            sensor = parts[1+ind]
            product = parts[2+ind]
            self.products.add(product)
            self.sensors.add(sensor)
            self.filenames[(product, sensor)] = f
        self.date = date

    def __getitem__(self, key):
        """ Return filename for product or (product, sensor) """
        if isinstance(key, tuple):
            return self.filenames[key]
        elif len(self.sensors) == 1:
            return self.filenames[(key, list(self.sensors)[0])]
        else:
            raise Exception('Need sensor key with multiple sensors')

    def __str__(self):
        return '%s: %s' % (self.date, ' '.join(self.products))

    @property
    def doy(self):
        """ Day of year """
        return self.date.strftime('%j')

    def open(self, product='', update=False):
        """ Open and return GeoImage """
        # TODO - assumes single sensor
        fname = self[product]
        if os.path.exists(fname):
            return gippy.GeoImage(fname, update)
        else:
            raise Exception('%s product does not exist' % product)

    @classmethod
    def discover(cls, files):
        """ Factory function returns instance for every date/sensor in filenames or directory """
        if not isinstance(files, list) and os.path.isdir(os.path.abspath(files)):
            files = glob.glob(os.path.join(files, '*.tif'))

        # files will have 3 or 4 parts, so ind is 0 or 1
        ind = len(basename(files[0]).split('_')) - 3

        instances = []
        for date, fnames in groupby(files, lambda x: dt.strptime(basename(x).split('_')[ind], '%Y%j').date()):
            instances.append(cls([f for f in fnames]))
        return instances


class Inventory(object):
    """ Base class for inventories """
    # TODO - add code for managing site(?) and common printing

    def __init__(self):
        pass

    def __getitem__(self, date):
        """ Indexing operator for class """
        # TODO - assumes single sensor for a date
        return self.data[date]

    @property
    def dates(self):
        """ Get sorted list of dates """
        return sorted(self.data.keys())

    @property
    def datestr(self):
        return '%s dates (%s - %s)' % (len(self.dates), self.dates[0], self.dates[-1])

    def _temporal_extent(self, dates, days):
        """ Temporal extent (define self.dates and self.days) """
        if dates is None:
            dates = '1984,2050'
        self.start_date, self.end_date = parse_dates(dates)
        if days:
            days = days.split(',')
        else:
            days = (1, 366)
        self.start_day, self.end_day = (int(days[0]), int(days[1]))


class ProjectInventory(Inventory):
    """ Inventory of project directory (collection of Products class) """

    def __init__(self, projdir='', products=[]):
        """ Create inventory of a GIPS project directory """
        self.projdir = os.path.abspath(projdir)
        if not os.path.exists(self.projdir):
            raise Exception('Directory %s does not exist!' % self.projdir)

        self.products = {}
        files = glob.glob(os.path.join(self.projdir, '*.tif'))
        self.product_set = set()
        for dat in Products.discover(files):
            self.products[dat.date] = dat
            self.product_set = self.product_set.union(dat.products)

        if not products:
            products = self.product_set
        self.requested_products = products

    @property
    def data(self):
        """ alias used by base class to self.products """
        return self.products

    def get_timeseries(self, product='', dates=None):
        """ Read all files as time series """
        if dates is None:
            dates = self.dates
        # TODO - multiple sensors
        filenames = [self.data[date][product] for date in dates]
        img = gippy.GeoImage(filenames)
        return img


class DataInventory(Inventory):
    """ Manager class for data inventories """
    _colors = [Colors.PURPLE, Colors.RED, Colors.GREEN, Colors.BLUE]

    def __init__(self, dataclass, site=None, tiles=None, dates=None, days=None, products=[], fetch=False, **kwargs):
        self.dataclass = dataclass
        Repository = dataclass.Asset.Repository

        self.site = site
        # default to all tiles
        if tiles is None and self.site is None:
            tiles = Repository.find_tiles()
        # if tiles provided, make coverage all 100%
        if tiles is not None:
            self.tiles = {}
            for t in tiles:
                self.tiles[t] = (1, 1)
        elif tiles is None and self.site is not None:
            self.tiles = Repository.vector2tiles(GeoVector(self.site), **kwargs)
        self._temporal_extent(dates, days)
        self.data = {}

        # Ensure product list if dictionary of products + arguments
        prod_dict = {}
        if isinstance(products, list):
            prod_dict = dict([p, [p]] for p in products)
        else:
            prod_dict = products

        # if no products specified and only 1 product available, use it
        if len(prod_dict) == 0 and len(dataclass._products) == 1:
            p = dataclass._products.keys()[0]
            prod_dict = {p: [p]}

        # seperate out standard (each tile processed) and composite products (using inventory)
        self.standard_products = {}
        self.composite_products = {}
        for p, val in prod_dict.items():
            if self.dataclass._products[val[0]].get('composite', False):
                self.composite_products[p] = val
            else:
                self.standard_products[p] = val
        if fetch:
            products = [val[0] for val in prod_dict.values()]
            try:
                dataclass.fetch(products, self.tiles, (self.start_date, self.end_date), (self.start_day, self.end_day))
            except:
                VerboseOut(traceback.format_exc(), 3)
            dataclass.Asset.archive(Repository.spath())

        # get all potential matching dates for tiles
        dates = []
        for t in self.tiles:
            try:
                for date in Repository.find_dates(t):
                    day = int(date.strftime('%j'))
                    if (self.start_date <= date <= self.end_date) and (self.start_day <= day <= self.end_day):
                        if date not in dates:
                            dates.append(date)
            except:
                VerboseOut(traceback.format_exc(), 3)

        self.numfiles = 0
        for date in sorted(dates):
            try:
                dat = Tiles(dataclass=dataclass, site=self.site, tiles=self.tiles,
                            date=date, products=self.standard_products, **kwargs)
                self.data[date] = dat
                self.numfiles = self.numfiles + len(dat.tiles)
            except Exception, e:
                pass
                #VerboseOut(traceback.format_exc(), 3)
        sensors = set([self.data[date].sensor for date in self.data])
        self.sensors = {}
        for i, s in enumerate(sensors):
            self.sensors[s] = self.dataclass.Asset._sensors[s]
            self.sensors[s]['color'] = Colors.BOLD + self._colors[i]
        if len(dates) == 0:
            raise Exception("No matching files in inventory!")

    def process(self, **kwargs):
        """ Process data in inventory """
        if len(self.standard_products) + len(self.composite_products) == 0:
            raise Exception('No products specified!')
        sz = self.numfiles
        if len(self.standard_products) > 0:
            start = dt.now()
            VerboseOut('Processing %s files: %s' % (sz, ' '.join(self.standard_products)), 1)
            for date in self.dates:
                try:
                    self.data[date].process(**kwargs)
                except:
                    pass
            VerboseOut('Completed processing in %s' % (dt.now()-start), 1)
        if len(self.composite_products) > 0:
            start = dt.now()
            VerboseOut('Processing %s files into composites: %s' % (sz, ' '.join(self.composite_products)), 1)
            self.dataclass.process_composites(self, self.composite_products, **kwargs)
            VerboseOut('Completed processing in %s' % (dt.now()-start), 1)

    def project(self, datadir=None, res=None, **kwargs):
        """ Create project files for data in inventory """
        self.process(**kwargs)
        start = dt.now()

        # formulate project directory name
        if res is None:
            res = self.dataclass.Asset._defaultresolution
        sitename = basename(self.site) if self.site else 'tiles'
        if datadir is None:
            datadir = '%s_%s_%sx%s' % (self.dataclass.name, sitename, res[0], res[1])
        if not os.path.exists(datadir):
            os.makedirs(datadir)

        VerboseOut('Creating GIPS project %s' % datadir)
        VerboseOut('  Dates: %s' % self.datestr)
        VerboseOut('  Products: %s' % ' '.join(self.standard_products))
        for date in self.dates:
            self.data[date].project(datadir=datadir, res=res, **kwargs)
        VerboseOut('Completed GIPS project in %s' % (dt.now()-start))
        return ProjectInventory(datadir)

    def print_inv(self, md=False, compact=False):
        """ Print inventory """
        self.print_tile_coverage()
        print
        if self.site is not None:
            print Colors.BOLD + 'Asset Coverage for site %s' % basename(self.site) + Colors.OFF
        else:
            print Colors.BOLD + 'Asset Holdings' + Colors.OFF
        # Header
        #datestr = ' Month/Day' if md else ' Day of Year'
        header = Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE')
        for a in sorted(self.dataclass.Asset._assets.keys()):
            header = header + ('{:^10}'.format(a if a != '' else 'Coverage'))
        header = header + '{:^10}'.format('Products')
        header = header + Colors.OFF
        print header
        dformat = '%m-%d' if md else '%j'
        # Asset inventory
        oldyear = 0
        formatstr = '\n{:<12}' if compact else '{:<12}\n'
        for date in self.dates:
            sensor = self.sensors[self.data[date].sensor]
            if date.year != oldyear:
                sys.stdout.write(Colors.BOLD + formatstr.format(date.year) + Colors.OFF)
            if compact:
                dstr = ('{:^%s}' % (7 if md else 4)).format(date.strftime(dformat))
                sys.stdout.write(sensor['color'] + dstr + Colors.OFF)
            else:
                self.data[date].print_assets(dformat, color=sensor['color'])
            oldyear = date.year
        if self.numfiles != 0:
            VerboseOut("\n\n%s files on %s dates" % (self.numfiles, len(self.dates)), 1)
            self.print_legend()

    def print_tile_coverage(self):
        if self.site is not None:
            print Colors.BOLD + '\nTile Coverage'
            print Colors.UNDER + '{:^8}{:>14}{:>14}'.format('Tile', '% Coverage', '% Tile Used') + Colors.OFF
            for t in sorted(self.tiles):
                print "{:>8}{:>11.1f}%{:>11.1f}%".format(t, self.tiles[t][0]*100, self.tiles[t][1]*100)

    def print_legend(self):
        print Colors.BOLD + '\nSENSORS' + Colors.OFF
        sensors = {}
        for key in sorted(self.sensors):
            scode = key+': ' if key != '' else ''
            print self.sensors[key]['color'] + '%s%s' % (scode, self.sensors[key]['description']) + Colors.OFF

    @staticmethod
    def main(cls):
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(description='GIPS %s Data Utility v%s' % (cls.name, __version__),
                                          formatter_class=argparse.RawTextHelpFormatter)
        subparser = parser0.add_subparsers(dest='command')

        # Archive
        parser = subparser.add_parser('archive', help='Move files from current directory to data archive')
        parser.add_argument('--keep', help='Keep files after adding to archive', default=False, action='store_true')
        parser.add_argument('--recursive', help='Iterate through subdirectories', default=False, action='store_true')
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

        invparser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = invparser.add_argument_group('inventory arguments')
        group.add_argument('-s', '--site', help='Vector file for region of interest', default=None)
        group.add_argument('-t', '--tiles', nargs='*', help='Tile designations', default=None)
        group.add_argument('-d', '--dates', help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
        group.add_argument('--days', help='Include data within these days of year (doy1,doy2)', default=None)
        group.add_argument('--sensors', help='Sensors to include', nargs='*', default=None)
        group.add_argument('--%cov', dest='pcov', help='Threshold of %% tile coverage over site', default=0, type=int)
        group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)
        group.add_argument('--fetch', help='Fetch any missing data (if supported)', default=False, action='store_true')
        group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
        extra = []
        for arg, kwargs in cls.extra_arguments().items():
            extra.append(kwargs['dest'])
            group.add_argument(arg, **kwargs)

        parents = [invparser, cls.arg_parser()]

        # Inventory
        parser = subparser.add_parser('inventory', help='Get Inventory', parents=parents, formatter_class=dhf)
        parser.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)
        parser.add_argument('--compact', help='Print only dates (no coverage)', default=False, action='store_true')

        # Processing
        parserp = subparser.add_parser('process', help='Process scenes', parents=parents, formatter_class=dhf)
        group = parserp.add_argument_group('Processing Options')
        group.add_argument('--suffix', help='Suffix on end of filename (before extension)', default='')
        group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')
        group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)

        # Project
        parser = subparser.add_parser('project', help='Create project', parents=parents, formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--suffix', help='Suffix on end of filename (before extension)', default='')
        group.add_argument('--crop', help='Crop output down to minimum bounding box (if warping)', default=False, action='store_true')
        group.add_argument('--nowarp', help='Mosaic, but do not warp', default=False, action='store_true')
        group.add_argument('--res', nargs=2, help='Resolution of (warped) output rasters', default=None, type=float)
        group.add_argument('--datadir', help='Directory to save project files (default auto-generated)', default=None)
        group.add_argument('--format', help='Format for output file', default="GTiff")
        group.add_argument('--chunksize', help='Chunk size in MB', type=float, default=512.0)

        args = parser0.parse_args()

        gippy.Options.SetVerbose(args.verbose)
        if 'format' in args:
            gippy.Options.SetDefaultFormat(args.format)

        VerboseOut('GIPS %s command line utility v%s' % (cls.name, __version__), 1)

        if args.command == 'archive':
            cls.Asset.archive(recursive=args.recursive, keep=args.keep)
            exit(1)

        try:
            suffix = '_' + args.suffix if args.suffix != '' else ''
        except:
            suffix = ''
        products = {}
        for p in cls._products:
            if p != '':
                val = eval('args.%s' % p)
                if val not in [None, False]:
                    if val is True:
                        products[p] = [p]
                    elif isinstance(val, list):
                        key = p
                        for i in val:
                            key = key + '_' + i
                        products[key+suffix] = [p] + val
                    else:
                        products[p+'_'+val+suffix] = [p, val]
        kwargs = dict(zip(extra, [eval('args.%s' % a) for a in extra]))

        try:
            inv = cls.inventory(
                site=args.site, dates=args.dates, days=args.days, tiles=args.tiles,
                products=products, pcov=args.pcov, ptile=args.ptile, fetch=args.fetch, sensors=args.sensors, **kwargs)
            if args.command == 'inventory':
                inv.print_inv(args.md, compact=args.compact)
            elif args.command == 'process':
                gippy.Options.SetChunkSize(args.chunksize)
                inv.process(overwrite=args.overwrite, **kwargs)
            elif args.command == 'project':
                gippy.Options.SetChunkSize(args.chunksize)
                projinv = inv.project(datadir=args.datadir, res=args.res, crop=args.crop, nowarp=args.nowarp)
            else:
                VerboseOut('Command %s not recognized' % cmd, 0)
        except Exception, e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
