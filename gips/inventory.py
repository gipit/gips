#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
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
from itertools import groupby
import numpy

import gippy
from gips.tiles import Tiles
from gips.utils import VerboseOut, parse_dates, Colors, basename
from gips.GeoVector import GeoVector
from gips.version import __version__
from pdb import set_trace


class Products(object):
    """ Collection of products (files) for a single date and sensor """

    def __init__(self, filenames):
        """ Create products with these filenames (should all be same date/sensor) """
        # TODO - check that all dates and sensors are the same for all filenames
        if not isinstance(filenames, list):
            raise TypeError('requires list (of filenames)')
        self.filenames = {}
        self.products = set()
        ind = len(os.path.basename(filenames[0]).split('_')) - 3
        for f in filenames:
            parts = basename(f).split('_')
            date = dt.strptime(parts[ind], '%Y%j').date()
            self.sensor = parts[1 + ind]
            product = parts[2 + ind]
            # is it a product or a mask
            self.products.add(product)
            self.filenames[product] = f
        self.date = date

    def __getitem__(self, key):
        """ Return filename for product or (product, sensor) """
        return self.filenames[key]

    def __str__(self):
        return '%s (%s): %s' % (self.date, self.sensor, ' '.join(sorted(self.products)))

    @property
    def numfiles(self):
        return len(self.filenames)

    @property
    def doy(self):
        """ Day of year """
        return self.date.strftime('%j')

    def masks(self, patterns=None):
        if not patterns:
            patterns = ['acca', 'fmask', 'mask']
        m = []
        for p in self.products:
            if any(pattern in p for pattern in patterns):
                m.append(p)
        return m

    def open(self, product='', update=False):
        """ Open and return GeoImage """
        # TODO - assumes single sensor
        fname = self[product]
        if os.path.exists(fname):
            return gippy.GeoImage(fname, update)
        else:
            raise Exception('%s product does not exist' % product)

    def pprint_header(self):
        return Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE') + '{:^10}'.format('Products') + Colors.OFF

    def pprint(self, dformat='%j', color=''):
        """ Print data products """
        sys.stdout.write(color)
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        sys.stdout.write('  '.join(sorted(self.products)))
        if color != '':
            sys.stdout.write(Colors.OFF)
        sys.stdout.write('\n')

    @classmethod
    def discover(cls, files):
        """ Factory function returns instance for every date/sensor in filenames or directory """
        if not isinstance(files, list) and os.path.isdir(os.path.abspath(files)):
            files = glob.glob(os.path.join(files, '*.tif'))
        # files will have 3 or 4 parts, so ind is 0 or 1
        ind = len(basename(files[0]).split('_')) - 3

        instances = []
        for date, fnames in groupby(sorted(files), lambda x: dt.strptime(basename(x).split('_')[ind], '%Y%j').date()):
            for sensor, fnames2 in groupby(sorted(fnames), lambda x: basename(x).split('_')[1]):
                instances.append(cls(list(fnames2)))
        return instances


class Inventory(object):
    """ Base class for inventories """
    _colors = [Colors.PURPLE, Colors.RED, Colors.GREEN, Colors.BLUE]

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
    def sensors(self):
        sensors = {}
        for i, s in enumerate(sorted(set([self.data[d].sensor for d in self.data.keys()]))):
            sensors[s] = {'color': self._colors[i]}
        return sensors

    @property
    def numfiles(self):
        """ Total number of files in inventory """
        return sum([d.numfiles for d in self.data.values()])

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

    def pprint(self, md=False, compact=False):
        print self.data[self.data.keys()[0]].pprint_header()
        dformat = '%m-%d' if md else '%j'
        oldyear = 0
        formatstr = '\n{:<12}' if compact else '{:<12}\n'
        for date in self.dates:
            scode = self.data[date].sensor
            if date.year != oldyear:
                sys.stdout.write(Colors.BOLD + formatstr.format(date.year) + Colors.OFF)
            if compact:
                dstr = ('{:^%s}' % (7 if md else 4)).format(date.strftime(dformat))
                sys.stdout.write(self.sensors[scode]['color'] + dstr + Colors.OFF)
            else:
                self.data[date].pprint(dformat, color=self.sensors[scode]['color'])
            oldyear = date.year
        if self.numfiles != 0:
            VerboseOut("\n\n%s files on %s dates" % (self.numfiles, len(self.dates)), 1)
            self.print_legend()

    def print_legend(self):
        print Colors.BOLD + '\nSENSORS' + Colors.OFF
        for key in sorted(self.sensors):
            try:
                desc = self.dataclass.Asset._sensors[key]['description']
                scode = key + ': ' if key != '' else ''
            except:
                desc = ''
                scode = key
            print self.sensors[key]['color'] + '%s%s' % (scode, desc) + Colors.OFF


class ProjectInventory(Inventory):
    """ Inventory of project directory (collection of Products class) """

    def __init__(self, projdir='', products=[]):
        """ Create inventory of a GIPS project directory """
        self.projdir = os.path.abspath(projdir)
        if not os.path.exists(self.projdir):
            raise Exception('Directory %s does not exist!' % self.projdir)

        self.products = {}
        files = glob.glob(os.path.join(self.projdir, '*.tif'))
        product_set = set()
        try:
            for dat in Products.discover(files):
                self.products[dat.date] = dat
                product_set = product_set.union(dat.products)
            if not products:
                products = product_set
            self.requested_products = products
        except:
            raise Exception("%s does not appear to be a GIPS project directory" % self.projdir)

    @property
    def data(self):
        """ alias used by base class to self.products """
        return self.products

    def product_list(self, date):
        """ Intersection of available products for this date and requested products """
        return self.products[date].products.intersection(self.requested_products)

    def new_image(self, filename, dtype=gippy.GDT_Byte, numbands=1, nodata=None):
        """ Create new image with the same template as the files in project """
        img = gippy.GeoImage(self.data[self.dates[0]].open(self.requested_products[0]))
        imgout = gippy.GeoImage(filename, img, dtype, numbands)
        img = None
        if nodata is not None:
            imgout.SetNoData(nodata)
        return imgout

    def get_data(self, dates):
        """ Read all files as time series, stacking all products """
        # TODO - change to absolute dates
        days = numpy.array([int(d.strftime('%j')) for d in dates])

        imgarr = []
        for p in self.requested_products:
            gimg = self.get_timeseries(p, dates=dates)
            # TODO - move numpy.squeeze into swig interface file?
            arr = numpy.squeeze(gimg.TimeSeries(days.astype('float64')))
            arr[arr == gimg[0].NoDataValue()] = numpy.nan

            if len(days) == 1:
                dims = arr.shape
                arr = arr.reshape(1, dims[0], dims[1])
            imgarr.append(arr)
        data = numpy.vstack(tuple(imgarr))
        return data

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

        if products is None:
            products = dataclass._products
        prod_dict = dict([p, p.split('-')] for p in products)

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

        for date in sorted(dates):
            try:
                dat = Tiles(dataclass=dataclass, site=self.site, tiles=self.tiles,
                            date=date, products=self.standard_products, **kwargs)
                self.data[date] = dat
            except Exception:
                pass
                #VerboseOut(traceback.format_exc(), 3)
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
                    VerboseOut(traceback.format_exc(), 3)
                    pass
            VerboseOut('Completed processing in %s' % (dt.now() - start), 1)
        if len(self.composite_products) > 0:
            start = dt.now()
            VerboseOut('Processing %s files into composites: %s' % (sz, ' '.join(self.composite_products)), 1)
            self.dataclass.process_composites(self, self.composite_products, **kwargs)
            VerboseOut('Completed processing in %s' % (dt.now() - start), 1)

    def project(self, datadir=None, suffix='', res=None, **kwargs):
        """ Create project files for data in inventory """
        self.process(**kwargs)
        start = dt.now()

        if suffix != '':
            suffix = '_' + suffix

        # formulate project directory name
        if res is None:
            res = self.dataclass.Asset._defaultresolution
        if res[0] == res[1]:
            resstr = str(res[0])
        else:
            resstr = '%sx%s' % (res[0], res[1])
        sitename = basename(self.site) if self.site else 'tiles'
        if datadir is None:
            datadir = '%s_%s_%s%s' % (sitename, resstr, self.dataclass.name, suffix)
        if not os.path.exists(datadir):
            os.makedirs(datadir)

        VerboseOut('Creating GIPS project %s' % datadir)
        VerboseOut('  Dates: %s' % self.datestr)
        VerboseOut('  Products: %s' % ' '.join(self.standard_products))
        for date in self.dates:
            self.data[date].project(datadir=datadir, res=res, **kwargs)
        VerboseOut('Completed GIPS project in %s' % (dt.now() - start))
        return ProjectInventory(datadir)

    def pprint(self, **kwargs):
        """ Print inventory """
        self.print_tile_coverage()
        print
        if self.site is not None:
            print Colors.BOLD + 'Asset Coverage for site %s' % basename(self.site) + Colors.OFF
        else:
            print Colors.BOLD + 'Asset Holdings' + Colors.OFF
        # Header
        #datestr = ' Month/Day' if md else ' Day of Year'
        super(DataInventory, self).pprint(**kwargs)

    def print_tile_coverage(self):
        if self.site is not None:
            print Colors.BOLD + '\nTile Coverage'
            print Colors.UNDER + '{:^8}{:>14}{:>14}'.format('Tile', '% Coverage', '% Tile Used') + Colors.OFF
            for t in sorted(self.tiles):
                print "{:>8}{:>11.1f}%{:>11.1f}%".format(t, self.tiles[t][0] * 100, self.tiles[t][1] * 100)

    @staticmethod
    def main(cls):
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(description='GIPS %s Data Utility v%s' % (cls.name, __version__))
        subparser = parser0.add_subparsers(dest='command')

        # Archive
        parser = subparser.add_parser('archive', help='Move files from current directory to data archive')
        parser.add_argument('--keep', help='Keep files after adding to archive', default=False, action='store_true')
        parser.add_argument('--recursive', help='Iterate through subdirectories', default=False, action='store_true')
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

        invparser = argparse.ArgumentParser(add_help=False, formatter_class=dhf, epilog='TEST')
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
        group.add_argument('-p', '--products', help='Requested Products (call products command to list)', nargs='*')
        extra = []
        for arg, kwargs in cls.extra_arguments().items():
            extra.append(kwargs['dest'])
            group.add_argument(arg, **kwargs)

        parents = [invparser, ]

        parser = subparser.add_parser('products', help='List available products')

        # Inventory
        parser = subparser.add_parser('inventory', help='Get Inventory', parents=parents, formatter_class=dhf)
        parser.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)
        parser.add_argument('--compact', help='Print only dates (no coverage)', default=False, action='store_true')

        # Processing
        parserp = subparser.add_parser('process', help='Process scenes', parents=parents, formatter_class=dhf)
        group = parserp.add_argument_group('Processing Options')
        # Useful option for debugging, not for general users
        #group.add_argument('--suffix', help='Suffix on end of filename (before extension)', default='')
        group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')
        group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)

        # Project
        parser = subparser.add_parser('project', help='Create project', parents=parents, formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--suffix', help='Suffix on end of project directory', default='')
        group.add_argument('--crop', help='Crop down to minimum bounding box', default=False, action='store_true')
        group.add_argument('--nowarp', help='Mosaic, but do not warp or crop', default=False, action='store_true')
        group.add_argument('--res', nargs=2, help='Resolution of (warped) output rasters', default=None, type=float)
        #group.add_argument('--datadir', help='Directory to save project files (default auto-generated)', default=None)
        group.add_argument('--format', help='Format for output file', default="GTiff")
        group.add_argument('--chunksize', help='Chunk size in MB', type=float, default=512.0)

        args = parser0.parse_args()

        VerboseOut('GIPS %s command line utility v%s' % (cls.name, __version__), 1)

        if args.command == 'products':
            cls.print_products()
            exit(1)

        gippy.Options.SetVerbose(args.verbose)
        if 'format' in args:
            gippy.Options.SetDefaultFormat(args.format)

        if args.command == 'archive':
            cls.Asset.archive(recursive=args.recursive, keep=args.keep)
            exit(1)

        kwargs = dict(zip(extra, [eval('args.%s' % a) for a in extra]))

        try:
            inv = cls.inventory(
                site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products,
                pcov=args.pcov, ptile=args.ptile, fetch=args.fetch, sensors=args.sensors, **kwargs)
            if args.command == 'inventory':
                inv.pprint(md=args.md, compact=args.compact)
            elif args.command == 'process':
                gippy.Options.SetChunkSize(args.chunksize)
                inv.process(overwrite=args.overwrite, **kwargs)
            elif args.command == 'project':
                gippy.Options.SetChunkSize(args.chunksize)
                inv.project(suffix=args.suffix, crop=args.crop, nowarp=args.nowarp, res=args.res, **kwargs)
            else:
                VerboseOut('Command %s not recognized' % args.command, 0)
        except Exception, e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
