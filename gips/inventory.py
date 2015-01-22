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
import argparse
from datetime import datetime as dt
import traceback
import numpy
from copy import deepcopy

import gippy
from gips import __version__, SpatialExtent, TemporalExtent
from gips.tiles import Tiles
from gips.utils import VerboseOut, Colors, basename, map_reduce
from gips.data.core import Data


class Inventory(object):
    """ Base class for inventories """
    _colors = [Colors.PURPLE, Colors.RED, Colors.GREEN, Colors.BLUE]

    def __init__(self):
        pass

    def __getitem__(self, date):
        """ Indexing operator for class """
        return self.data[date]

    def get_subset(self, dates):
        """ Return subset of inventory """
        inv = deepcopy(self)
        for d in inv.dates:
            if d not in dates:
                del inv.data[d]
        return inv

    @property
    def sensor_set(self):
        sset = set()
        for date in self.dates:
            sset.update(self.data[date].sensor_set)
        return sorted(sset)

    @property
    def dates(self):
        """ Get sorted list of dates """
        return sorted(self.data.keys())

    @property
    def numfiles(self):
        """ Total number of files in inventory """
        return sum([len(dat) for dat in self.data.values()])

    @property
    def datestr(self):
        return '%s dates (%s - %s)' % (len(self.dates), self.dates[0], self.dates[-1])

    def color(self, sensor):
        """ Return color for sensor """
        return self._colors[list(self.sensor_set).index(sensor)]

    def pprint(self, md=False, compact=False):
        print self.data[self.data.keys()[0]].pprint_header()
        dformat = '%m-%d' if md else '%j'
        oldyear = 0
        formatstr = '\n{:<12}' if compact else '{:<12}\n'
        colors = {k: self.color(k) for k in self.sensor_set}
        for date in self.dates:
            if date.year != oldyear:
                sys.stdout.write(Colors.BOLD + formatstr.format(date.year) + Colors.OFF)
            if compact:
                color = self.color(self.data[date].sensor_set[0])
                dstr = color + ('{:^%s}' % (7 if md else 4)).format(date.strftime(dformat)) + Colors.OFF
                sys.stdout.write(dstr)
            else:
                self.data[date].pprint(dformat, colors)
            oldyear = date.year
        if self.numfiles != 0:
            VerboseOut("\n\n%s files on %s dates" % (self.numfiles, len(self.dates)), 1)
            self.print_legend()

    def print_legend(self):
        print Colors.BOLD + '\nSENSORS' + Colors.OFF
        for key in sorted(self.sensor_set):
            try:
                desc = self.dataclass.Asset._sensors[key]['description']
                scode = key + ': ' if key != '' else ''
            except:
                desc = ''
                scode = key
            print self.color(key) + '%s%s' % (scode, desc) + Colors.OFF


class ProjectInventory(Inventory):
    """ Inventory of project directory (collection of Data class) """

    def __init__(self, projdir='', products=[]):
        """ Create inventory of a GIPS project directory """
        self.projdir = os.path.abspath(projdir)
        if not os.path.exists(self.projdir):
            raise Exception('Directory %s does not exist!' % self.projdir)

        self.data = {}
        product_set = set()
        sensor_set = set()
        try:
            for dat in Data.discover(self.projdir):
                self.data[dat.date] = dat
                # All products and sensors used across all dates
                product_set = product_set.union(dat.product_set)
                sensor_set = sensor_set.union(dat.sensor_set)

            if not products:
                products = list(product_set)
            self.requested_products = products
            self.sensors = sensor_set
        except:
            VerboseOut(traceback.format_exc(), 4)
            raise Exception("%s does not appear to be a GIPS project directory" % self.projdir)

    def products(self, date):
        """ Intersection of available products and requested products for this date """
        return set(self.data[date].products).intersection(set(self.requested_products))

    def new_image(self, filename, dtype=gippy.GDT_Byte, numbands=1, nodata=None):
        """ Create new image with the same template as the files in project """
        img = gippy.GeoImage(self.data[self.dates[0]].open(self.requested_products[0]))
        imgout = gippy.GeoImage(filename, img, dtype, numbands)
        img = None
        if nodata is not None:
            imgout.SetNoData(nodata)
        return imgout

    def data_size(self):
        """ Get 'shape' of inventory: #products x rows x columns """
        img = gippy.GeoImage(self.data[self.dates[0]].open(self.requested_products[0]))
        sz = (len(self.requested_products), img.YSize(), img.XSize())
        return sz

    def get_data(self, dates=None, products=None, chunk=0):
        """ Read all files as time series, stacking all products """
        # TODO - change to absolute dates
        days = numpy.array([int(d.strftime('%j')) for d in dates])
        imgarr = []
        if dates is None:
            dates = self.dates
        if products is None:
            products = self.requested_products
        for p in products:
            gimg = self.get_timeseries(p, dates=dates)
            # TODO - move numpy.squeeze into swig interface file?
            arr = numpy.squeeze(gimg.TimeSeries(days.astype('float64'), chunk))
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

    def map_reduce(self, func, numbands=1, products=None, readfunc=None, **kwargs):
        """ Apply func to inventory to generate an image with numdim output bands """
        if products is None:
            products = self.requested_products
        if readfunc is None:
            readfunc = lambda x: self.get_data(products=products, chunk=x)
        sz = self.data_size()
        return map_reduce(sz[1:3], readfunc, func, numbands=numbands, **kwargs)


class DataInventory(Inventory):
    """ Manager class for data inventories (collection of Tiles class) """

    def __init__(self, dataclass, site=None, tiles=None, dates=None, days=None,
                 products=None, fetch=False, pcov=0.0, ptile=0.0, **kwargs):
        """ Create a new inventory
        :dataclass: The Data class to use (e.g., LandsatData, ModisData)
        :site: The site shapefile
        :tiles: List of tile ids
        :dates: tuple of begin and end date
        :days: tuple of begin and end day of year
        :products: List of requested products of interest
        :fetch: bool indicated if missing data should be downloaded
        """
        self.dataclass = dataclass
        Repository = dataclass.Asset.Repository

        try:
            self.spatial = SpatialExtent(dataclass, site, tiles, pcov, ptile)
            self.temporal = TemporalExtent(dates, days)
            self.products = dataclass.RequestedProducts(products)
        except Exception, e:
            import traceback
            VerboseOut(traceback.format_exc(), 4)
            raise Exception('Illformed parameters: %s' % e)

        if fetch:
            try:
                dates = self.temporal.datebounds
                days = self.temporal.daybounds
                dataclass.fetch(self.products.base, self.spatial.tiles, dates, days)
            except Exception:
                #VerboseOut(traceback.format_exc(), 3)
                raise Exception('DataInventory: Error downloading')
            dataclass.Asset.archive(Repository.spath())
        self.data = {}
        for date in self.temporal.prune_dates(self.spatial.available_dates):
            try:
                dat = Tiles(dataclass, self.spatial, date, self.products, **kwargs)
                if len(dat) > 0:
                    self.data[date] = dat
            except Exception, e:
                raise Exception('DataInventory: Error accessing tiles in %s repository' % dataclass.name)

        if len(self.data) == 0:
            raise Exception("No matching files in inventory!")

    @property
    def sensor_set(self):
        return sorted(self.dataclass.Asset._sensors.keys())

    def process(self, **kwargs):
        """ Process data in inventory """
        if len(self.products.standard) + len(self.products.composite) == 0:
            raise Exception('No products specified!')
        sz = self.numfiles
        if len(self.products.standard) > 0:
            start = dt.now()
            VerboseOut('Processing %s files: %s' % (sz, ' '.join(self.products.standard)), 1)
            for date in self.dates:
                try:
                    self.data[date].process(**kwargs)
                except:
                    VerboseOut(traceback.format_exc(), 3)
                    pass
            VerboseOut('Completed processing in %s' % (dt.now() - start), 1)
        if len(self.products.composite) > 0:
            start = dt.now()
            VerboseOut('Processing %s files into composites: %s' % (sz, ' '.join(self.products.composite)), 1)
            self.dataclass.process_composites(self, self.products.composite, **kwargs)
            VerboseOut('Completed processing in %s' % (dt.now() - start), 1)

    def project(self, datadir=None, suffix='', res=None, nomosaic=False, **kwargs):
        """ Create project files for data in inventory """
        self.process(**kwargs)
        start = dt.now()
        sitename = 'tiles' if self.spatial.site is None else basename(self.spatial.site)
        if res is None:
            res = self.dataclass.Asset._defaultresolution

        suffix = '' if suffix is None else '_' + suffix
        if datadir is None:
            # formulate project directory name
            if res[0] == res[1]:
                resstr = str(res[0])
            else:
                resstr = '%sx%s' % (res[0], res[1])
            datadir = '%s_%s_%s%s' % (sitename, resstr, self.dataclass.name, suffix)
        if self.spatial.site is None or nomosaic:
            datadir = os.path.join(datadir, "TILEID")

        VerboseOut('Creating GIPS project %s' % datadir)
        VerboseOut('  Dates: %s' % self.datestr)
        VerboseOut('  Products: %s' % ' '.join(self.products.standard))

        [self.data[d].project(datadir=datadir, res=res, nomosaic=nomosaic, **kwargs) for d in self.dates]

        VerboseOut('Completed GIPS project in %s' % (dt.now() - start))
        if not nomosaic and self.spatial.site is not None:
            inv = ProjectInventory(datadir)
            inv.pprint()
            return inv

    def pprint(self, **kwargs):
        """ Print inventory """
        print
        if self.spatial.site is not None:
            print Colors.BOLD + 'Asset Coverage for site %s' % basename(self.spatial.site) + Colors.OFF
            self.spatial.print_tile_coverage()
            print
        else:
            print Colors.BOLD + 'Asset Holdings' + Colors.OFF
        # Header
        #datestr = ' Month/Day' if md else ' Day of Year'
        super(DataInventory, self).pprint(**kwargs)

    @staticmethod
    def main(cls):
        """ Command line intrepreter """
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(description='GIPS %s Data Utility v%s' % (cls.name, __version__))
        subparser = parser0.add_subparsers(dest='command')

        # Archive
        parser = subparser.add_parser('archive', help='Move files from current directory to data archive')
        parser.add_argument('--keep', help='Keep files after adding to archive', default=False, action='store_true')
        parser.add_argument('--recursive', help='Iterate through subdirectories', default=False, action='store_true')
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

        invparser = argparse.ArgumentParser(add_help=False, formatter_class=dhf)
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
        group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')
        group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)

        # Project
        parser = subparser.add_parser('project', help='Create project', parents=parents, formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--res', nargs=2, help='Resolution of (warped) output rasters', default=None, type=float)
        group.add_argument('--crop', help='Crop down to minimum bounding box', default=False, action='store_true')
        group.add_argument('--nowarp', help='Do not warp (or crop)', default=False, action='store_true')
        group.add_argument('--interpolation', help='Interpolate using: 0-NN, 1-Bilinear, 2-Cubic', default=0, type=int)
        group.add_argument('--nomosaic', help='Do not mosaic (keep as tiles)', default=False, action='store_true')
        group.add_argument('--datadir', help='Directory to store output', default=None)
        group.add_argument('--suffix', help='Suffix on end of project directory', default=None)
        group.add_argument('--format', help='Format for output file', default="GTiff")
        group.add_argument('--chunksize', help='Chunk size in MB', type=float, default=512.0)

        args = parser0.parse_args()

        VerboseOut(Colors.BOLD + 'GIPS %s utility v%s' % (cls.name, __version__) + Colors.OFF, 1)

        if args.command == 'products':
            cls.print_products()
            exit(1)

        # Set GIPPY options
        gippy.Options.SetVerbose(args.verbose)
        if 'format' in args:
            gippy.Options.SetDefaultFormat(args.format)
        if 'chunksize' in args:
            gippy.Options.SetChunkSize(args.chunksize)

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
                inv.process(overwrite=args.overwrite, **kwargs)
            elif args.command == 'project':
                inv.project(datadir=args.datadir, suffix=args.suffix, res=args.res, crop=args.crop,
                            nowarp=args.nowarp, interpolation=args.interpolation, nomosaic=args.nomosaic, **kwargs)
            else:
                VerboseOut('Command %s not recognized' % args.command, 0)
        except Exception, e:
            VerboseOut(traceback.format_exc(), 4)
            VerboseOut('Error in %s: %s' % (args.command, e))
