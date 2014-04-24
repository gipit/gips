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
from datetime import datetime
import traceback

import gippy
from gipif.utils import VerboseOut, parse_dates
from gipif.GeoVector import GeoVector


class Tiles(object):
    """ Collection of tiles for a single date """

    def __init__(self, dataclass, site=None, tiles=None, date=None, products=None,
                 sensors=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a Tile instance
        self.products - dictionary of product name and final product filename
        """
        self.dataclass = dataclass
        self.site = site
        # Calculate spatial extent
        if tiles is not None:
            self.tile_coverage = dict((t, 1) for t in tiles)
        elif site is not None:
            self.tile_coverage = self.Repository.vector2tiles(GeoVector(site), **kwargs)
        else:
            self.tile_coverage = dict((t, (1, 1)) for t in self.Repository.find_tiles())
        self.date = date
        self.products = {}
        self.requested_products = products
        #for p in products:
        #    key = p
        #    for arg in products[p]:
        #        key = key + '_' + arg
        #    self.products[key] = ''

        # For each tile locate files/products
        if sensors is None:
            sensors = dataclass.Asset._sensors.keys()
        self.used_sensors = {s: dataclass.Asset._sensors.get(s, None) for s in sensors}

        # TODO - expand verbose text: tiles, date, etc.
        VerboseOut('%s: searching %s tiles for products and assets' % (self.date, len(self.tile_coverage)), 4)
        self.tiles = {}
        for t in self.tile_coverage.keys():
            try:
                tile = dataclass(t, self.date)
                # Custom filter based on dataclass
                good = tile.filter(**kwargs)
                if good and tile.sensor in sensors:
                    self.tiles[t] = tile
                # check all tiles - should be same sensor - MODIS?
                self.sensor = tile.sensor
            except:
                VerboseOut(traceback.format_exc(), 5)
                continue
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    def process(self, overwrite=False):
        """ Determines what products need to be processed for each tile and calls Data.process """
        for tileid, tile in self.tiles.items():
            toprocess = {}
            for pname, args in self.requested_products.items():
                if pname not in tile.products or overwrite:
                    toprocess[pname] = args
            if len(toprocess) != 0:
                VerboseOut('Processing products for tile %s: %s' % (tileid, ' '.join(toprocess.keys())), 2)
                self.tiles[tileid].process(toprocess)

    def project(self, res=None, datadir='', mask=None, nowarp=False):
        """ Create image of final product (reprojected/mosaiced) """
        if datadir == '':
            datadir = self.dataclass.name+'_data'
        self.process()
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        datadir = os.path.abspath(datadir)
        if res is None:
            res = self.dataclass.Asset._defaultresolution
        if not hasattr(res, "__len__"):
            res = [res, res]
        #elif len(res) == 1: res = [res[0],res[0]]
        start = datetime.now()
        if self.site is None:
            for t in self.tiles:
                self.tiles[t].link(products=self.requested_products.keys(), path=datadir, copy=True if mask else False)
                filenames = [self.tiles[t].products[p] for p in self.requested_products]
                if mask is not None:
                    self._applymask(filenames, self.tiles[t].products[mask])
        else:
            bname = os.path.splitext(os.path.basename(self.site))[0] + '_' + self.date.strftime('%Y%j')
            sensor = self.sensor if self.sensor != '' else ''
            for product in self.requested_products:
                filename = os.path.join(datadir, bname + ('_%s_%s.tif' % (sensor, product)))
                if not os.path.exists(filename):
                    filenames = [self.tiles[t].products[product] for t in self.tiles]
                    # TODO - cookiecutter should validate pixels in image.  Throw exception if not
                    if nowarp:
                        imgout = self._mosaic(filenames, filename, self.site)
                    else:
                        imgout = gippy.CookieCutter(filenames, filename, self.site, res[0], res[1])
                    imgout = None
                self.products[product] = filename
            if mask is not None:
                self._applymask(self.products.values(), self.products[mask])
        t = datetime.now() - start
        VerboseOut('%s: created project files for %s tiles in %s' % (self.date, len(self.tiles), t), 3)

    def _mosaic(self, infiles, outfile, vectorfile):
        """ Mosaic multple files together, but do not warp """
        from osgeo import gdal, osr
        import fiona
        from fiona.crs import to_string
        from pyproj import Proj, transform
        COMMAND = 'gdal_merge.py -o %s -ul_lr %s %s'
        fp = gdal.Open(infiles[0])
        crs = osr.SpatialReference()
        crs.ImportFromWkt(fp.GetProjection())
        crs = crs.ExportToProj4()
        with fiona.open(vectorfile, 'r') as source:
            proj_in = Proj(to_string(source.crs))
            proj_out = Proj(crs)
            for row in source:
                assert row['geometry']['type'] == "Polygon"
                xs = []
                ys = []
                for ring in row['geometry']['coordinates']:
                    x, y = transform(proj_in, proj_out, *zip(*ring))
                    xs.extend(x)
                    ys.extend(y)
            ullr = "%f %f %f %f" % (min(xs), max(ys), max(xs), min(ys))
            infiles = " ".join(infiles)
            command = COMMAND % (outfile, ullr, infiles)
            result = commands.getstatusoutput(command)
        return gippy.GeoImage(outfile)

    def _applymask(self, filenames, mask):
        mimg = gippy.GeoImage(mask)
        for f in filenames:
            if f != mask:
                img = gippy.GeoImage(f)
                img.AddMask(mimg[0]).Process()
                img = None
        mimg = None

    def open(self, product='', update=True):
        """ Open and return final product GeoImage """
        fname = self.products[product]
        if os.path.exists(fname):
            print 'open', product, fname
            return gippy.GeoImage(fname, update)
        else:
            raise Exception('%s product does not exist' % product)


class DataInventory(object):
    """ Manager class for data inventories """
    # redo color, combine into ordered dictionary
    _colororder = ['purple', 'bright red', 'bright green', 'bright blue', 'bright purple']
    _colorcodes = {
        'bright yellow':   '1;33',
        'bright red':      '1;31',
        'bright green':    '1;32',
        'bright blue':     '1;34',
        'bright purple':   '1;35',
        'bright cyan':     '1;36',
        'red':             '0;31',
        'green':           '0;32',
        'blue':            '0;34',
        'cyan':            '0;36',
        'yellow':          '0;33',
        'purple':          '0;35',
    }

    def _colorize(self, txt, color):
        return "\033["+self._colorcodes[color]+'m' + txt + "\033[0m"

    @property
    def dates(self):
        """ Get sorted list of dates """
        return [k for k in sorted(self.data)]

    @property
    def numdates(self):
        """ Get number of dates """
        return len(self.data)

    def __getitem__(self, date):
        return self.data[date]

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

        self.temporal_extent(dates, days)

        self.data = {}
        # Create product dictionary of requested products and filename
        if isinstance(products, list):
            self.requested_products = dict([p, [p]] for p in products)
        else:
            self.requested_products = products

        # if no products specified and only 1 product available, use it
        if len(self.requested_products) == 0 and len(dataclass._products) == 1:
            p = dataclass._products.keys()[0]
            self.requested_products = {p: [p]}
        #if self.products is None:
        #    self.products = dataclass.Tile._products.keys()
        #if len(self.products) == 0:
        #    self.products = dataclass.Tile._products.keys()

        if fetch:
            products = [val[0] for val in self.requested_products.values()]
            dataclass.fetch(products, self.tiles, (self.start_date, self.end_date), (self.start_day, self.end_day))
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
                VerboseOut(traceback.format_exc(), 4)

        self.numfiles = 0
        for date in sorted(dates):
            try:
                dat = Tiles(dataclass=dataclass, site=self.site, tiles=self.tiles.keys(),
                            date=date, products=self.requested_products, **kwargs)
                self.data[date] = dat
                self.numfiles = self.numfiles + len(dat.tiles)
            except Exception, e:
                VerboseOut(traceback.format_exc(), 4)
                #VerboseOut('Inventory error %s' % e)

    def temporal_extent(self, dates, days):
        """ Temporal extent (define self.dates and self.days) """
        if dates is None:
            dates = '1984,2050'
        self.start_date, self.end_date = parse_dates(dates)
        if days:
            days = days.split(',')
        else:
            days = (1, 366)
        self.start_day, self.end_day = (int(days[0]), int(days[1]))

    def process(self, *args, **kwargs):
        """ Process data in inventory """
        if self.requested_products is None:
            raise Exception('No products specified for processing')
        start = datetime.now()
        VerboseOut('Requested products (%s) for %s files' % (' '.join(self.requested_products), self.numfiles))
        for date in self.dates:
            self.data[date].process(*args, **kwargs)
        VerboseOut('Completed processing in %s' % (datetime.now()-start))

    def project(self, *args, **kwargs):
        start = datetime.now()
        pstr = ' '.join(self.requested_products)
        dstr = '%s dates (%s - %s)' % (len(self.dates), self.dates[0], self.dates[-1])
        VerboseOut('Creating project files (%s) for %s' % (pstr, dstr))
        # res should default to data?
        for date in self.dates:
            self.data[date].project(*args, **kwargs)
        VerboseOut('Completed creating project files in %s' % (datetime.now()-start))

    # TODO - check if this is needed
    def get_products(self, date):
        # Get list of products for given date
        # this doesn't handle different tiles (if prod exists for one tile, it lists it)
        prods = []
        dat = self.data[date]
        for t in dat.tiles:
            for p in dat.tiles[t].products:
                prods.append(p)
            #for prod in data.products.keys(): prods.append(prod)
        return sorted(set(prods))

    def printcalendar(self, md=False, products=False):
        """ print calendar for raw original datafiles """
        #import calendar
        #cal = calendar.TextCalendar()
        oldyear = ''

        # print tile coverage
        if self.site is not None:
            print '\nTILE COVERAGE'
            print '{:^8}{:>14}{:>14}'.format('Tile', '% Coverage', '% Tile Used')
            for t in sorted(self.tiles):
                print "{:>8}{:>11.1f}%{:>11.1f}%".format(t, self.tiles[t][0]*100, self.tiles[t][1]*100)
        # print inventory
        print '\nINVENTORY'
        for date in self.dates:
            dat = self.data[date]
            if md:
                daystr = str(date.month) + '-' + str(date.day)
            else:
                daystr = str(date.timetuple().tm_yday)
                if len(daystr) == 1:
                    daystr = '00' + daystr
                elif len(daystr) == 2:
                    daystr = '0' + daystr
            if date.year != oldyear:
                sys.stdout.write('\n{:>5}: '.format(date.year))
                if products:
                    sys.stdout.write('\n ')
            colors = {}
            for i, s in enumerate(self.dataclass.Asset.sensor_names()):
                colors[s] = self._colororder[i]
            #for dat in self.data[date]:
            col = colors[self.dataclass.Asset._sensors[dat.sensor]['description']]

            sys.stdout.write(self._colorize('{:<6}'.format(daystr), col))
            if products:
                sys.stdout.write('        ')
                #prods = [p for p in self.get_products(date) if p.split('_')[0] in self.products]
                #for p in prods:
                for p in self.get_products(date):
                    sys.stdout.write(self._colorize('{:<12}'.format(p), col))
                sys.stdout.write('\n ')
            oldyear = date.year
        sys.stdout.write('\n')
        if self.numfiles != 0:
            VerboseOut("\n%s files on %s dates" % (self.numfiles, self.numdates))
            self.legend()
        else:
            VerboseOut('No matching files')

    def legend(self):
        print '\nSENSORS'
        #sensors = sorted(self.dataclass.Tile.Asset._sensors.values())
        for i, s in enumerate(self.dataclass.Asset.sensor_names()):
            print self._colorize(s, self._colororder[i])
            #print self._colorize(self.dataclass.sensors[s], self._colororder[s])

    def get_timeseries(self, product=''):
        """ Read all files as time series """
        # assumes only one sensor row for each date
        filenames = [self.data[date].products[product] for date in self.dates]
        img = gippy.GeoImage(filenames)
        return img

    @staticmethod
    def main(cls):
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(description='%s Data Utility' % cls.name,
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
        group.add_argument('--fetch', help='Fetch any missing data (if supported)', default=False, action='store_true')
        #group.add_argument('-p', '--products', nargs='*', help='Process/filter these products', default=None)
        group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
        group.add_argument('--%cov', dest='pcov', help='Threshold of %% tile coverage over site', default=0, type=int)
        group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)
        group.add_argument('--suffix', help='Suffix on end of filename (before extension)', default='')
        extra = []
        for arg, kwargs in cls.extra_arguments().items():
            extra.append(kwargs['dest'])
            group.add_argument(arg, **kwargs)

        parents = [invparser, cls.arg_parser()]

        # Inventory
        parser = subparser.add_parser('inventory', help='Get Inventory', parents=parents, formatter_class=dhf)
        parser.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)
        parser.add_argument('-p', '--products', help='Show products', default=False, action='store_true')

        # Processing
        parserp = subparser.add_parser('process', help='Process scenes', parents=parents, formatter_class=dhf)
        group = parserp.add_argument_group('Processing Options')
        group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')
        group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)

        # Project
        parser = subparser.add_parser('project', help='Create project', parents=parents, formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--res', nargs=2, help='Resolution of output rasters', default=None, type=float)
        group.add_argument('--nowarp', help='Mosaic, but do not warp to site', default=False, action='store_true')
        group.add_argument('--mask', nargs='?', help='Apply this product to all products', const='acca')
        group.add_argument('--datadir', help='Directory to save project files', default=cls.name+'_data')
        group.add_argument('--format', help='Format for output file', default="GTiff")
        group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)

        args = parser0.parse_args()

        gippy.Options.SetVerbose(args.verbose)
        if 'format' in args:
            gippy.Options.SetDefaultFormat(args.format)

        VerboseOut('GIPPY %s command line utility' % cls.name)

        if args.command == 'archive':
            # TODO - take in path argument
            cls.Asset.archive(recursive=args.recursive, keep=args.keep)
            exit(1)

        suffix = '_' + args.suffix if args.suffix != '' else ''
        #products = [p for p in cls.Tile._products if eval('args.%s' % p) not in [None, False]]
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
        if args.command == 'project':
            if args.mask:
                m = args.mask.split('_')
                if len(m) == 1:
                    products[args.mask] = m
                else:
                    products[args.mask] = m
        #print 'Requested Products: ', products
        kwargs = dict(zip(extra, [eval('args.%s' % a) for a in extra]))

        try:
            inv = cls.inventory(
                site=args.site, dates=args.dates, days=args.days, tiles=args.tiles,
                products=products, pcov=args.pcov, ptile=args.ptile, fetch=args.fetch, **kwargs)
            if args.command == 'inventory':
                inv.printcalendar(args.md, products=args.products)
            elif args.command == 'process':
                gippy.Options.SetChunkSize(args.chunksize)
                inv.process(overwrite=args.overwrite)
            elif args.command == 'project':
                gippy.Options.SetChunkSize(args.chunksize)
                inv.project(args.res, datadir=args.datadir, mask=args.mask, nowarp=args.nowarp)
            else:
                VerboseOut('Command %s not recognized' % cmd)
        except Exception, e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
