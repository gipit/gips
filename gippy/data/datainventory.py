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
import argparse
import gippy
from gippy.utils import VerboseOut, parse_dates

from datetime import datetime
import traceback

from pdb import set_trace


class Tiles(object):
    """ Base class for groups of tiles """
    # Classes for assets and tiles of this Data
    #Tile = Tile

    def __init__(self, dataclass, site=None, tiles=None, date=None, products=None,
                 suffix='', sensors=None, fetch=False, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a Tile instance
        self.products - dictionary of product name and final product filename
        """

        self.site = site
        # Calculate spatial extent
        if tiles is not None:
            self.tile_coverage = dict((t, 1) for t in tiles)
        elif site is not None:
            self.tile_coverage = self.Repository.vector2tiles(gippy.GeoVector(site), **kwargs)
        else:
            self.tile_coverage = dict((t, (1, 1)) for t in self.Repository.find_tiles())
        self.date = date

        #VerboseOut("Locating matching data for %s" % self.date, 3)

        # Create product dictionary of requested products and filename
        self.products = {}
        for p in products:
            self.products[p] = [''] + products[p]
        self.suffix = suffix

        # For each tile locate files/products
        if sensors is None:
            sensors = dataclass.Asset._sensors.keys()
        self.used_sensors = {s: dataclass.Asset._sensors.get(s, None) for s in sensors}

        # TODO - expand verbose text: tiles, date, etc.
        #VerboseOut('Finding products for %s tiles ' % (len(self.tile_coverage)), 4)
        self.tiles = {}
        for t in self.tile_coverage.keys():
            #VerboseOut("Tile %s" % t, 4)
            try:
                tile = dataclass(t, self.date)
                # Custom filter based on dataclass
                #good = self.filter(t,filename, **kwargs)
                #if good == False:
                #    empty_tiles.append(t)
                self.tiles[t] = tile
                # check all tiles - should be same sensor - MODIS?
                self.sensor = tile.sensor
            except:
                VerboseOut(traceback.format_exc(), 5)
                continue
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    def open(self, product='', update=True):
        """ Open and return final product GeoImage """
        if product == '':
            # default to initial product?
            product = self.products.keys()[0]
        fname = self.products[product][0]
        if os.path.exists(fname):
            return gippy.GeoImage(fname, update)
        else:
            raise Exception('%s product does not exist' % product)

    def process(self, overwrite=False):
        """ Determines what products need to be processed for each tile and calls processtile """
        if self.suffix != '':
            self.suffix = '_' + self.suffix if self.suffix[0] != '_' else self.suffix
        for tileid, tile in self.tiles.items():
            # Determine what needs to be processed
            toprocess = {}
            for p in self.products:
                fout = os.path.join(tile.path, tile.basename+'_'+p)
                print fout
                for i in range(1, len(self.products[p])):
                    fout = fout + '_' + self.products[p][i]
                fout = fout + self.suffix
                # need to figure out extension properly
                # TODO - this is after inventory, just check the inventory for matching filename?
                if len(glob.glob(fout+'*')) == 0 or overwrite:
                    toprocess[p] = self.products[p]
                    toprocess[p][0] = fout
                    # TODO - this needs tile name
            if len(toprocess) != 0:
                VerboseOut(['Processing products for tile %s' % tileid, toprocess], 3)
                self.tiles[tileid].process(toprocess)

    def project(self, res=None, datadir=''):
        """ Create image of final product (reprojected/mosaiced) """
        if datadir == '':
            datadir = self.name+'_data'
        self.process()
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        datadir = os.path.abspath(datadir)
        if res is None:
            res = self._defaultresolution
        if not hasattr(res, "__len__"):
            res = [res, res]
        #elif len(res) == 1: res = [res[0],res[0]]
        if self.site is None:
            # Link files instead
            for t in self.tiles:
                self.tiles[t].link(products=self.products, path=datadir)
        else:
            # TODO - better file naming
            for product in self.products:
                if self.products[product][0] == '':
                    start = datetime.now()
                    sitebase = os.path.splitext(os.path.basename(self.site))[0] + '_'
                    bname = sitebase + self.date.strftime('%Y%j') + '_%s_%s.tif' % (product, self.sensor)
                    filename = os.path.join(datadir, bname)
                    if not os.path.exists(filename):
                        filenames = [self.tiles[t].products[product] for t in self.tiles]
                        # cookiecutter should validate pixels in image.  Throw exception if not
                        imgout = gippy.CookieCutter(filenames, filename, self.site, res[0], res[1])
                        VerboseOut('Projected and cropped %s files -> %s in %s' % (len(filenames),
                                   imgout.Basename(), datetime.now() - start))
                    self.products[product][0] = filename


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

    def __init__(self, dataclass, site=None, tiles=None, dates=None, days=None, products=None, fetch=False, **kwargs):
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
            self.tiles = Repository.vector2tiles(gippy.GeoVector(self.site), **kwargs)

        self.temporal_extent(dates, days)

        self.data = {}
        self.products = products
        #if self.products is None:
        #    self.products = dataclass.Tile._products.keys()
        #if len(self.products) == 0:
        #    self.products = dataclass.Tile._products.keys()

        if fetch and products is not None:
            dataclass.fetch(products, self.tiles, (self.start_date, self.end_date), (self.start_day, self.end_day))
            Repository.archive(Repository.spath())

        # get all potential matching dates for tiles
        dates = []
        for t in self.tiles:
            #VerboseOut('locating matching dates', 5)
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
                dat = Tiles(dataclass=dataclass, site=self.site, tiles=self.tiles.keys(), date=date, products=self.products, **kwargs)
                self.data[date] = dat
                self.numfiles = self.numfiles + len(dat.tiles)
            except Exception, e:
                VerboseOut('Inventory error %s' % e, 3)
                VerboseOut(traceback.format_exc(), 4)

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
        if self.products is None:
            raise Exception('No products specified for processing')
        start = datetime.now()
        VerboseOut('Requested %s products for %s files' % (len(self.products), self.numfiles))
        for date in self.dates:
            self.data[date].process(*args, **kwargs)
        VerboseOut('Completed processing in %s' % (datetime.now()-start))

    def project(self, *args, **kwargs):
        start = datetime.now()
        VerboseOut('Projecting data for %s dates (%s - %s)' % (len(self.dates), self.dates[0], self.dates[-1]))
        # res should default to data?
        for date in self.dates:
            self.data[date].project(*args, **kwargs)
        VerboseOut('Completed projecting in %s' % (datetime.now()-start))

    # TODO - check if this is needed
    def get_products(self, date):
        """ Get list of products for given date """
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
        img = self.data[self.dates[0]][0].open(product=product)
        for i in range(1, len(self.dates)):
            img.AddBand(self.data[self.dates[i]][0].open(product=product)[0])
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

        parents = [invparser, cls.arg_parser()]

        # Inventory
        parser = subparser.add_parser('inventory', help='Get Inventory', parents=parents, formatter_class=dhf)
        parser.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)
        parser.add_argument('-p', '--products', help='Show products', default=False, action='store_true')

        # Processing
        parserp = subparser.add_parser('process', help='Process scenes', parents=parents, formatter_class=dhf)
        group = parserp.add_argument_group('Processing Options')
        group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')

        # Project
        parser = subparser.add_parser('project', help='Create project', parents=parents, formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--res', nargs=2, help='Resolution of output rasters', default=None, type=float)
        group.add_argument('--datadir', help='Directory to save project files', default=cls.name+'_data')
        group.add_argument('--format', help='Format for output file', default="GTiff")

        args = parser0.parse_args()

        gippy.Options.SetVerbose(args.verbose)
        # TODO - replace with option
        gippy.Options.SetChunkSize(128.0)
        if 'format' in args:
            gippy.Options.SetDefaultFormat(args.format)

        VerboseOut('GIPPY %s command line utility' % cls.name)

        if args.command == 'archive':
            # TODO - take in path argument
            cls.Asset.archive(recursive=args.recursive, keep=args.keep)
            exit(1)

        #products = [p for p in cls.Tile._products if eval('args.%s' % p) not in [None, False]]
        products = {}
        for p in cls._products:
            if p != '':
                val = eval('args.%s' % p)
                if val not in [None, False]:
                    products[p] = val

        try:
            inv = cls.inventory(
                site=args.site, dates=args.dates, days=args.days, tiles=args.tiles,
                products=products, pcov=args.pcov, ptile=args.ptile, fetch=args.fetch, suffix=args.suffix)
            if args.command == 'inventory':
                inv.printcalendar(args.md, products=args.products)
            elif args.command == 'link':
                inv.links(args.hard)
            elif args.command == 'process':
                inv.process(overwrite=args.overwrite)
            elif args.command == 'project':
                inv.project(args.res, datadir=args.datadir)
            else:
                VerboseOut('Command %s not recognized' % cmd)
        except Exception, e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
