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
import shutil
import textwrap
from pprint import pprint

import gippy
from gipif.utils import VerboseOut, parse_dates, Colors
from gipif.GeoVector import GeoVector
import commands
import tempfile
from gipif.version import __version__
from pdb import set_trace


def project_inventory(datadir=''):
    """ Inventory a directory of project files """
    files = glob.glob(os.path.join(datadir, '*.tif'))
    inv = {}
    for f in files:
        basename = os.path.splitext(os.path.basename(f))[0]
        parts = basename.split('_')
        date = datetime.strptime(parts[0], '%Y%j').date()
        sensor = parts[1]
        product = basename[len(parts[0])+len(parts[1])+2:]
        if date in inv.keys():
            inv[date][product] = f
        else:
            inv[date] = {product: f}
    return inv


class Tiles(object):
    """ Collection of tiles for a single date """

    def __init__(self, dataclass, site=None, tiles=None, date=None,
                 products=None, sensors=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a Tile instance
        self.products - dictionary of product name and final product filename
        """
        self.dataclass = dataclass
        self.site = site
        # Calculate spatial extent
        if isinstance(tiles, list):
            tiles = dict((t, 1) for t in tiles)
        self.tile_coverage = tiles
        #if tiles is not None:
        #    self.tile_coverage = dict((t, 1) for t in tiles)
        #elif site is not None:
        #    self.tile_coverage = self.Repository.vector2tiles(GeoVector(site), **kwargs)
        #else:
        #    self.tile_coverage = dict((t, (1, 1)) for t in self.Repository.find_tiles())
        self.date = date
        self.products = {}
        self.requested_products = products
        if sensors is None:
            sensors = dataclass.Asset._sensors.keys()

        # For each tile locate files/products
        VerboseOut('%s: searching %s tiles for products and assets' % (self.date, len(self.tile_coverage)), 4)
        self.tiles = {}
        for t in self.tile_coverage.keys():
            try:
                tile = dataclass(t, self.date)
                # Custom filter based on dataclass
                good = tile.filter(**kwargs)
                if good and tile.sensor in sensors:
                    self.tiles[t] = tile
                # TODO - check all tiles - should be same sensor?
                self.sensor = tile.sensor
            except:
                #VerboseOut(traceback.format_exc(), 5)
                continue
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    def coverage(self):
        """ Calculates % coverage of site for each asset """
        asset_coverage = {}
        for a in self.dataclass.Asset._assets:
            cov = 0.0
            norm = float(len(self.tile_coverage)) if self.site is None else 1.0
            for t in self.tiles:
                if a in self.tiles[t].assets:
                    cov = cov + (self.tile_coverage[t][0]/norm)
            asset_coverage[a] = cov*100
        return asset_coverage

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

    def project(self, res=None, datadir='', nowarp=False):
        """ Create image of final product (reprojected/mosaiced) """
        if res is None:
            res = self.dataclass.Asset._defaultresolution
        if not hasattr(res, "__len__"):
            res = [res, res]
        start = datetime.now()
        bname = self.date.strftime('%Y%j')
        sensor = self.sensor if self.sensor != '' else ''
        if self.site is None:
            for t in self.tiles:
                datadir = self._datadir(t) if datadir == '' else datadir
                for p in self.requested_products:
                    fout = os.path.join(datadir, bname + ('_%s_%s.tif' % (sensor, p)))
                    if not os.path.exists(fout):
                        try:
                            VerboseOut("Creating %s" % os.path.basename(fout))
                            shutil.copy(self.tiles[t].products[p], fout)
                        except:
                            VerboseOut("Problem copying %s" % filename, 3)
        else:
            sitename = os.path.splitext(os.path.basename(self.site))[0]
            resstr = '_%sx%s' % (res[0], res[1])
            datadir = self._datadir(sitename + resstr) if datadir == '' else datadir
            for product in self.requested_products:
                filename = os.path.join(datadir, bname + ('_%s_%s.tif' % (sensor, product)))
                if not os.path.exists(filename):
                    try:
                        filenames = [self.tiles[t].products[product] for t in self.tiles]
                        # TODO - cookiecutter should validate pixels in image.  Throw exception if not
                        if nowarp:
                            imgout = self._mosaic(filenames, filename, self.site)
                        else:
                            imgout = gippy.CookieCutter(filenames, filename, self.site, res[0], res[1])
                        imgout = None
                    except:
                        VerboseOut("Problem projecting %s" % filename, 3)
                self.products[product] = filename
        t = datetime.now() - start
        VerboseOut('%s: created project files for %s tiles in %s' % (self.date, len(self.tiles), t), 2)

    def _datadir(self, name):
        """ Determine data directory and created if needed """
        if name != '':
            name = '_' + name
        datadir = '%sData%s' % (self.dataclass.name, name)
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        return os.path.abspath(datadir)

    def _mosaic(self, infiles, outfile, vectorfile):
        """ Mosaic multiple files together, but do not warp """
        img = gippy.GeoImage(infiles[0])
        nd = img[0].NoDataValue()
        srs = img.Projection()
        img = None
        for f in range(1, len(infiles)):
            img = gippy.GeoImage(infiles[f])
            if img.Projection() != srs:
                raise Exception("Input files have non-matching projections and must be warped")
        # transform vector to image projection
        vector = GeoVector(vectorfile)
        vsrs = vector.proj()
        from gipif.GeoVector import transform_shape
        geom = transform_shape(vector.union(), vsrs, srs)
        extent = geom.bounds
        ullr = "%f %f %f %f" % (extent[0], extent[3], extent[2], extent[1])
        # run command
        nodatastr = '-n %s -a_nodata %s -init %s' % (nd, nd, nd)
        cmd = 'gdal_merge.py -o %s -ul_lr %s %s %s' % (outfile, ullr, nodatastr, " ".join(infiles))
        result = commands.getstatusoutput(cmd)
        imgout = gippy.GeoImage(outfile)
        # warp and rasterize vector
        vec1 = vector.transform(srs)
        vec1name = vec1.filename
        td = tempfile.mkdtemp()
        mask = gippy.GeoImage(os.path.join(td, vector.layer.GetName()), imgout, gippy.GDT_Byte, 1)
        maskname = mask.Filename()
        mask = None
        cmd = 'gdal_rasterize -at -burn 1 -l %s %s %s' % (vec1.layer.GetName(), vec1.filename, maskname)
        result = commands.getstatusoutput(cmd)
        mask = gippy.GeoImage(maskname)
        imgout.AddMask(mask[0]).Process().ClearMasks()
        vec1 = None
        mask = None
        shutil.rmtree(os.path.dirname(maskname))
        shutil.rmtree(os.path.dirname(vec1name))
        return imgout

    def open(self, product='', update=True):
        """ Open and return final product GeoImage """
        fname = self.products[product]
        if os.path.exists(fname):
            print 'open', product, fname
            return gippy.GeoImage(fname, update)
        else:
            raise Exception('%s product does not exist' % product)

    def print_assets(self, dformat='%j', color=''):
        """ Print coverage for each and every asset """
        #assets = [a for a in self.dataclass.Asset._assets]
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        asset_coverage = self.coverage()
        for a in sorted(asset_coverage):
            sys.stdout.write(color + '  {:>4.1f}%   '.format(asset_coverage[a]) + Colors.OFF)
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
            sys.stdout.write('  '+p)
        sys.stdout.write('\n')

    #def print_products(self, dformat='%j'):
    #    """ Print products that have been processed """
    #    sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))


class DataInventory(object):
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

    @property
    def dates(self):
        """ Get sorted list of dates """
        return [k for k in sorted(self.data)]

    def __getitem__(self, date):
        return self.data[date]

    def get_timeseries(self, product=''):
        """ Read all files as time series """
        filenames = [self.data[date].products[product] for date in self.dates]
        img = gippy.GeoImage(filenames)
        return img

    def process(self, **kwargs):
        """ Process data in inventory """
        if len(self.standard_products) + len(self.composite_products) == 0:
            raise Exception('No products specified!')
        sz = self.numfiles
        if len(self.standard_products) > 0:
            start = datetime.now()
            VerboseOut('Processing %s files: %s' % (sz, ' '.join(self.standard_products)), 1)
            for date in self.dates:
                try:
                    self.data[date].process(**kwargs)
                except:
                    pass
            VerboseOut('Completed processing in %s' % (datetime.now()-start), 1)
        if len(self.composite_products) > 0:
            start = datetime.now()
            VerboseOut('Processing %s files into composites: %s' % (sz, ' '.join(self.composite_products)), 1)
            self.dataclass.process_composites(self, self.composite_products, **kwargs)
            VerboseOut('Completed processing in %s' % (datetime.now()-start), 1)

    def project(self, *args, **kwargs):
        """ Create project files for data in inventory """
        self.process(**kwargs)
        start = datetime.now()
        pstr = ' '.join(self.standard_products)  # .join(self.composite_products)
        dstr = '%s dates (%s - %s)' % (len(self.dates), self.dates[0], self.dates[-1])
        VerboseOut('Creating project files (%s) for %s' % (pstr, dstr), 1)
        for date in self.dates:
            self.data[date].project(*args, **kwargs)
        VerboseOut('Completed creating project files in %s' % (datetime.now()-start), 1)

    def print_inv(self, md=False, compact=False):
        """ Print inventory """
        self.print_tile_coverage()
        print
        if self.site is not None:
            site_name = os.path.splitext(os.path.basename(self.site))[0]
            print Colors.BOLD + 'Asset Coverage for site %s' % site_name + Colors.OFF
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
        parser0 = argparse.ArgumentParser(description='GIPIF %s Data Utility v%s' % (cls.name, __version__),
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
        group.add_argument('--nowarp', help='Mosaic, but do not warp', default=False, action='store_true')
        group.add_argument('--res', nargs=2, help='Resolution of (warped) output rasters', default=None, type=float)
        group.add_argument('--datadir', help='Directory to save project files (default auto-generated)', default='')
        group.add_argument('--format', help='Format for output file', default="GTiff")
        group.add_argument('--chunksize', help='Chunk size in MB', type=float, default=512.0)

        args = parser0.parse_args()

        gippy.Options.SetVerbose(args.verbose)
        if 'format' in args:
            gippy.Options.SetDefaultFormat(args.format)

        VerboseOut('GIPIF %s command line utility v%s' % (cls.name, __version__), 1)

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
                inv.project(args.res, datadir=args.datadir, nowarp=args.nowarp)
            else:
                VerboseOut('Command %s not recognized' % cmd, 0)
        except Exception, e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
