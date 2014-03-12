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
import sys
import shutil
import errno
import argparse
import ogr
from datetime import datetime
import glob
from shapely.wkb import loads
from shapely.geometry import shape
import tarfile
import traceback
from pdb import set_trace

import gippy
from gippy.utils import VerboseOut, RemoveFiles, File2List, List2File
from gippy.data.datainventory import DataInventory
from gippy.settings import DATABASES


class Asset(object):
    """ Class for a single file asset """
    # root directory to data
    _rootpath = ''
    # Format code of date directories in repository
    _datedir = '%Y%j'
    # Sensors
    _sensors = {
        '': {'description': ''},
    }
    # dictionary of assets
    _assets = {
        '': {
            'pattern': '*',
        }
    }

    def __init__(self, filename):
        """ Inspect a single file and populate variables. Needs to be extended """
        _tilespath = os.path.join(self._rootpath, self._tilesdir)
        _stagepath = os.path.join(self._rootpath, self._stagedir)
        _qpath = os.path.join(self._rootpath, self._qdir)
        # full filename to asset
        self.filename = filename
        # the asset code
        self.asset = ''
        # if filename is archive, index of datafiles in archive...needed?
        #self.datafiles = []
        # tile designation
        self.tile = ''
        # full date
        self.date = datetime.now()
         # base/root name of this tile/date (used for product naming)
        self.basename = os.path.basename(filename)
        # sensor code (key used in cls.sensors dictionary)
        self.sensor = ''
        # dictionary of existing products in asset {'product name': filename}
        self.products = {}

    @classmethod
    def asset_dates(cls, asset, tile, dates, days):
        """ For a given asset get all dates possible (in repo or not) - used for fetch """
        from dateutil.rrule import rrule, DAILY
        # default assumes daily regardless of asset or tile
        datearr = rrule(DAILY, dtstart=dates[0], until=dates[1])
        dates = [dt for dt in datearr if days[0] <= int(dt.strftime('%j')) <= days[1]]
        return dates

    @classmethod
    def fetch(cls, asset, tile, date):
        """ Get this asset for this tile and date """
        raise Exception("Fetch not implemented for %s" % cls.__name__)

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    @property
    def path(self):
        """ Return repository path to this asset """
        return os.path.join(self._tilespath, self.tile, str(date.strftime(cls._datedir)))

    @classmethod
    def path(cls, tile, date):
        return os.path.join(cls._rootpath, cls._tilesdir, tile, str(date.strftime(cls._datedir)))

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    _tilesdir = 'tiles'
    _qdir = 'quarantine'
    _stagedir = 'stage'
    _vectordir = 'vectors'

    def sensor_meta(self):
        return self._sensors[self.sensor]

    @classmethod
    def find_assets(cls, tile, date, asset=None):
        """ Find assets in this path and return list of Asset objects """
        path = cls.path(tile, date)
        if asset is not None:
            assets = [asset]
        else:
            assets = cls._assets.keys()
        found = []
        for a in assets:
            files = glob.glob(os.path.join(path, cls._assets[a]['pattern']))
            # more than 1 asset??
            if len(files) > 1:
                VerboseOut(files, 3)
                raise Exception("Duplicate(?) assets found")
            if len(files) == 1:
                found.append(cls(files[0]))
        return found

    @classmethod
    def archive(cls, path='', recursive=False, keep=False):
        """ Move assets from directory to archive location """
        start = datetime.now()

        try:
            os.makedirs(cls._qpath)
        except:
            pass

        if recursirve:
            for root, subdirs, files in os.walk(path):
                fnames = []
                for a in cls._assets.values():
                    fnames.extend(glob.glob(os.path.join(root, a['pattern'])))
        else:
            fnames = glob.glob(os.path.join(path, a['pattern']))

        numlinks = 0
        numfiles = 0
        for f in fnames:
            links = cls._archivefile(f)
            if links >= 0:
                numfiles = numfiles + 1
                numlinks = numlinks + added
                if not keep:
                    RemoveFiles(f, ['.index', '.aux.xml'])

        # Summarize
        VerboseOut('%s files (%s links) from %s added to archive in %s' %
                  (numfiles, numlinks, os.path.abspath(path), datetime.now()-start))
        if numfiles != len(fnames):
            VerboseOut('%s files not added to archive' % (len(fnames)-numfiles))

    @classmethod
    def _archivefile(cls, filename):
        """ archive specific file """
        bname = os.path.basename(filename)
        try:
            asset = cls(filename)
        except Exception, e:
            # if problem with inspection, move to quarantine
            qname = os.path.join(qdir, bname)
            if not os.path.exists(qname):
                os.link(os.path.abspath(filename), qname)
            VerboseOut('%s -> quarantine (file error)' % filename, 2)
            VerboseOut(traceback.format_exc(), 4)
            return 0

        if not hasattr(asset.date, '__len__'):
            asset.date = asset.date
        numlinks = 0
        otherversions = False
        for d in asset.date:
            tpath = cls.path(asset.tile, d)
            newfilename = os.path.join(tpath, bname)
            if not os.path.exists(newfilename):
                # check if another asset exists
                existing = cls.find_assets(asset.tile, d, asset.asset)
                if len(existing) > 0:
                    VerboseOut('%s: other version(s) already exists:' % bname, 2)
                    for ef in existing:
                        VerboseOut('\t%s' % os.path.basename(ef), 2)
                    otherversions = True
                else:
                    try:
                        os.makedirs(tpath)
                    except OSError as exc:
                        if exc.errno == errno.EEXIST and os.path.isdir(tpath):
                            pass
                        else:
                            raise Exception('Unable to make data directory %s' % tpath)
                    os.link(os.path.abspath(filename), newfilename)
                    #shutil.move(os.path.abspath(f),newfilename)
                    VerboseOut(bname + ' -> ' + newfilename, 2)
                    numlinks = numlinks + 1
            else:
                VerboseOut('%s already in archive' % filename, 2)
        if otherversions and numlinks == 0:
            return -1
        else:
            return numlinks
        # should return asset instance

    def datafiles(self):
        """ Get list of datafiles from asset (if archive file) """
        if tarfile.is_tarfile(self.filename):
            tfile = tarfile.open(self.filename)
        else:
            raise Exception('%s is not a valid tar file' % self.filename)
        path = os.path.dirname(self.filename)
        indexfile = os.path.join(path, self.filename+'.index')
        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:
            tfile = tarfile.open(self.filename)
            datafiles = tfile.getnames()
            List2File(datafiles, indexfile)
        return datafiles

    def extract(self, filenames=[]):
        """ Extract filenames from asset (if archive file) """
        if tarfile.is_tarfile(self.filename):
            tfile = tarfile.open(self.filename)
        else:
            raise Exception('%s is not a valid tar file' % self.filename)
        path = os.path.dirname(self.filename)
        if len(filenames) == 0:
            filenames = self.datafiles()
        extracted_files = []
        for f in filenames:
            fname = os.path.join(path, f)
            if not os.path.exists(fname):
                tfile.extract(f, path)
            try:
                # this ensures we have permissions on extracted files
                if not os.path.isdir(fname):
                    os.chmod(fname, 0664)
            except:
                pass
            extracted_files.append(fname)
        return extracted_files

    def __str__(self):
        return os.path.basename(self.filename)


class Tile(object):
    """ A tile of data (multiple assets/products) on a single date """

    # pattern of created products
    _prodpattern = '*.tif'
    # dictionary of available products for this dataset
    _products = {
    }

    Asset = Asset

    def process(self, products):
        """ Make sure all products exist and process if needed """
        pass

    def meta(self):
        """ Retrieve metadata for this tile """
        print '%s metadata!' % self.__name__
        #meta = self.Asset(filename)
        # add metadata to dictionary
        return {}

    #def filter(self, **kwargs):
    #    """ Check if tile passes filter """
    #    return True

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    @property
    def path(self):
        """ Return repository path to this tile dir """
        return os.path.join(self.Asset._rootpath, self.Asset._tilesdir,
                            self.id, str(self.date.strftime(self.Asset._datedir)))

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    def open(self, product=''):
        if product != '':
            return gippy.GeoImage(self.products[product])
        elif len(self.products) == 1:
            return gippy.GeoImage(self.products[self.products.keys()[0]])
        else:
            # return filename of a tile from self.tiles ?
            raise Exception('Invalid product %s' % product)

    def link(self, products, path='', copy=False):
        """ Create links in path to tile products """
        for p in products:
            fname = self.products[p]
            bname = os.path.basename(fname)
            fullbname = os.path.join(path, bname)
            if copy:
                try:
                    # TODO - copying doesn't seem to work
                    os.copy(fname, fullbname)
                    VerboseOut('%s: copying' % bname, 2)
                    return
                except:
                    VerboseOut('%s: Problem copying file' % bname, 2)
            # try hard link first, if it fails, soft link
            try:
                os.link(fname, fullbname)
                VerboseOut('%s: hard linking' % bname, 2)
            except:
                try:
                    os.symlink(fname, fullbname)
                    VerboseOut('%s: soft linking' % bname, 2)
                except:
                    VerboseOut('%s: Problem creating link' % bname, 2)

    @classmethod
    def products2assets(cls, products):
        """ Get list of assets needed for these products """
        assets = []
        for p in products:
            if 'assets' in cls._products[p]:
                assets.extend(cls._products[p]['assets'])
            else:
                assets.append('')
        return set(assets)

    def __init__(self, tile, date):
        """ Find all data and assets for this tile and date """
        self.id = tile
        self.date = date
        self.assets = {}
        #self.basename = ''
        self.products = {}
        for asset in self.Asset.find_assets(tile, date):
            self.assets[asset.asset] = asset
            # sensor and basename assumes same value every time ?
            self.sensor = asset.sensor
            # should this be property of tile date, sensor?  same for all assets ?
            self.basename = asset.basename
            # products that come automatically with assets
            self.products.update(asset.products)

            # find additional products...for each asset ?
            files = glob.glob(os.path.join(self.path, asset.basename+self._prodpattern))
            files2 = []
            for f in files:
                ext = os.path.splitext(f)[1]
                # bad....don't check gz, this should check that is not asset
                exts = ['.hdr', '.xml', 'gz', '.index']
                if ext != exts:
                    files2.append(f)
            files = files2
            # TODO - BASENAME IS DIFFERENT FOR EVERY ASSET !??  Because of sensor issue
            for f in files:
                fname, ext = os.path.splitext(os.path.split(f)[1])
                self.products[fname[len(asset.basename)+1:]] = f

        if len(self.assets) == 0:
            raise Exception('no assets')

        VerboseOut('%s %s: assets and products found' % (tile, date), 5)
        VerboseOut(self.assets, 5)
        VerboseOut(self.products, 5)


class Data(object):
    """ Base class for data objects """
    # name of this dataset
    name = ''

    # best resolution of data (default used in Data.project)
    # better off in Asset...per asset , or per product?
    _defaultresolution = [30.0, 30.0]

    # Classes for assets and tiles of this Data
    Tile = Tile

    # vector of tiling/grid system
    _tiles_vector = 'tiles.shp'
    # column name in _tiles_vector holding tile designation
    _vectoratt = 'tile'

    @classmethod
    def feature2tile(cls, feature):
        """ Get tile designaation from a geospatial feature (i.e. a row) """
        fldindex = feature.GetFieldIndex(cls._vectoratt)
        return str(feature.GetField(fldindex))

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    @classmethod
    def find_tiles(cls):
        """ Get list of all available tiles """
        return os.listdir(os.path.join(cls.Tile.Asset._rootpath, cls.Tile.Asset._tilesdir))

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates available in repository for a tile """
        # TODO - ugly getting path info from cls.Tile.Asset
        tdir = os.path.join(cls.Tile.Asset._rootpath, cls.Tile.Asset._tilesdir, tile)
        if os.path.exists(tdir):
            return [datetime.strptime(os.path.basename(d), cls.Tile.Asset._datedir).date() for d in os.listdir(tdir)]
        else:
            return []

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    @classmethod
    def inventory(cls, **kwargs):
        return DataInventory(cls, **kwargs)

    # TODO - not sure if this is needed
    @classmethod
    def sensor_names(cls):
        """ All possible sensor names """
        return sorted([s['description'] for s in cls.Tile.Asset._sensors.values()])

    def open(self, product='', update=True):
        """ Open and return final product GeoImage """
        if product != '':
            return gippy.GeoImage(self.products[product], update)
        elif len(self.products) == 1:
            return gippy.GeoImage(self.products[self.products.keys()[0]], update)
        else:
            # return filename of a tile from self.tiles ?
            raise Exception('No product provided')

    @classmethod
    def fetch(cls, products, tiles, dates, days):
        """ Download data for tiles and add to archive """
        assets = cls.Tile.products2assets(products)
        for a in assets:
            for t in tiles:
                asset_dates = cls.Tile.Asset.asset_dates(a, t, dates, days)
                for d in asset_dates:
                    if not find_assets(t, d, a):
                        status = cls.Tile.Asset.fetch(a, t, d)
                        VerboseOut("Fetch status: %s" % status, 2)
                        # what to do if status is nonzero?

                        # move files as you get them
                        VerboseOut('Copying data to archive', 2)
                        try:
                            print "calling archive moving to %s" % cls._stagedir
                            cls.archive(cls._stagedir)
                        except:
                            VerboseOut('archive of downloaded files was unsuccessful', 2)

    def process(self, overwrite=False, suffix=''):
        """ Determines what products need to be processed for each tile and calls processtile """
        if suffix != '' and suffix[:1] != '_':
            suffix = '_' + suffix
        for tileid, tile in self.tiles.items():
            # Determine what needs to be processed
            toprocess = {}
            prods = [p for p in self.products if p in self.Tile._products.keys()]
            for p in prods:
                fout = os.path.join(tile.path, tile.basename+'_'+p+suffix)
                # need to figure out extension properly
                if len(glob.glob(fout+'*')) == 0 or overwrite:
                    toprocess[p] = fout
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
            for product in self.products:
                if self.products[product] == '':
                    start = datetime.now()
                    filename = os.path.join(datadir, self.date.strftime('%Y%j') + '_%s_%s.tif' % (product, self.sensor))
                    if not os.path.exists(filename):
                        filenames = [self.tiles[t].products[product] for t in self.tiles]
                        # cookiecutter should validate pixels in image.  Throw exception if not
                        imgout = gippy.CookieCutter(filenames, filename, self.site, res[0], res[1])
                        VerboseOut('Projected and cropped %s files -> %s in %s' % (len(filenames),
                                   imgout.Basename(), datetime.now() - start))
                    self.products[product] = filename

    @classmethod
    def get_tiles_vector(cls):
        """ Get GeoVector of sensor grid """
        fname = os.path.join(cls.Tile.Asset._rootpath, cls.Tile.Asset._vectordir, cls._tiles_vector)
        if os.path.isfile(fname):
            VerboseOut('%s: tiles vector %s' % (cls.__name__, fname), 4)
            tiles = gippy.GeoVector(fname)
        else:
            try:
                VerboseOut('%s: tiles vector %s' % (cls.__name__, cls._tiles_vector), 4)
                db = DATABASES['tiles']
                dbstr = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                        (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
                print dbstr
                tiles = gippy.GeoVector(dbstr, layer=cls._tiles_vector)
            except:
                raise Exception('unable to access %s tiles (file or database)' % cls.__name__)
        return tiles

    @classmethod
    def vector2tiles(cls, vector, pcov=0.0, ptile=0.0):
        """ Return matching tiles and coverage % for provided vector """
        start = datetime.now()
        import osr
        geom = vector.union()
        ogrgeom = ogr.CreateGeometryFromWkb(geom.wkb)
        tvector = cls.get_tiles_vector()
        tlayer = tvector.layer
        trans = osr.CoordinateTransformation(vector.layer.GetSpatialRef(), tlayer.GetSpatialRef())
        ogrgeom.Transform(trans)
        geom = loads(ogrgeom.ExportToWkb())
        tlayer.SetSpatialFilter(ogrgeom)
        tiles = {}
        tlayer.ResetReading()
        feat = tlayer.GetNextFeature()
        #fldindex = feat.GetFieldIndex(cls._tiles_attribute)
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkb())
            area = geom.intersection(tgeom).area
            if area != 0:
                tile = cls.feature2tile(feat)
                tiles[tile] = (area/geom.area, area/tgeom.area)
            feat = tlayer.GetNextFeature()
        remove_tiles = []
        for t in tiles:
            if (tiles[t][0] < pcov/100.0) or (tiles[t][1] < ptile/100.0):
                remove_tiles.append(t)
        for t in remove_tiles:
            tiles.pop(t, None)
        VerboseOut('%s: vector2tiles completed in %s' % (cls.__name__, datetime.now() - start), 4)
        return tiles

    def __init__(self, site=None, tiles=None, date=None, products=None, sensors=None, fetch=False, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a Tile instance
        self.products - dictionary of product name and final product filename
        """

        # shapefile name or PostGIS layer name describing sensor tiles
        #self.tiles_vector = os.path.join(self.Tile.Asset._rootpath, self.Tile.Asset._vectordir, self._tiles_vector)

        self.site = site
        # Calculate spatial extent
        if tiles is not None:
            self.tile_coverage = dict((t, 1) for t in tiles)
        elif site is not None:
            self.tile_coverage = self.vector2tiles(gippy.GeoVector(site), **kwargs)
        else:
            self.tile_coverage = dict((t, (1, 1)) for t in self.find_tiles())
        self.date = date

        #VerboseOut("Locating matching data for %s" % self.date, 3)

        # Create product dictionaries for use by child class
        self.tiles = {}
        if products is None:
            products = self.Tile._products.keys()
        if len(products) == 0:
            products = self.Tile._products.keys()

        self.products = {}
        for p in products:
            self.products[p] = ''

        # For each tile locate files/products
        if sensors is None:
            sensors = self.Tile.Asset._sensors.keys()
        self.used_sensors = {s: self.Tile.Asset._sensors.get(s, None) for s in sensors}

        VerboseOut('Finding products for %s tiles ' % (len(self.tile_coverage)), 4)
        for t in self.tile_coverage.keys():
            VerboseOut("Tile %s" % t, 4)
            try:
                tile = self.Tile(t, self.date)
                # Custom filter based on dataclass
                #good = self.filter(t,filename, **kwargs)
                #if good == False:
                #    empty_tiles.append(t)
                self.tiles[t] = tile
                # check all tiles - should be same sensor - MODIS?
                self.sensor = tile.sensor
            except:
                #print traceback.format_exc()
                continue
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    @staticmethod
    def args_inventory():
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('inventory arguments')
        group.add_argument('-s', '--site', help='Vector file for region of interest', default=None)
        group.add_argument('-t', '--tiles', nargs='*', help='Tile designations', default=None)
        group.add_argument('-d', '--dates', help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
        group.add_argument('--days', help='Include data within these days of year (doy1,doy2)', default=None)
        group.add_argument('--fetch', help='Fetch any missing data (if supported)', default=False, action='store_true')
        group.add_argument('-p', '--products', nargs='*', help='Process/filter these products', default=None)
        group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
        group.add_argument('--%cov', dest='pcov', help='Threshold of %% tile coverage over site', default=0, type=int)
        group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)
        return parser

    #@staticmethod
    #def add_subparsers(parser):
    #    return parser

    def __str__(self):
        return self.sensor + ': ' + str(self.date)

    @classmethod
    def main(cls):
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(description='%s Data Utility' % cls.name,
                                          formatter_class=argparse.RawTextHelpFormatter)
        subparser = parser0.add_subparsers(dest='command')

        parser = subparser.add_parser('help', help='Print extended help', formatter_class=dhf)

        invparser = cls.args_inventory()

        # Inventory
        parser = subparser.add_parser('inventory', help='Get Inventory', parents=[invparser], formatter_class=dhf)
        parser.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)

        # Processing
        parser = subparser.add_parser('process', help='Process scenes', parents=[invparser], formatter_class=dhf)
        group = parser.add_argument_group('Processing Options')
        group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')
        group.add_argument('--suffix', help='Append string to end of filename (before extension)', default='')
        #group.add_argument('--nooverviews', help='Do not add overviews to output', default=False, action='store_true')

        # Project
        parser = subparser.add_parser('project', help='Create project', parents=[invparser], formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--res', nargs=2, help='Resolution of output rasters', default=None, type=float)
        group.add_argument('--datadir', help='Directory to save project files', default=cls.name+'_data')

        # Misc
        parser = subparser.add_parser('archive', help='Move files from current directory to data archive')
        parser.add_argument('--keep', help='Keep files after adding to archive', default=False, action='store_true')
        parser.add_argument('--recursive', help='Iterate through subdirectories', default=False, action='store_true')
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

        #cls.add_subparsers(subparser)
        # Pull in cls options here
        #dataparser = subparser.add_parser('data',help='', parents=[invparser],formatter_class=dhf)

        args = parser0.parse_args()
        if args.command == 'help':
            parser0.print_help()
            print '\navailable products:'
            for key, val in cls.Tile._products.items():
                print '    {:<20}{:<100}'.format(key, val['description'])
            exit(1)

        gippy.Options.SetVerbose(args.verbose)
        # TODO - replace with option
        gippy.Options.SetChunkSize(128.0)

        VerboseOut('GIPPY %s command line utility' % cls.name)

        try:
            if args.command == 'archive':
                # TODO - take in path argument
                cls.archive(recursive=args.recursive, keep=args.keep)
                exit(1)
            inv = cls.inventory(
                site=args.site, dates=args.dates, days=args.days, tiles=args.tiles,
                products=args.products, pcov=args.pcov, ptile=args.ptile, fetch=args.fetch)
            if args.command == 'inventory':
                inv.printcalendar(args.md)
            elif args.command == 'link':
                inv.links(args.hard)
            elif args.command == 'process':
                inv.process(overwrite=args.overwrite, suffix=args.suffix)
            elif args.command == 'project':
                inv.project(args.res, datadir=args.datadir)
            else:
                VerboseOut('Command %s not recognized' % cmd)
        except Exception, e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
