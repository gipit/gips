#!/usr/bin/env python

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
import traceback
import tarfile
import re

from pdb import set_trace

import gippy
from gippy.utils import VerboseOut, RemoveFiles
from gippy.data.datainventory import DataInventory


class Data(object):
    """ Base class for data objects """
    # name of this dataset
    name = ''
    # dictionary of codes and full names of all sensors
    sensors = {}
    # best resolution of data (default used in Data.project)
    _defaultresolution = [30.0, 30.0]

    # root directory to data
    _rootdir = ''
    _tiledir = os.path.join(_rootdir, 'tiles')
    _stagedir = os.path.join(_rootdir, '.stage')

    # Format code of date directories in repository
    _datedir = '%Y%j'
    # pattern of created products
    _prodpattern = '*.tif'
    # pattern of a metadata or header file
    _metapattern = ''
    # shapefile name or PostGIS layer name describing sensor tiles
    _tiles_vector = os.path.join(_rootdir, 'vectors', 'tiles.shp')
    # column name in _tiles_vector holding tile designation
    _tiles_attribute = 'tile'
    # dictionary of available assets for this dataset
    _assets = {}
    # dictionary of available products for this dataset
    _products = {}

    @classmethod
    def inspect(cls, filename):
        """ Inspect a single asset and get info - Needs to be overridden by child """
        return {
            'asset': '',    # the asset code
            'filename': '',  # full filename to asset
            'datafiles': [], # if filename is archive, index of datafiles in archive
            'tile': '',      # tile designation
            'date': '',     # full date
            #'path': '',     # path in repository where asset belongs
            'basename': '',  # base/root name of this tile/date (used for product naming)
            'sensor': '',   # sensor code (key used in cls.sensors dictionary)
            'products': {}   # dictionary of existing products in asset {'product name': filename}
        }

    def meta(self, tile):
        """ Retrieve metadata for this tile """
        filename = self.tiles[tile]['filename']
        meta = self.inspect(filename)
        # add additional metadata to dictionary
        return meta

    def processtile(self, tile, products):
        """ Make sure all products exist and process if needed """
        pass

    def filter(self, tile, **kwargs):
        """ Check if tile passes filter """
        return True

    @classmethod
    def feature2tile(cls, feature):
        """ Get tile designaation from a geospatial feature (i.e. a row) """
        fldindex = feature.GetFieldIndex(cls._tiles_attribute)
        return str(feature.GetField(fldindex))

    @classmethod
    def fetch_asset(cls, asset, tile, date):
        """ Get this asset for this tile and date """
        raise Exception("Fetch not implemented for %s" % cls.name)

    @classmethod
    def asset_dates(cls, asset, tile, dates, days):
        """ For a given asset get all dates possible (in repo or not) """
        from dateutil.rrule import rrule, DAILY
        # default assumes daily regardless of asset or tile
        dates = [dt for dt in rrule(DAILY, dtstart=dates[0], until=dates[1]) if days[0] <= int(dt.strftime('%j')) <= days[1]]
        return dates

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    def find_assets(self, tile):
        """ Find assets (raw/original data) for this tile """
        assets = []
        for key in self._assets:
            files = glob.glob(os.path.join(self._tiledir, tile, self.date.strftime(self._datedir), self._assets[key]['pattern']))
            assets = assets + files
        return assets

    @classmethod
    def find_tiles(cls):
        """ Get list of all available tiles """
        return os.listdir(cls._tiledir)

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates available in repository for a tile """
        tdir = os.path.join(cls._tiledir, tile)
        if os.path.exists(tdir):
            return [datetime.strptime(os.path.basename(d), cls._datedir).date() for d in os.listdir(tdir)]
        else:
            return []

    @classmethod
    def path(cls, tile, date, filename=''):
        """ Return path in repository for this tile and date """
        return os.path.join(cls._tiledir, tile, str(date.strftime(cls._datedir)), filename)

    def opentile(self, tile, product=''):
        if product != '':
            return gippy.GeoImage(self.products[product])
        elif len(self.products) == 1:
            return gippy.GeoImage(self.products[self.products.keys()[0]])
        else:
            # return filename of a tile from self.tiles ?
            raise Exception('Invalid product %s' % product)

    def discover(self, tile):
        """ This analyzes a tile and date directory for assets and products """
        assets = self.find_assets(tile)
        if len(assets) == 0:
            raise Exception('no assets')
        dat = {'assets': [], 'datafiles': {}, 'path': self.path(tile,self.date), 'basename': '', 'products': {}}
        for a in assets:
            info = self.inspect(a)
            dat['assets'].append(info['filename'])
            dat['datafiles'][info['filename']] = info['datafiles']
            # should this be property of tile date, sensor?
            dat['basename'] = info['basename']
            dat['sensor'] = info['sensor']
            dat['products'].update(info['products'])
            # find additional products
            files = glob.glob(os.path.join(dat['path'],info['basename']+self._prodpattern))
            files2 = []
            for f in files:
                ext = os.path.splitext(f)[1]
                # bad....don't check gz, this should check that is not asset
                if ext != '.hdr' and ext != '.xml' and ext != '.gz' and ext != '.index': files2.append(f)
            files = files2
            # TODO - BASENAME IS DIFFERENT FOR EVERY ASSET !??  Because of sensor issue
            for f in files:
                fname, ext = os.path.splitext(os.path.split(f)[1])
                dat['products'][fname[len(info['basename'])+1:]] = f
        #if any(assets):
        #    VerboseOut('%s %s: assets and products found' % (tile, self.date),3)
        #    VerboseOut(assets,3)
        #    VerboseOut(dat['products'],3)
        return dat

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    @classmethod
    def inventory(cls, **kwargs):
        return DataInventory(cls, **kwargs)

    @classmethod
    def sensor_names(cls):
        """ All possible sensor names """
        return sorted(cls.sensors.values())

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
    def archive(cls, path='', recursive=False, keep=False):
        """ Move files from directory to archive location """
        start = datetime.now()

        qdir = os.path.join(cls._rootdir, 'quarantine')
        try:
            os.makedirs(qdir)
        except:
            pass

        to_remove = []
        fnames = []
        for a in cls._assets.values():
            fnames.extend(glob.glob(os.path.join(path, a['pattern'])))

        numlinks = 0
        numfiles = 0
        for f in fnames:
            try:
                info = cls.inspect(f)
            except Exception, e:
                # if problem with inspection, move to quarantine
                qname = os.path.join(qdir, os.path.basename(f))
                if not os.path.exists(qname):
                    os.link(os.path.abspath(f), qname)
                to_remove.append(f)
                VerboseOut('%s -> quarantine (file error)' % f, 2)
                VerboseOut(traceback.format_exc(), 4)
                continue
            if not hasattr(info['date'],'__len__'):
                info['date'] = [info['date']]
            for d in info['date']:
                added = 0
                bname = os.path.basename(f)
                tpath = cls.path(info['tile'], d)
                pattern = cls._assets[info['asset']]['pattern']
                newfilename = os.path.join(tpath, bname)
                if not os.path.exists(newfilename):
                    # check if another asset exists
                    existing_assets = glob.glob(os.path.join(tpath, pattern))
                    if len(existing_assets) > 0:
                        VerboseOut('%s: other version(s) already exists:' % bname, 2)
                        for ef in existing_assets: VerboseOut('\t%s' % os.path.basename(ef),2)
                    else:
                        try:
                            os.makedirs(tpath)
                        except OSError as exc:
                            if exc.errno == errno.EEXIST and os.path.isdir(tpath):
                                pass
                            else:
                                raise Exception('Unable to make data directory %s' % tpath)
                        os.link(os.path.abspath(f), newfilename)
                        #shutil.move(os.path.abspath(f),newfilename)
                        VerboseOut(bname + ' -> ' + newfilename, 2)
                        added = 1
                        numlinks = numlinks + added
                        to_remove.append(f)
                else:
                    VerboseOut('%s already in archive' % f, 2)
                    to_remove.append(f)
            numfiles = numfiles + added

        if not keep: RemoveFiles(to_remove,['.index', '.aux.xml'])
        # Summarize
        VerboseOut('%s files (%s links) from %s added to archive in %s' %
                    (numfiles, numlinks, os.path.abspath(path), datetime.now()-start) )
        if numfiles != len(fnames):
            VerboseOut('%s files not added to archive' % (len(fnames)-numfiles))

    @classmethod
    def fetch(cls, products, tiles, dates, days):
        """ Download data for tile and add to archive """
        assets = cls.products2assets(products)

        for a in assets:
            for t in tiles:

                asset_dates = cls.asset_dates(a, t, dates, days)

                for d in asset_dates:

                    if not glob.glob(cls.path(t,d,cls._assets[a]['pattern'])):

                        status = cls.fetch_asset(a,t,d)

                        print "status:", status
                        # what to do if status is nonzero?

                        # move files as you get them
                        VerboseOut('copying data to archive', 2)

                        try:
                            print "calling archive moving to %s" % cls._stagedir
                            cls.archive(cls._stagedir)
                        except:
                            VerboseOut('archive of downloaded files was unsuccessful', 2)


    @classmethod
    def products2assets(cls,products):
        """ Get list of assets needed for these products """
        assets = []
        for p in products:
            if 'assets' in cls._products[p]:
                assets.extend(cls._products[p]['assets'])
            else: assets.append('')
        return set(assets)

    def process(self, overwrite=False, suffix=''):
        """ Determines what products need to be processed for each tile and calls processtile """
        if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
        for tile, info in self.tiles.items():
            # Determine what needs to be processed
            toprocess = {}
            prods = [p for p in self.products if p in self._products.keys()]
            for p in prods:
                fout = os.path.join(info['path'],info['basename']+'_'+p+suffix)
                # need to figure out extension properly
                if len(glob.glob(fout+'*')) == 0 or overwrite: toprocess[p] = fout
            if len(toprocess) != 0:
                VerboseOut(['Processing products for tile %s' % tile, toprocess],3)
                self.processtile(tile,toprocess)

    def project(self, res=None, datadir=''):
        """ Create image of final product (reprojected/mosaiced) """
        if datadir == '': datadir = self.name+'_data'
        self.process()
        if not os.path.exists(datadir): os.makedirs(datadir)
        datadir = os.path.abspath(datadir)
        if res is None: res = self._defaultresolution
        if not hasattr(res, "__len__"): res = [res,res]
        #elif len(res) == 1: res = [res[0],res[0]]
        if self.site is None:
            raise Exception("No site file supplied")
        for product in self.products:
            if self.products[product] == '':
                start = datetime.now()
                filename = os.path.join(datadir, self.date.strftime('%Y%j') + '_%s_%s.tif' % (product,self.sensor))
                if not os.path.exists(filename):
                    filenames = [self.tiles[t]['products'][product] for t in self.tiles]
                    # cookiecutter should validate pixels in image.  Throw exception if not
                    imgout = gippy.CookieCutter(filenames, filename, self.site, res[0], res[1])
                    VerboseOut('Projected and cropped %s files -> %s in %s' % (len(filenames),imgout.Basename(),datetime.now() - start))
                self.products[product] = filename

    @classmethod
    def get_tiles_vector(cls):
        """ Get GeoVector of sensor grid """
        if os.path.isfile(cls._tiles_vector):
            tiles = gippy.GeoVector(cls._tiles_vector)
        else:
            try:
                tiles = gippy.GeoVector("PG:dbname=geodata host=congo port=5432 user=ags", layer=cls._tiles_vector)
            except:
                raise Exception('unable to access %s tiles (file or database)' % cls.name)
        return tiles

    @classmethod
    def vector2tiles(cls, vector, pcov=0.0, ptile=0.0):
        """ Return matching tiles and coverage % for provided vector """
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
            if (tiles[t][0] < pcov/100.0) or(tiles[t][1] < ptile/100.0):
                remove_tiles.append(t)
        for t in remove_tiles: tiles.pop(t,None)
        return tiles

    @classmethod
    def extracthdr(cls,filename):
        """ extract metadata header file """
        if tarfile.is_tarfile(filename):
            tfile = tarfile.open(filename)
        else: raise Exception('%s is not a valid tar file' % filename)
        index = tfile.getnames()
        dirname = os.path.dirname(filename)
        # need error handling
        hdrfname = ([f for f in index if cls._metapattern in f])[0]
        hdrfname2 = os.path.join(dirname,hdrfname)
        if not os.path.exists(os.path.join(dirname,hdrfname2)):
            tfile.extract(hdrfname,dirname)
        return hdrfname2

    @classmethod
    def extractdata(cls,filename):
        """ Extract data from original datafile, if tarfile """
        if tarfile.is_tarfile(filename):
            tfile = tarfile.open(filename)
        else: raise Exception('%s is not a valid tar file' % filename)
        index = tfile.getnames()
        dirname = os.path.dirname(filename)
        datafiles = []
        tfile.extractall(dirname)
        for f in index:
            if cls._metapattern in f:
                hdrfile = os.path.join(dirname,f)
            else: datafiles.append(os.path.join(dirname,f))
            # make sure file is readable and writable
            try:
                ff = os.path.join(dirname,f)
                if not os.path.isdir(ff):
                    os.chmod(ff,0664)
            except:
                pass
        return {'headerfile': hdrfile, 'datafiles': datafiles}

    def __init__(self, site=None, tiles=None, date=None, products=None, sensors=None, fetch=False, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a tile dictionary (see next)
        self.tiles[tile] - dictionary of tile data {path, baename, sensor, products}
        self.products - dictionary of product name and final product filename
        """
        self.site = site
        # Calculate spatial extent
        if tiles is not None:
            self.tile_coverage = dict((t,1) for t in tiles)
        elif site is not None:
            self.tile_coverage = self.vector2tiles(gippy.GeoVector(site),**kwargs)
        else:
            self.tile_coverage = dict((t,(1,1)) for t in self.find_tiles())
        self.date = date

        #VerboseOut("Locating matching data for %s" % self.date, 3)

        # Create product dictionaries for use by child class
        self.tiles = {}
        if products is None: products = self._products.keys()
        if len(products) == 0: products = self._products.keys()

        self.products = {}
        for p in products: self.products[p] = ''

        # For each tile locate files/products
        if sensors is None: sensors = self.sensors.keys()
        self.used_sensors = {s: self.sensors.get(s,None) for s in sensors}

        #VerboseOut('Finding products for %s tiles ' % (len(self.tile_coverage)),3)
        for t in self.tile_coverage.keys():
            try:
                self.tiles[t] = self.discover(t)
                # Custom filter based on dataclass
                #good = self.filter(t,filename, **kwargs)
                #if good == False:
                #    empty_tiles.append(t)
                # check all tiles - should be same sensor - MODIS
                self.sensor = self.tiles[t]['sensor']
            except:
                continue
        if len(self.tiles) == 0: raise Exception('No valid data found')

    @staticmethod
    def args_inventory():
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('inventory arguments')
        group.add_argument('-s','--site',help='Vector file for region of interest', default=None)
        group.add_argument('-t','--tiles', nargs='*', help='Tile designations', default=None)
        group.add_argument('-d','--dates',help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
        group.add_argument('--days',help='Include data within these days of year (doy1,doy2)',default=None)
        group.add_argument('--fetch',help='Fetch any missing data (if supported)',default=False, action='store_true')
        group.add_argument('-p','--products', nargs='*', help='Process/filter these products', default=None)
        group.add_argument('-v','--verbose',help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
        group.add_argument('--%cov',dest='pcov',help='Threshold of %% coverage of tile over site', default=0, type=int)
        group.add_argument('--%tile',dest='ptile',help='Threshold of %% tile used', default=0, type=int)
        return parser

    #@staticmethod
    #def add_subparsers(parser):
    #    return parser

    def __str__(self):
        return self.sensor + ': ' + str(self.date)

    @classmethod
    def main(cls):
        print 'Ian\'s working copy!'
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(description='%s Data Utility' % cls.name, formatter_class=argparse.RawTextHelpFormatter)
        subparser = parser0.add_subparsers(dest='command')

        invparser = cls.args_inventory()

        # Help
        parser = subparser.add_parser('help',help='Print extended help', formatter_class=dhf)

        # Inventory
        parser = subparser.add_parser('inventory',help='Get Inventory', parents=[invparser], formatter_class=dhf)
        parser.add_argument('--md',help='Show dates using MM-DD',action='store_true',default=False)

        # Processing
        parser = subparser.add_parser('process',help='Process scenes', parents=[invparser],formatter_class=dhf)
        group = parser.add_argument_group('Processing Options')
        group.add_argument('--overwrite', help='Overwrite output files if they exist', default=False, action='store_true')
        group.add_argument('--suffix', help='Append string to end of filename (before extension)', default='')
        #group.add_argument('--nooverviews', help='Do not add overviews to output', default=False, action='store_true')
        #pparser.add_argument('--link', help='Create links in current directory to output', default=False, action='store_true')
        #pparser.add_argument('--multi', help='Use multiple processors', default=False, action='store_true')

        # Project
        parser = subparser.add_parser('project',help='Create project', parents=[invparser], formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--res',nargs=2,help='Resolution of output rasters', default=None, type=float)
        group.add_argument('--datadir', help='Directory to save project files', default=cls.name+'_data')

        # Links
        parser = subparser.add_parser('link',help='Link to Products', parents=[invparser], formatter_class=dhf)
        parser.add_argument('--hard',help='Create hard links instead of symbolic', default=False,action='store_true')

        # Misc
        parser = subparser.add_parser('archive',help='Move files from current directory to data archive')
        #parser.add_argument('--link',help='Create symbolic links instead of moving', default=False,action='store_true')
        parser.add_argument('--keep',help='Keep files after adding to archive', default=False,action='store_true')
        parser.add_argument('--recursive',help='Iterate through subdirectories', default=False,action='store_true')
        parser.add_argument('-v','--verbose',help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

        parser = subparser.add_parser('fetch',help='Fetch products from remote location')

        #cls.add_subparsers(subparser)

        # Pull in cls options here
        #dataparser = subparser.add_parser('data',help='', parents=[invparser],formatter_class=dhf)

        args = parser0.parse_args()
        if args.command == 'help':
            parser0.print_help()
            print '\navailable products:'
            for key,val in cls._products.items():
                print '    {:<20}{:<100}'.format(key, val['description'])
            exit(1)

        gippy.Options.SetVerbose(args.verbose)
        gippy.Options.SetChunkSize(128.0)   # replace with option

        VerboseOut('GIPPY %s command line utility' % cls.name)

        try:
            if args.command == 'archive':
                if args.recursive:
                    for root, subdirs, files in os.walk('.'):
                        cls.archive(path=root, keep=args.keep)
                else:
                    cls.archive(keep=args.keep)
                exit(1)
            inv = cls.inventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles,
                products=args.products, pcov=args.pcov, ptile=args.ptile, fetch=args.fetch)
            if args.command == 'inventory':
                inv.printcalendar(args.md)
            elif args.command == 'link':
                inv.links(args.hard)
            elif args.command == 'process':
                inv.process(overwrite=args.overwrite,suffix=args.suffix) #, nooverviews=args.nooverviews)
            elif args.command == 'project':
                inv.project(args.res, datadir=args.datadir)
            else:
                VerboseOut('Command %s not recognized' % cmd)
        except Exception,e:
            VerboseOut('Error in %s: %s' % (args.command, e))
            VerboseOut(traceback.format_exc(), 4)
