#!/usr/bin/env python

import os, sys
import shutil, errno
import argparse
import ogr
from datetime import datetime
import glob
from shapely.wkb import loads
from shapely.geometry import shape
import traceback
import tarfile

from pdb import set_trace

import gippy
from gippy.utils import VerboseOut
from gippy.data.datainventory import DataInventory

class Data(object):
    """ Base class for data objects """
    # name of this dataset
    name = ''
    # dictionary of codes and full names of all sensors
    sensors = {}
    # best resolution of data (default used in Data.project)
    _defaultresolution = [30.0,30.0]

    # root directory to data
    _rootdir = ''

    # Format code of date directories in repository
    _datedir = '%Y%j'
    # pattern of original raw data file
    _pattern = ''
    # pattern of created products
    _prodpattern = '*.tif'
    # pattern of a metadata or header file
    _metapattern = ''
    # shapefile name or PostGIS layer name describing sensor tiles
    _tiles_vector = ''
    # column name in _tiles_vector holding tile designation
    _tiles_attribute = 'tile'
    # dictionary of available products for this dataset
    _products = {}
    
    @classmethod
    def inspect(cls, filename):
        """ Inspect a single file and derive some info (for archiving) - Needs to be overridden by child """
        return {
            'filename':'',  # full filename
            'datafiles':[], # if filename is archive, index of datafiles in archive
            'tile':'',      # tile designation
            'date': '',     # full date
            'basename':'',  # base/root name of this tile/date (used for product naming)
            'sensor': '',   # sensor code (key used in cls.sensors dictionary)
            'products':{}   # dictionary {'product name': filename}
    }

    def meta(self, tile):
        """ Retrieve metadata for this tile """
        filename = self.tiles[tile]['filename']
        meta = self.inspect(filename)
        # add additional metadata to dictionary
        return meta

    def processtile(self,tile,products):
        """ Make sure all products exist and process if needed """
        pass

    def filter(self, tile, **kwargs):
        """ Check if tile passes filter """
        return True

    @classmethod
    def feature2tile(cls,feature):
        """ Get tile designaation from a geospatial feature (i.e. a row) """
        fldindex = feature.GetFieldIndex(cls._tiles_attribute)
        return str(feature.GetField(fldindex))

    def fetch(self,*args,**kwargs):
        """ Download data and add to archive """
        raise Exception("Fetch not implemented for %s" % self.name)

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    def find_raw(self, tile):
        """ Find raw/original data for this tile """
        return glob.glob(os.path.join(self._rootdir, tile, self.date.strftime(self._datedir), self._pattern))

    @classmethod
    def find_tiles(cls):
        """ Get list of all available tiles """
        return os.listdir(cls._rootdir)

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates available for a tile """
        tdir = os.path.join(cls._rootdir,tile)
        if os.path.exists(tdir):
            return [datetime.strptime(os.path.basename(d),cls._datedir).date() for d in os.listdir( os.path.join(cls._rootdir,tile) )]
        else: return []

    @classmethod
    def path(cls, tile, date):
        """ Return path in repository for this tile and date """
        return os.path.join(cls._rootdir, tile, str(date.strftime(cls._datedir)))

    def opentile(self, tile, product=''):
        if product != '':
            return gippy.GeoImage(self.products[product])
        elif len(self.products) == 1:
            return gippy.GeoImage(self.products[self.products.keys()[0]])
        else:
            # return filename of a tile from self.tiles ?
            raise Exception('Invalid product %s' % product)

    def find_data(self, tile):
        """ Find all data for tile, add products dictionary to info dict from inspect function """
        # find original datafile
        filenames = self.find_raw(tile)
        if len(filenames) == 0: return {}
        # Some datasets will have more then one file
        if len(filenames) > 1:
            raise Exception('More than 1 file found for same tile/date')
        info = self.inspect(filenames[0])

        # find additional products named basename_product
        # TODO replace with regular expressions
        path = self.path(info['tile'],info['date'])
        files = glob.glob(os.path.join(path,info['basename']+self._prodpattern))
        files2 = []
        for f in files:
            ext = os.path.splitext(f)[1]
            if ext != '.hdr' and ext != '.xml': files2.append(f)
        files = files2
        products = info['products']
        for f in files:
            fname,ext = os.path.splitext(os.path.split(f)[1])
            products[ fname[len(info['basename'])+1:]  ] = f
        # extend info['products'] instead of replace ?
        info['products'] = products
        return info

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    @classmethod
    def inventory(cls, **kwargs): return DataInventory(cls, **kwargs)

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
    def archive(cls, path='', link=True):
        """ Move files from directory to archive location """
        start = datetime.now()
        fnames = glob.glob(os.path.join(path,cls._pattern))
        qdir = os.path.join(cls._rootdir,'../quarantine')
        try:
            os.makedirs(qdir)
        except:
            pass
        numlinks = 0
        numfiles = 0
        for f in fnames:
            try:
                meta = cls.inspect(f)
            except Exception,e:
                # if problem with inspection, move to quarantine
                qname = os.path.join(qdir,f)
                if not os.path.exists(qname):
                    os.link(os.path.abspath(f),os.path.join(qdir,f))
                VerboseOut('%s -> quarantine (file error)' % f,2)
                VerboseOut(traceback.format_exc(), 4)
                continue
            if not hasattr(meta['date'],'__len__'): meta['date'] = [meta['date']]
            for d in meta['date']:
                added = cls._move2archive(tile=meta['tile'],date=d,filename=f)
                numlinks = numlinks + added
            numfiles = numfiles + added
        # Summarize
        VerboseOut( '%s files (%s links) added to archive in %s' % (numfiles, numlinks, datetime.now()-start) )
        if numfiles != len(fnames):
            VerboseOut( '%s files not added to archive' % (len(fnames)-numfiles) )

    @classmethod
    def _move2archive(cls, tile, date, filename):
        path = cls.path(tile, date)
        try:
            os.makedirs(path)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise Exception('Unable to make data directory %s' % path)
        # Move or link file
        origpath,fname = os.path.split(filename)
        newfilename = os.path.join(path,fname)
        if not os.path.exists(newfilename):
            # Check for older versions
            existing_files = glob.glob(os.path.join(path,cls._pattern))
            if len(existing_files) > 0:
                VerboseOut('Other version of %s already exists:' % fname,2)
                for ef in existing_files: VerboseOut('\t%s' % os.path.basename(ef),2)
            else:
                os.link(os.path.abspath(filename),newfilename)
                #shutil.move(os.path.abspath(f),newfilename)
            VerboseOut(fname + ' -> ' + newfilename,2)
            return 1
        else:
            VerboseOut('%s already in archive' % fname, 2)
        return 0

    def process(self, overwrite=False, suffix=''):
        """ Determines what products need to be processed for each tile and calls processtile """
        if suffix != '' and suffix[:1] != '_': suffix = '_' + suffix
        for tile, info in self.tiles.items():
            # Determine what needs to be processed
            path = self.path(tile,info['date'])
            toprocess = {}
            prods = [p for p in self.products if p in self._products.keys()]
            for p in prods:
                fout = os.path.join(path,info['basename']+'_'+p+suffix)
                # need to figure out extension properly
                if len(glob.glob(fout+'*')) == 0 or overwrite: toprocess[p] = fout
            if len(toprocess) != 0:
                VerboseOut(['Processing products for tile %s' % tile, toprocess],3)
                self.processtile(tile,toprocess)

    def project(self, res=None, datadir='gipdata'):
        """ Create image of final product (reprojected/mosaiced) """
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
            tiles = gippy.GeoVector("PG:dbname=geodata host=congo port=5432 user=ags", layer=cls._tiles_vector)
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

        if fetch: self.fetch()
        # Find products
        for t in self.tile_coverage.keys():
            tile = self.find_data(t)
            # Custom filter based on dataclass
            #good = self.filter(t,filename, **kwargs)
            #if good == False:
            #    empty_tiles.append(t)
            if any(tile):
                self.tiles[t] = tile
        self.sensor = self.tiles[self.tiles.keys()[0]]['sensor']
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
        group.add_argument('--suffix', help='Append string to end of filename (before extension)',default='')
        #group.add_argument('--nooverviews', help='Do not add overviews to output', default=False, action='store_true')
        #pparser.add_argument('--link', help='Create links in current directory to output', default=False, action='store_true')
        #pparser.add_argument('--multi', help='Use multiple processors', default=False, action='store_true')

        # Project
        parser = subparser.add_parser('project',help='Create project', parents=[invparser], formatter_class=dhf)
        group = parser.add_argument_group('Project options')
        group.add_argument('--res',nargs=2,help='Resolution of output rasters', default=None, type=float)
        group.add_argument('--datadir', help='Directory to save project files', default='gipdata')

        # Links
        parser = subparser.add_parser('link',help='Link to Products', parents=[invparser], formatter_class=dhf)
        parser.add_argument('--hard',help='Create hard links instead of symbolic', default=False,action='store_true')

        # Misc
        parser = subparser.add_parser('archive',help='Move files from current directory to data archive')
        parser.add_argument('--link',help='Create symbolic links instead of moving', default=False,action='store_true')
        parser.add_argument('-v','--verbose',help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

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
        gippy.Options.SetChunkSize(256.0)   # replace with option

        VerboseOut('GIPPY %s command line utility' % cls.name)

        try:
            if args.command == 'archive':
                cls.archive(link=args.link)
                exit(1)
            inv = cls.inventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, 
                products=args.products, pcov=args.pcov, ptile=args.ptile)
            if args.fetch: inv.fetch()
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