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
import gippy
import agspy.utils.dateparse as dateparse
from pdb import set_trace
import pprint

def VerboseOut(obj, level=1):
    if gippy.Options.Verbose() >= level: 
        #pprint.PrettyPrinter().pprint(obj)
        if not isinstance(obj,(list,tuple)): obj = [obj]
        for o in obj: print o

def File2List(filename):
    f = open(filename)
    txt = f.readlines()
    txt2 = []
    for t in txt: txt2.append( t.rstrip('\n') )
    return txt2

def List2File(lst,filename):
    f = open(filename,'w')
    f.write('\n'.join(lst)+'\n')
    f.close()

def RemoveFiles(filenames):
    for f in filenames:
        try:
            os.remove(f)
        except OSError as e:
            if e.errno != errno.ENOENT: raise
            continue

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

    def find_data(self, tile):
        """ Find all data for tile, add products dictionary to info dict from inspect function """
        # find original datafile
        filenames = self.find_raw(tile)
        if len(filenames) == 0: return {}
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
        products = {}
        for f in files:
            fname,ext = os.path.splitext(os.path.split(f)[1])
            products[ fname[len(info['basename'])+1:]  ] = f
        # extend info['products'] instead of replace ?
        info['products'] = products
        return info

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

    @classmethod
    def fetch(cls):
        """ Download data and add to archive """
        raise Exception("Fetch not implemented for %s" % cls.name)

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
                VerboseOut(traceback.format_exc(), 4)
                qname = os.path.join(qdir,f)
                if not os.path.exists(qname):
                    os.link(os.path.abspath(f),os.path.join(qdir,f))
                VerboseOut('%s -> quarantine (file error)' % f,2)
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
                os.chmod(os.path.join(dirname,f),0664)
            except:
                pass
        return {'headerfile': hdrfile, 'datafiles': datafiles}

    def __init__(self, site=None, tiles=None, date=None, products=None, sensors=None, **kwargs):
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

        try:
            if args.command == 'archive':
                cls.archive(link=args.link)
                exit(1)
            inv = cls.inventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, 
                products=args.products, pcov=args.pcov, ptile=args.ptile)
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


class DataInventory(object):
    """ Manager class for data inventories """
    # redo color, combine into ordered dictionary
    _colororder = ['purple', 'bright red', 'bright green', 'bright blue','bright purple']
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

    def _colorize(self,txt,color):
        return "\033["+self._colorcodes[color]+'m' + txt + "\033[0m"

    @property
    def dates(self):
        """ Get sorted list of dates """
        return [k for k in sorted(self.data)]

    @property
    def numdates(self):
        """ Get number of dates """
        return len(self.data)

    def __getitem__(self,date):
        return self.data[date]

    def __init__(self, dataclass, site=None, tiles=None, dates=None, days=None, products=None, **kwargs):
        self.dataclass = dataclass
        self.site = site
        self.tiles = tiles
        self.temporal_extent(dates, days)
        self.data = {}
        if products is not None:
            if len(products) == 0: products = dataclass._products.keys()
        self.products = products
        self.AddData(dataclass, **kwargs)

    def AddData(self, dataclass, **kwargs):
        """ Add additional data to this inventory (usually from different sensors """
        if self.tiles is None and self.site is None:
            self.tiles = dataclass.find_tiles()
            #raise Exception('No shapefile or tiles provided for inventory')
        if self.tiles is None and self.site is not None:
            self.tiles = dataclass.vector2tiles(gippy.GeoVector(self.site),**kwargs)
        # get all potential matching dates for tiles
        dates = []
        for t in self.tiles:
            try:
                for date in dataclass.find_dates(t):
                    day = int(date.strftime('%j'))
                    if (self.start_date <= date <= self.end_date) and (self.start_day <= day <= self.end_day):
                        if date not in dates: dates.append(date)
            except: 
                VerboseOut(traceback.format_exc(),4)

        self.numfiles = 0
        for date in sorted(dates):
            try:
                dat = dataclass(site=self.site, tiles=self.tiles, date=date, products=self.products, **kwargs)
                self.data[date] = [ dat ]
                self.numfiles = self.numfiles + len(dat.tiles)
            except: 
                VerboseOut(traceback.format_exc(),4)

    def temporal_extent(self, dates, days):
        """ Temporal extent (define self.dates and self.days) """
        if dates is None: dates='1984,2050'
        self.start_date,self.end_date = dateparse.range(dates)
        if days:
            days = days.split(',')
        else: days = (1,366)
        self.start_day,self.end_day = ( int(days[0]), int(days[1]) )

    def process(self, *args, **kwargs):
        """ Process data in inventory """
        if self.products is None:
            raise Exception('No products specified for processing')
        start = datetime.now()
        VerboseOut('Requested %s products for %s files' % (len(self.products), self.numfiles))
        for date in self.dates:
            for data in self.data[date]:
                data.process(*args, **kwargs)
        VerboseOut('Completed processing in %s' % (datetime.now()-start))

    def project(self, *args, **kwargs):
        self.process()
        start = datetime.now()
        VerboseOut('Projecting data for %s dates (%s - %s)' % (len(self.dates),self.dates[0],self.dates[-1]))
        # res should default to data?
        for date in self.dates:
            for data in self.data[date]:
                data.project(*args, **kwargs)
        VerboseOut('Completed projecting in %s' % (datetime.now()-start))

    def links(self,hard=False):
        """ Create links to tiles - move linking to core """
        for date in self.data:
            for data in self.data[date]:
                for t in data.tiles:
                    for p in data.tiles[t]['products']:
                        fname = data.tiles[t]['products'][p]
                        print p,fname
                        if hard:
                            f = os.link
                        else: f = os.symlink
                        try:
                            f( fname, os.path.basename(fname) )
                        except:
                            print 'problem'

    # TODO - check if this is needed
    def get_products(self, date):
        """ Get list of products for given date """
        # this doesn't handle different tiles (if prod exists for one tile, it lists it)
        prods = []
        for data in self.data[date]:
            for t in data.tiles:
                for p in data.tiles[t]['products']:
                    prods.append(p)
                #for prod in data.products.keys(): prods.append(prod)
        return sorted(set(prods))

    def printcalendar(self,md=False):
        """ print calendar for raw original datafiles """
        #import calendar
        #cal = calendar.TextCalendar()
        oldyear = ''

        print '%s INVENTORY' % self.dataclass.name

        # print tile coverage
        if self.site is not None:
            print '{:^8}{:>14}{:>14}'.format('Tile','% Coverage','% Tile Used')
            for t in sorted(self.tiles): 
                print "{:>8}{:>11.1f}%{:>11.1f}%".format(t,self.tiles[t][0]*100,self.tiles[t][1]*100)
        # print inventory
        for date in self.dates:
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
                if self.products: sys.stdout.write('\n ')
            colors = {}
            for i,s in enumerate(self.dataclass.sensor_names()): colors[s] = self._colororder[i]

            for dat in self.data[date]:
                sys.stdout.write(self._colorize('{:<6}'.format(daystr), colors[self.dataclass.sensors[dat.sensor]] ))
            if self.products:
                sys.stdout.write('        ')
                prods = [p for p in self.get_products(date) if p in self.products]
                for p in prods:
                    sys.stdout.write(self._colorize('{:<12}'.format(p), colors[self.dataclass.sensors[dat.sensor]] ))
                sys.stdout.write('\n ')
            oldyear = date.year
        sys.stdout.write('\n')
        if self.numfiles != 0:
            self.legend()
            VerboseOut("%s files on %s dates" % (self.numfiles, self.numdates))
        else:
            VerboseOut('No matching files')

    def legend(self):
        sensors = sorted(self.dataclass.sensors.values())
        for i,s in enumerate(sensors):
            print self._colorize(s, self._colororder[i])
            #print self._colorize(self.dataclass.sensors[s], self._colororder[s])

    def get_timeseries(self,product=''):
        """ Read all files as time series """
        # assumes only one sensor row for each date
        img = self.data[self.dates[0]][0].open(product=product)
        for i in range(1,len(self.dates)):
            img.AddBand(self.data[self.dates[i]][0].open(product=product)[0])
        return img