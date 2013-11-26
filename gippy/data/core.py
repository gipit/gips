#!/usr/bin/env python

import os, sys
import argparse
import ogr 
import datetime
import glob
from shapely.wkb import loads
from shapely.geometry import shape
import traceback
import gippy
import agspy.utils.dateparse as dateparse
from pdb import set_trace

def VerboseOut(txt, level=1):
    if gippy.Options.Verbose() >= level: print txt

class Data(object):
    """ Base class for data objects """
    name = ''
    sensors = {}
    rootdir = ''
    _tiles_vector = ''
    _tiles_attribute = ''

    @classmethod
    def find_tiles(cls):
        """ Get list of all available tiles """
        return os.listdir(cls.rootdir)

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates for a tile """
        return [datetime.datetime.strptime(os.path.basename(d),'%Y%j').date() for d in os.listdir(cls.path(tile))]

    @classmethod
    def find_products(cls, tile, date, products):
        """ Find all products for specified tile dictionary. Create tile dictionary containing:
             path - path to files/products
             basename - base/root name of this tile/date
             products - dictionary of product name and filename
             sensor - name of sensor
        """
        return {'path': '', 'basename': '', 'sensor': '', 'products': {} }

    @classmethod
    def filter(cls, tile, filename):
        """ Check if tile passes filter """
        # TODO - remove filename parameter
        return True

    @classmethod
    def path(cls,tile,date=''):
        """ Path to date or tile directory (assuming tiledir/datedir structure """
        if date == '':
            return os.path.join(cls.rootdir, tile)
        else:
            return os.path.join(cls.rootdir, tile, date)

    @classmethod
    def fetch(cls):
        """ Download data and add to archive """
        raise Exception("Fetch not implemented for %s" % cls.name)

    @classmethod
    def archive(cls, path=''):
        """ Move files from directory to archive location """
        raise Exception("Archive not implemented for %s" % cls.name)
        pass

    def process(self):
        """ Make sure all products exist and process if needed """
        pass

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    @classmethod
    def inventory(cls,site=None, tiles=None, dates=None, days=None, products=None, **kwargs):
        return DataInventory(cls, site, tiles, dates, days, products, **kwargs)

    @classmethod
    def sensor_names(cls):
        """ All possible sensor names """
        return sorted(cls.sensors.values())

    @classmethod 
    def get_tiles_vector(cls):
        """ Get GeoVector of sensor grid """
        return gippy.GeoVector("PG:dbname=geodata host=congo port=5432 user=ags", layer=cls._tiles_vector)

    @classmethod
    def vector2tiles(cls, vector, mincoverage=5):
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
        fldindex = feat.GetFieldIndex(cls._tiles_attribute)
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkb())
            area = geom.intersection(tgeom).area
            if area != 0: 
                tile = str(feat.GetField(fldindex))
                # TODO - THIS IS LANDSAT SPECIFIC !
                if len(tile) == 5: tile = '0' + tile
                tiles[tile] = area/geom.area
            feat = tlayer.GetNextFeature()
        remove_tiles = []
        for t in tiles:
            if tiles[t] < mincoverage/100.0: remove_tiles.append(t)
        for t in remove_tiles: tiles.pop(t,None)
        return tiles

    def open(self, product=''):
        """ Open and return final product GeoImage """
        if product != '':
            return gippy.GeoImage(self.products[product])
        elif len(self.products) == 1:
            return gippy.GeoImage(self.products[self.products.keys()[0]])
        else:
            # return filename of a tile from self.tiles ?
            raise Exception('No product provided')

    def project(self, res, datadir='gipdata'):
        """ Create image of final product (reprojected/mosaiced) """
        self.process()
        if not os.path.exists(datadir): os.makedirs(datadir)
        datadir = os.path.abspath(datadir)
        if len(res) == 1: res = [res,res]
        if self.site is None:
            raise Exception("No site file supplied")
        for product in self.products:
            if self.products[product] == '':
                start = datetime.datetime.now()
                filename = os.path.join(datadir, self.date.strftime('%Y%j') + '_%s_%s.tif' % (product,self.sensor))
                if not os.path.exists(filename):
                    filenames = [self.tiles[t]['products'][product] for t in self.tiles]
                    imgout = gippy.CookieCutter(filenames, filename, self.site, res[0], res[1])
                    print 'Projected and cropped %s files -> %s in %s' % (len(filenames),imgout.Basename(),datetime.datetime.now() - start)
                self.products[product] = filename

    def __init__(self, site=None, tiles=None, date=None, products=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a tile dictionary (see next)
        self.tiles[tile] - dictionary of tile data
        self.products - dictionary of product name and final product filename
        """
        self.site = site
        # Calculate spatial extent
        if tiles is not None:
            self.tile_coverage = dict((t,1) for t in tiles)
            self.tiles = tiles
        elif site is not None:
            self.tile_coverage = self.vector2tiles(gippy.GeoVector(site))
            self.tiles = self.tile_coverage.keys()
        else:
            self.tile_coverage = dict((t,1) for t in self.find_tiles())
        self.date = date
        # Create tile and product dictionaries for use by child class
        self.tiles = {}
        for t in self.tile_coverage.keys(): self.tiles[t] = {}
        if products is None: products = self._products.keys()
        if len(products) == 0: products = self._products.keys()
        self.products = {}
        for p in products: self.products[p] = ''
        # For each tile locate files
        empty_tiles = []
        for t in self.tiles:
            try:
                self.tiles[t] = self.find_products(t,date,products)
                self.sensor = self.tiles[t]['sensor']
                #print self.tiles[t]
            except Exception,e:
                #print 'Error: %s' % (e)
                #VerboseOut(traceback.format_exc(), 3)
                empty_tiles.append(t)
                continue

            # Custom filter based on dataclass
            #good = self.filter(t,filename, **kwargs)
            #if good == False:
            #    empty_tiles.append(t)

        for t in empty_tiles: self.tiles.pop(t,None)
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    def __str__(self):
        return self.sensor + ': ' + str(self.date)


class DataInventory(object):
    """ Base class for data inventories """
    # redo color, combine into ordered dictionary
    _colororder = ['bright yellow', 'bright red', 'bright green', 'bright blue']
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

    def __init__(self, dataclass, site=None, tiles=None, dates=None, days=None, products=None, **kwargs):
        self.dataclass = dataclass
        self.site = site
        self.tiles = tiles
        self.temporal_extent(dates, days)
        self.AddData(dataclass, products=products, **kwargs)

    def __getitem__(self,date):
        return self.data[date]

    #def filenames(self,product):
    #    """ Return dictionary (date keys) of filenames for given sensor and product if supplied """
    #    filenames = {}
    #    for date in self.dates:
    #        filenames[date] = self.data[date][0].filename(product)
    #    return filenames

    @property
    def dates(self):
        """ Get sorted list of dates """
        return [k for k in sorted(self.data)]

    @property
    def numdates(self): 
        """ Get number of dates """
        return len(self.data)

    def get_timeseries(self,product=''):
        """ Read all files as time series """
        # assumes only one sensor row for each date
        img = self.data[self.dates[0]][0].open(product=product)
        for i in range(1,len(self.dates)):
            img.AddBand(self.data[self.dates[i]][0].open(product=product)[0])
        return img

    #@property
    #def sensor_names(self):
    #    """ Get list of all sensors """
    #    return [self._colorize(k, self._colors[k]) for k in sorted(self._colors)]

    def _colorize(self,txt,color): 
        return "\033["+self._colorcodes[color]+'m' + txt + "\033[0m"

    def AddData(self, dataclass, products=None, **kwargs):
        """ Add additional data to this inventory (usually from different sensors """
        if self.tiles is None and self.site is None:
            raise Exception('No shapefile or tiles provided for inventory')
        if self.tiles is None and self.site is not None:
            self.tiles = dataclass.vector2tiles(gippy.GeoVector(self.site))
        # get all potential matching dates for tiles
        self.products = products
        dates = []
        for t in self.tiles:
            try:
                for date in dataclass.find_dates(t):
                    day = int(date.strftime('%j'))
                    if (self.start_date <= date <= self.end_date) and (self.start_day <= day <= self.end_day):
                        if date not in dates: dates.append(date)
            except: pass
        self.numfiles = 0
        self.data = {}
        for date in sorted(dates):
            try:
                dat = dataclass(site=self.site, tiles=self.tiles, date=date, products=products, **kwargs)
                self.data[date] = [ dat ]
                self.numfiles = self.numfiles + len(dat.tiles)
            except: pass
        
    def temporal_extent(self, dates, days):
        """ Temporal extent (define self.dates and self.days) """
        if dates is None: dates='1984,2050'
        self.start_date,self.end_date = dateparse.range(dates)
        if days: 
            days = days.split(',')
        else: days = (1,366)
        self.start_day,self.end_day = ( int(days[0]), int(days[1]) )

    def process(self, overwrite=False, suffix='', overviews=False):
        """ Process all data in inventory """
        VerboseOut('Requested %s products for %s files' % (len(self.products), self.numfiles))
        # TODO only process if don't exist
        for date in self.dates:
            for data in self.data[date]:
                data.process(overwrite, suffix, overviews)
                # TODO - add completed product(s) to inventory         
        VerboseOut('Completed processing')

    def project(self, res=None, datadir='gipdata'):
        VerboseOut('Preparing data for %s dates (%s - %s)' % (len(self.dates),self.dates[0],self.dates[-1]))
        # res should default to data?
        for date in self.dates:
            for data in self.data[date]:
                data.project(res, datadir=datadir)

    # TODO - check if this is needed
    def get_products(self, date):
        """ Get list of products for given date """
        # this doesn't handle different tiles (if prod exists for one tile, it lists it)
        prods = []
        for data in self.data[date]:
            for prod in data.products.keys():
                prods.append(prod)
        return sorted(prods)

    def createlinks(self,hard=False):
        """ Create product links """
        for date in self.data:
            for data in self.data[date]:
                for t in data.tiles:
                    for p in data.tiles[t]['products']:
                        link( data.tiles[t]['products'][p], hard )

    def printcalendar(self,md=False, products=False):
        """ print calendar for raw original datafiles """
        #import calendar
        #cal = calendar.TextCalendar()
        oldyear = ''
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
                if products: sys.stdout.write('\n ')
            colors = {}
            for i,s in enumerate(self.dataclass.sensor_names()): colors[s] = self._colororder[i]

            for dat in self.data[date]:
                sys.stdout.write(self._colorize('{:<6}'.format(daystr), colors[self.dataclass.sensors[dat.sensor]] ))
            if products:
                sys.stdout.write('        ')
                prods = self.get_products(date)
                for p in prods:
                    sys.stdout.write(self._colorize('{:<12}'.format(p), colors[self.dataclass.sensors[dat.sensor]] ))
                sys.stdout.write('\n ')
            oldyear = date.year
        sys.stdout.write('\n')
        self.legend()
        print self
        if self.site is not None:
            print 'Tile Coverage:'
            for t in sorted(self.tiles): print '%s: %2.0f%%' % (t,self.tiles[t]*100)

    def legend(self):
        sensors = sorted(self.dataclass.sensors.values())
        for i,s in enumerate(sensors):
            print self._colorize(s, self._colororder[i])
            #print self._colorize(self.dataclass.sensors[s], self._colororder[s])
        
    def __str__(self):
        if self.numfiles != 0:
            s = "Data Inventory: %s files on %s dates" % (self.numfiles,self.numdates)
        else: 
            s = 'Data Inventory: No matching files'
        return s

def link(f,hard=False):
    """ Create link to file in current directory """
    faux = f + '.aux.xml'
    if hard:
        try:
            os.link(f,os.path.basename(f))
            os.link(faux,os.path.basename(faux))
        except:
            pass
    else: 
        try:
            os.symlink(f,os.path.basename(f))
            if os.path.isfile(faux):
                os.symlink(faux,os.path.basename(faux))
        except:
            pass

def main(dataclass):
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='%s Data Utility' % dataclass.name, formatter_class=argparse.RawTextHelpFormatter)
    subparser = parser0.add_subparsers(dest='command')

    invparser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = invparser.add_argument_group('inventory arguments')
    group.add_argument('-s','--site',help='Vector file for region of interest', default=None)
    group.add_argument('-t','--tiles', nargs='*', help='Tile designations', default=None)
    group.add_argument('-d','--dates',help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    group.add_argument('--days',help='Include data within these days of year (doy1,doy2)',default=None)
    group.add_argument('-p','--products', nargs='*', help='Process/filter these products') #default=False)
    group.add_argument('-v','--verbose',help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)

    # Help
    parser = subparser.add_parser('help',help='Print extended help', parents=[invparser], formatter_class=dhf)

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
    group.add_argument('--res',nargs=2,help='Resolution of output rasters', default=[30,30], type=float)
    group.add_argument('--datadir', help='Directory to save project files', default='gipdata')

    # Links
    parser = subparser.add_parser('link',help='Link to Products', parents=[invparser], formatter_class=dhf)
    parser.add_argument('--hard',help='Create hard links instead of symbolic', default=False,action='store_true')

    # Misc
    parser_archive = subparser.add_parser('archive',help='Move files from current directory to data archive')

    # Pull in dataclass options here
    #dataparser = subparser.add_parser('data',help='', parents=[invparser],formatter_class=dhf)

    args = parser0.parse_args()
    if args.command == 'help':
        parser0.print_help()
        print '\navailable products:'
        for key,val in dataclass._products.items(): 
            print '    {:<20}{:<100}'.format(key, val['description'])
        exit(1)

    if args.command == 'archive':
        dataclass.archive()
        exit(1)

    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetChunkSize(128.0)   # replace with option

    try:
        inv = dataclass.inventory(site=args.site, dates=args.dates, days=args.days, tiles=args.tiles, products=args.products)
    except Exception,e:
        print 'Error getting inventory: %s' % (e)
        VerboseOut(traceback.format_exc(), 3)
        exit(1)

    if args.command == 'inventory':
        if args.products is None:
            inv.printcalendar(args.md)
        else: inv.printcalendar(args.md,True)
        
    elif args.command == 'link':
        inv.createlinks(args.hard)

    elif args.command == 'process':
        try:
            #merrafname = fetchmerra(meta['datetime'])
            inv.process(overwrite=args.overwrite,suffix=args.suffix) #, nooverviews=args.nooverviews)
        except Exception,e:
            print 'Error processing: %s' % e
            VerboseOut(traceback.format_exc(), 3)

    elif args.command == 'project': inv.project(args.res, datadir=args.datadir)

    else:
        print 'Command %s not recognized' % cmd
