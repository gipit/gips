#!/usr/bin/env python

import os, sys
import ogr 
import datetime
import glob
from shapely.wkb import loads
from shapely.geometry import shape
import gippy
import agspy.utils.dateparse as dateparse
from pdb import set_trace

#from gippy.data.landsat import LandsatData

class DataInventory(object):

    _colorcodes = {
        'black':    '0;30',     'bright gray':  '0;37',
        'blue':     '0;34',     'white':        '1;37',
        'green':    '0;32',     'bright blue':  '1;34',
        'cyan':     '0;36',     'bright green': '1;32',
        'red':      '0;31',     'bright cyan':  '1;36',
        'purple':   '0;35',     'bright red':   '1;31',
        'yellow':   '0;33',     'bright purple':'1;35',
        'dark gray':'1;30',     'bright yellow':'1;33',
        'normal':   '0'
    }

    def __init__(self, site=None, tiles=None, dates=None, days=None, products=None):
        self.spatial_extent(site, tiles)
        self.temporal_extent(dates, days)

    def spatial_extent(self, site, tiles):
        """ Spatial extent (define self.tiles) """
        #set_trace()
        if tiles is None: tiles = os.listdir(self.path())
        self.site = site
        if self.site is not None: 
            self.site_to_tiles(gippy.GeoVector(self.site))
        else: self.tiles = dict((t,1) for t in tiles)

    def temporal_extent(self, dates, days):
        """ Temporal extent (define self.dates and self.days) """
        if dates is None: dates='1984,2050'
        self.start_date,self.end_date = dateparse.range(dates)
        if days: 
            days = days.split(',')
        else: days = (1,366)
        self.start_day,self.end_day = ( int(days[0]), int(days[1]) )

    def __getitem__(self,date):
        return self.data[date]

    def path(self,tile='',date=''):
        """ Path to date or tile directory """
        if tile == '':
            return self.rootdir
        elif date == '':
            return os.path.join(self.rootdir, tile)
        else:
            return os.path.join(self.rootdir, tile, date)

    def _colorize(self,txt,color): 
        return "\033["+self._colorcodes[color]+'m' + txt + "\033[0m"

    @property
    def dates(self):
        """ Get list of all dates """
        return [k for k in sorted(self.data)]

    @property
    def numdates(self): 
        """ Get number of dates """
        return len(self.data)

    @property
    def tile_names(self):
        """ Get list of all tiles """
        return [k for k,v in self.tiles.items()]

    @property
    def sensor_names(self):
        """ Get list of all sensors """
        return [self._colorize(k, self._colors[k]) for k in sorted(self._colors)]

    def sensor_color(self,sensor):
        return self._colors[self.sensors[sensor]]

    #def get_tiles(self, date):
    #    """ Get list of tile dictionaries for given date """    
    #    fnames = []
    #    for f in self.tiles[date]:
    #        fnames.append(f['filename'])
    #    return fnames

    def get_products(self, date):
        """ Get list of products for given date """
        # this doesn't handle different tiles (if prod exists for one tile, it lists it)
        prods = []
        for sensor in self.data[date]:
            for datafile in self.data[date][sensor]:
                for prod in datafile.products:
                    prods.append(prod)
        return sorted(prods)

    def site_to_tiles(self,vector):
        """ Identify sensor tiles that fall within vector """
        import osr
        geom = vector.union()
        ogrgeom = ogr.CreateGeometryFromWkb(geom.wkb)
        tvector = self.get_tile_vector()
        tlayer = tvector.layer
        trans = osr.CoordinateTransformation(vector.layer.GetSpatialRef(), tlayer.GetSpatialRef())
        ogrgeom.Transform(trans)
        geom = loads(ogrgeom.ExportToWkb())
        tlayer.SetSpatialFilter(ogrgeom)
        tiles = {}
        tlayer.ResetReading()
        feat = tlayer.GetNextFeature()
        fldindex = feat.GetFieldIndex(self._tile_attribute)
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkb())
            area = geom.intersection(tgeom).area
            if area != 0: 
                tile = str(feat.GetField(fldindex))
                if len(tile) == 5: tile = '0' + tile
                tiles[tile] = area/geom.area
            feat = tlayer.GetNextFeature()  
        self.tiles = tiles
        return tiles

    def createlinks(self,hard=False):
        """ Create product links """
        for date in self.data:
            for sensor in self.data[date]:
                for f in self.data[date][sensor]:
                    for p in f['products']:
                        link( f['products'][p] )

    #def printprodcal(self,md=False):
    #    """ print calendar for each tile displaying products """
    #    self.printcalendar(md,True)

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
            for s in self.data[date]:
                sys.stdout.write(self._colorize('{:<6}'.format(daystr), self.sensor_color(s) ))
            if products:
                sys.stdout.write('        ')
                prods = self.get_products(date)
                for p in prods:
                    sys.stdout.write(self._colorize('{:<9}'.format(p), self.sensor_color(s) ))
                sys.stdout.write('\n ')
            oldyear = date.year
        sys.stdout.write('\n')
        for s in self.sensor_names: print s
        print self
        if self.site is not None:
            print 'Tile Coverage:'
            for t in sorted(self.tiles): print '%s: %2.0f%%' % (t,self.tiles[t]*100)
        
    def __str__(self):
        if self.numfiles != 0:
            s = "Data Inventory: %s files on %s dates" % (self.numfiles,self.numdates)
        else: s = 'Data Inventory: No matching files'
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
