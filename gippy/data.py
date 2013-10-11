#!/usr/bin/env python

import os, sys
import ogr 
from shapely.wkb import loads

class DataInventory(object):
    _rootdir = ''
    _origdir = ''
    _proddir = ''

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

    def __getitem__(self,date):
        return self.tiles[date]

    def origpath(self):
        return os.path.join(self._rootdir,self._origdir)

    def prodpath(self):
        return os.path.join(self._rootdir,self._proddir)

    def tilepath(self,tile): 
        return os.path.join(self.origpath(),tile)

    def _colorize(self,txt,color): 
        return "\033["+self._colorcodes[color]+'m' + txt + "\033[0m"

    @property
    def tile_names(self):
        """ Get list of all tiles """
        return [k for k,v in self.tile_coverage.items()]

    def dates(self):
        """ Get list of all dates """
        return [k for k in sorted(self.tiles)]

    @property
    def numdates(self): 
        """ Get number of dates """
        return len(self.tiles)

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
        for datafile in self.tiles[date]:
            for prod in datafile['products']:
                prod = os.path.splitext(prod)[0]
                prods.append( prod[len(datafile['basename'])+1:] )
        return (prods)

    def site_to_tiles(self,vector):
        """ Identify sensor tiles that fall within vector """
        geom = vector.union()
        ogrgeom = ogr.CreateGeometryFromWkb(geom.wkb)
        tvector = self.get_tile_vector()
        tlayer = tvector.layer
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
        self.tile_coverage = tiles
        return tiles

    def createlinks(self,hard=False):
        """ Create product links """
        for date in self.tiles:
            for f in self.tiles[date]:
                for p in f['products']:
                    link(os.path.join(f['productpath'], p), hard)

    def printprodcal(self,md=False):
        """ print calendar for each tile displaying products """
        self.printcalendar(md,True)

    def printcalendar(self,md=False, products=False):
        """ print calendar for raw original datafiles """
        #import calendar
        #cal = calendar.TextCalendar()
        oldyear = ''
        for date in self.dates():        
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
            sys.stdout.write(self._colorize('{:<6}'.format(daystr),self.tiles[date][0]['color']))
            if products:
                sys.stdout.write('        ')
                prods = self.get_products(date)
                for p in prods:
                    sys.stdout.write(self._colorize('{:<9}'.format(p),self.tiles[date][0]['color']))
                sys.stdout.write('\n ')
            oldyear = date.year
        sys.stdout.write('\n')
        for s in sorted(self._colors): print self._colorize(s,self._colors[s])
        print self
        
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