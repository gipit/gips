#!/usr/bin/env python

import sys
import os
import gippy
from gippy.utils import VerboseOut, parse_dates

from datetime import datetime
import traceback

from pdb import set_trace

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

    def __init__(self, dataclass, site=None, tiles=None, dates=None, days=None, products=None, fetch=False, **kwargs):
        self.dataclass = dataclass

        self.site = site
        # default to all tiles
        if tiles is None and self.site is None:
            tiles = dataclass.find_tiles()
        # if tiles provided, make coverage all 100%
        if tiles is not None:
            self.tiles = {}
            for t in tiles: self.tiles[t] = (1,1)
        elif tiles is None and self.site is not None:
            self.tiles = dataclass.vector2tiles(gippy.GeoVector(self.site),**kwargs)

        self.temporal_extent(dates, days)

        self.data = {}
        if products is not None:
            if len(products) == 0: products = dataclass._products.keys()
        self.products = products

        if fetch and products is not None:
            dataclass.fetch(products,self.tiles,(self.start_date,self.end_date),(self.start_day,self.end_day))
            exit(1)
            # archive
            #dataclass.archive(os.path.join(dataclass._rootdir),'staging')

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
                dat = dataclass(site=self.site, tiles=self.tiles.keys(), date=date, products=self.products, **kwargs)
                self.data[date] = [ dat ]
                self.numfiles = self.numfiles + len(dat.tiles)
            except Exception,e:
                VerboseOut('Inventory error %s' % e,3)
                VerboseOut(traceback.format_exc(),4)

    #def fetch(self):
        """ fetch data from remote location """
        # call fetch with range of dates, then archive 
    #    for t in self.tile_coverage.keys():
    #        set_trace()
    #        self.fetch(t)
        # archive files
        #self.archive(os.path.join(self._rootdir,'staging'))

    def temporal_extent(self, dates, days):
        """ Temporal extent (define self.dates and self.days) """
        if dates is None: dates='1984,2050'
        self.start_date,self.end_date = parse_dates(dates)
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
                        bname = os.path.basename(fname)
                        if hard:
                            f = os.link
                        else: f = os.symlink
                        try:
                            f( fname, bname )
                            VerboseOut('%s: linking' % bname)
                        except: 
                            VerboseOut('%s: Problem creating link' % bname,2)

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

        # print tile coverage
        if self.site is not None:
            print '\nTILE COVERAGE'
            print '{:^8}{:>14}{:>14}'.format('Tile','% Coverage','% Tile Used')
            for t in sorted(self.tiles): 
                print "{:>8}{:>11.1f}%{:>11.1f}%".format(t,self.tiles[t][0]*100,self.tiles[t][1]*100)
        # print inventory
        print '\nINVENTORY'
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
            VerboseOut("\n%s files on %s dates" % (self.numfiles, self.numdates))
            self.legend()
        else:
            VerboseOut('No matching files')

    def legend(self):
        print '\nSENSORS'
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