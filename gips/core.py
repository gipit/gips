#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  mhanson@ags.io
#
#    Copyright (C) 2014 Applied Geosolutions
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

import datetime
import calendar

from gips.utils import Colors
from gips.GeoVector import GeoVector


class DataNotFoundException(Exception):
    """ Error thrown when data is not found """
    pass


class Products(object):
    """ Collection of product names and options """
    # TODO - move Products to dataclass specific class in data.core?
    def __init__(self, dataclass, products=None):
        """ Create product object """
        products = products if products is not None else dataclass._products.keys()
        self.products = {p: p.split('-') for p in products}
        self.standard = {}
        self.composite = {}
        for p, val in self.products.items():
            if val[0] not in dataclass._products:
                raise Exception('Invalid product %s' % val[0])
            if dataclass._products[val[0]].get('composite', False):
                self.composite[p] = val
            else:
                self.standard[p] = val

    @property
    def requested(self):
        """ Return list of requested products """
        return [val[0] for val in self.products.values()]

    def __str__(self):
        return ' '.join(self.products.keys())


class SpatialExtent(object):
    """ Description of spatial extent """

    def __init__(self, dataclass, site=None, tiles=None, pcov=0.0, ptile=0.0):
        """ Create spatial extent object """
        self.repo = dataclass.Asset.Repository

        # Spatial extent
        if tiles is None and site is None:
            tiles = self.repo.find_tiles()
        if tiles is not None:
            # if tiles only provided, coverage of each is 100%
            tiles = {t: (1, 1) for t in tiles}
        elif tiles is None and site is not None:
            tiles = self.repo.vector2tiles(GeoVector(site), pcov, ptile)

        self.site = site
        self.coverage = tiles

    @property
    def tiles(self):
        return self.coverage.keys()

    @property
    def available_dates(self):
        """ Get list of all dates for these tiles """
        dates = []
        for t in self.tiles:
            dates.extend(self.repo.find_dates(t))
        return dates

    def print_tile_coverage(self):
        """ Print tile coverage info """
        if self.site is not None:
            print Colors.BOLD + '\nTile Coverage'
            print Colors.UNDER + '{:^8}{:>14}{:>14}'.format('Tile', '% Coverage', '% Tile Used') + Colors.OFF
            for t in sorted(self.coverage):
                print "{:>8}{:>11.1f}%{:>11.1f}%".format(t, self.coverage[t][0] * 100, self.coverage[t][1] * 100)

    def __str__(self):
        return '%s: %s' % (self.site, ' '.join(self.tiles))

    def __len__(self):
        return len(self.coverage)


class TemporalExtent(object):
    """ Description of temporal extent """

    def __init__(self, dates=None, days=None):
        """ Create temporal extent object from string input """
        if dates is None:
            dates = '1950,2050'
        if days is None:
            days = (1, 366)
        else:
            days = days.split(',')
            days = (int(days[0]), int(days[1]))
        try:
            (d1, d2) = dates.replace(',', ' ').split()
            dates = (self._parse_date(d1), self._parse_date(d2, True))
        except:
            dates = (self._parse_date(dates), self._parse_date(dates, True))
        self.datebounds = dates
        self.daybounds = days
        self.datearray = []
        #self.dates = [d for t in tiles for d in repo.find_dates(t) if datecheck(d) and daycheck(d)]

    def prune_dates(self, dates):
        """ Prune down given list of dates to those that meet temporal extent """
        dates = list(set(dates))
        datecheck = lambda d: self.datebounds[0] <= d <= self.datebounds[1]
        daycheck = lambda d: self.daybounds[0] <= int(d.strftime('%j')) <= self.daybounds[1]
        return sorted([d for d in dates if datecheck(d) and daycheck(d)])

    @classmethod
    def _parse_date(cls, dstring, last=False):
        """ Parses string of YYYY or YYYY-MM or YYYY-MM-DD or YYYY-DOY and returns date object """
        d = dstring.split('-')
        if len(d) == 2 and len(d[1]) == 3:
            dttmp = datetime.datetime(int(d[0]), 1, 1) + datetime.timedelta(days=int(d[1]) - 1)
            d[1] = dttmp.month
            d.append(dttmp.day)
        if (not last):
            if (len(d) == 1):
                d.append('1')
            if (len(d) == 2):
                d.append('1')
        else:
            if (len(d) == 1):
                d.append('12')
            if (len(d) == 2):
                d.append(calendar.monthrange(int(d[0]), int(d[1]))[1])
        return datetime.date(int(d[0]), int(d[1]), int(d[2]))

    def __str__(self):
        return '%s - %s (days %s-%s)' % (self.datebounds[0], self.datebounds[1], self.daybounds[0], self.daybounds[1])
