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
import errno
import gippy
import datetime
import calendar


def VerboseOut(obj, level=1):
    if gippy.Options.Verbose() >= level:
        #pprint.PrettyPrinter().pprint(obj)
        if not isinstance(obj, (list, tuple)):
            obj = [obj]
        for o in obj:
            print o


def File2List(filename):
    f = open(filename)
    txt = f.readlines()
    txt2 = []
    for t in txt:
        txt2.append(t.rstrip('\n'))
    return txt2


def List2File(lst, filename):
    f = open(filename, 'w')
    f.write('\n'.join(lst)+'\n')
    f.close()


def RemoveFiles(filenames, extensions=['']):
    for f in filenames:
        for ext in ([''] + extensions):
            try:
                os.remove(f+ext)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise
                continue


def atmospheric_model(doy, lat):
    """ Determine atmospheric model
    1 - Tropical
    2 - Mid-Latitude Summer
    3 - Mid-Latitude Winter
    4 - Sub-Arctic Summer
    5 - Sub-Arctic Winter
    6 - US Standard Atmosphere
    """
    # Determine season
    if doy < 121 or doy > 274:
        if lat < 0:
            summer = True
        else:
            summer = False
    else:
        if lat < 0:
            summer = False
        else:
            summer = True
    # Determine model
    if abs(lat) <= 15:
        model = 1
    elif abs(lat) >= 60:
        if summer:
            model = 4
        else:
            model = 5
    else:
        if summer:
            model = 2
        else:
            model = 3
    return model


def _parse_date(dstring, last=False):
    """ Parses string of YYYY or YYYY-MM or YYYY-MM-DD or YYYY-DOY and returns date object """
    d = dstring.split('-')
    if len(d) == 2 and len(d[1]) == 3:
        dttmp = datetime.datetime(int(d[0]), 1, 1) + datetime.timedelta(days=int(d[1])-1)
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


def parse_dates(dstring):
    """ Parses string of 1 or 2 dates separated by a comma.  Valid formats: YYYY, YYYY-MM, YYYY-MM-DD, YYYY-DOY """
    try:
        (d1, d2) = dstring.replace(',', ' ').split()
        return (_parse_date(d1), _parse_date(d2, True))
    except:
        return (_parse_date(dstring), _parse_date(dstring, True))
