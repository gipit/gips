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

import sys
import argparse
import traceback

import gips.settings as settings
from gips.data.core import repository_class


class GIPSParser(argparse.ArgumentParser):
    """ Extends argparser parser to print help on error """
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def add_data_sources(parser):
    """ This adds available data sources as subparsers """
    subparser = parser.add_subparsers(dest='command')
    for key in sorted(settings.REPOS.keys()):
        # get description
        try:
            repo = repository_class(settings.REPOS[key]['class'])
            subparser.add_parser(key, help=repo.description)
        except:
            print traceback.format_exc()
    return parser


def add_inventory_parser(parser):
    """ This adds inventory arguments to an argument parser """
    group = parser.add_argument_group('inventory arguments')
    group.add_argument('data', help='Data type to use (e.g. landsat)')
    group.add_argument('-s', '--site', help='Vector file for region of interest', default=None)
    group.add_argument('-t', '--tiles', nargs='*', help='Tile designations', default=None)
    group.add_argument('-d', '--dates', help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    group.add_argument('--days', help='Include data within these days of year (doy1,doy2)', default=None)
    group.add_argument('--sensors', help='Sensors to include', nargs='*', default=None)
    group.add_argument('--%cov', dest='pcov', help='Threshold of %% tile coverage over site', default=0, type=int)
    group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)
    group.add_argument('--fetch', help='Fetch any missing data (if supported)', default=False, action='store_true')
    group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
    group.add_argument('-p', '--products', help='Requested Products (call products command to list)', nargs='*')
