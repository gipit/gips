#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  matt.a.hanson@gmail.com
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

import os
import sys
import argparse
import traceback

from gips.utils import settings, data_sources
from gips.data.core import repository_class
import gippy


class GIPSParser(argparse.ArgumentParser):
    """ Extends argparser parser to print help on error """

    def __init__(self, datasources=True, **kwargs):
        super(GIPSParser, self).__init__(**kwargs)
        self.datasources = datasources
        self.formatter_class = argparse.ArgumentDefaultsHelpFormatter
        self.parent_parsers = []

    def parse_args(self, **kwargs):
        if self.datasources:
            self.add_data_sources()
        args = super(GIPSParser, self).parse_args(**kwargs)
        set_gippy_options(args)
        return args

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

    def add_parser(self, parser):
        self.parent_parsers.append(parser)

    def add_default_parser(self):
        """ This adds a parser with default options """
        if self.datasources:
            parser = GIPSParser(add_help=False)
        else:
            parser = self
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
        self.parent_parsers.append(parser)
        return parser

    def add_inventory_parser(self, site_required=False):
        """ This adds a parser with inventory options """
        if self.datasources:
            parser = GIPSParser(add_help=False)
        else:
            parser = self
        group = parser.add_argument_group('inventory options')
        h = 'Vector layer (file or db) for region of interest'
        group.add_argument('-s', '--site', help=h, default=None, required=site_required)
        h = 'Attribute to use as lookup in in vector file (defaults to index)'
        group.add_argument('-k', '--key', help=h, default="")
        group.add_argument('-w', '--where', help="attribute=value pairs to limit features", nargs='*')
        group.add_argument('-t', '--tiles', nargs='*', help='Tile designations', default=None)
        group.add_argument('-d', '--dates', help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
        group.add_argument('--days', help='Include data within these days of year (doy1,doy2)', default=None)
        group.add_argument('--sensors', help='Sensors to include', nargs='*', default=None)
        group.add_argument('--%cov', dest='pcov', help='Threshold of %% tile coverage over site', default=0, type=int)
        group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)
        group.add_argument('--fetch', help='Fetch any missing data (if supported)', default=False, action='store_true')
        group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
        group.add_argument('-p', '--products', help='Requested Products', nargs='*')
        self.parent_parsers.append(parser)
        return parser

    def add_process_parser(self):
        """ This adds a parser with processing options """
        if self.datasources:
            parser = GIPSParser(add_help=False)
        else:
            parser = self
        group = parser.add_argument_group('processing options')
        group.add_argument('--overwrite', help='Overwrite existing output file(s)', default=False, action='store_true')
        group.add_argument('--chunksize', help='Chunk size in MB', default=128.0, type=float)
        group.add_argument('--numprocs', help='Desired number of processors (if allowed)', default=2, type=int)
        group.add_argument('--format', help='Format for output file', default="GTiff")
        self.parent_parsers.append(parser)
        return parser

    def add_project_parser(self):
        """ This adds a parser with project options """
        if self.datasources:
            parser = GIPSParser(add_help=False)
        else:
            parser = self
        group = parser.add_argument_group('project directory options')
        h = 'Directory to store project(s) (default to current directory)'
        group.add_argument('--outdir', help=h, default='')
        group.add_argument('--suffix', help='Suffix to add to auto generated output directory', default='')
        h = 'Do not create a top-level directory to hold project directories'
        group.add_argument('--notld', help=h, default=False, action='store_true')
        h = 'Create project directories in tree form'
        group.add_argument('--tree', help=h, default=False, action='store_true')
        self.parent_parsers.append(parser)
        return parser

    def add_warp_parser(self):
        """ This adds a parser with warping options """
        if self.datasources:
            parser = GIPSParser(add_help=False)
        else:
            parser = self
        group = parser.add_argument_group('warp options')
        h = 'Output resolution in site projected coordinates (no warping done if not provided)'
        group.add_argument('--res', nargs=2, help=h, default=None, type=float)
        h = 'If warping interpolate using: 0-NN, 1-Bilinear, 2-Cubic'
        group.add_argument('--interpolation', help=h, choices=[0, 1, 2], default=0, type=int)
        group.add_argument('--crop', help='Crop down to minimum bounding box', default=False, action='store_true')
        self.parent_parsers.append(parser)
        return parser

    def add_projdir_parser(self):
        """ This adds a parser with options for reading a project output directory """
        if self.datasources:
            parser = GIPSParser(add_help=False)
        else:
            parser = self
        group = parser.add_argument_group('input project options')
        group.add_argument('projdir', help='GIPS Project directory', nargs='*')
        group.add_argument('-p', '--products', help='Products to operate on', nargs='*')
        self.parent_parsers.append(parser)
        return parser

    def add_data_sources(self):
        """ Adds data sources to parser """
        subparser = self.add_subparsers(dest='command')
        for src, desc in data_sources().items():
            subparser.add_parser(src, help=desc, parents=self.parent_parsers)


def set_gippy_options(args):
    """ Set gippy options from parsed command line arguments """
    if 'verbose' in args:
        gippy.Options.SetVerbose(args.verbose)
    if 'format' in args:
        gippy.Options.SetDefaultFormat(args.format)
    if 'chunksize' in args:
        gippy.Options.SetChunkSize(args.chunksize)
    if 'numprocs' in args:
        gippy.Options.SetNumCores(args.numprocs)
