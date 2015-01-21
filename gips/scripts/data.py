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

import argparse
from gips import __version__ as gipsversion
from gips.data.core import data_class
from gips.parsers import GIPSParser, add_data_sources


def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    h0 = 'GIPS v%s Data Repositories' % (gipsversion)
    parser0 = GIPSParser(description=h0, formatter_class=dhf)

    add_data_sources(parser0)

    args = parser0.parse_args()

    print h0
    cls = data_class(args.command)
    cls.print_products()
