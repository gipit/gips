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

from gips import __version__ as gipsversion
from gips.parsers import GIPSParser, inventory_parser
from gips.data.core import data_class
from gips.utils import Colors, VerboseOut


def main():
    title = Colors.BOLD + 'GIPS Data Processing Utility v%s' % gipsversion + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser = inventory_parser()
    group = parser.add_argument_group('Processing Options')
    group.add_argument('--overwrite', help='Overwrite exiting output file(s)', default=False, action='store_true')
    group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)
    parser0.add_data_sources(parents=[parser])
    args = parser0.parse_args()

    try:
        print title
        cls = data_class(args.command)
        inv = cls.inventory(**vars(args))
        inv.process(overwrite=args.overwrite)
    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Data processing error: %s' % e


if __name__ == "__main__":
    main()
