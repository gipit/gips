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
from gips.parsers import GIPSParser, add_data_sources, add_inventory_parser, set_gippy_options
from gips.data.core import data_class
from gips.utils import Colors


def main():
    #dhf = argparse.ArgumentDefaultsHelpFormatter
    h = Colors.BOLD + 'GIPS Data Inventory Utility v%s' % gipsversion + Colors.OFF
    parser0 = GIPSParser(description=h)

    add_data_sources(parser0)

    add_inventory_parser(parser0)

    # inventory display options
    group = parser0.add_argument_group('inventory display')
    group.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)
    group.add_argument('--compact', help='Print only dates (no coverage)', default=False, action='store_true')

    args = parser0.parse_args()

    set_gippy_options(args)

    try:
        print h
        cls = data_class(args.command)
        inv = cls.inventory(**vars(args))
        #inv = DataInventory(eval(args.data), site=args.site, tiles=args.tiles, dates=args.dates,
        #                    pcov=args.pcov, ptile=args.ptile, fetch=args.fetch, sensors=args.sensors, **kwargs)
        inv.pprint(md=args.md, compact=args.compact)
    except Exception, e:
        import traceback
        print traceback.format_exc()
        print e
        #VerboseOut(traceback.format_exc(), 4)

if __name__ == "__main__":
    main()
