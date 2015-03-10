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

from gips import __version__ as gipsversion
from gips.parsers import GIPSParser
from gips.core import SpatialExtent, TemporalExtent
from gips.data.core import data_class
from gips.utils import Colors, VerboseOut, open_vector
from gips.inventory import DataInventory


def main():
    title = Colors.BOLD + 'GIPS Data Inventory (v%s)' % gipsversion + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser = parser0.add_inventory_parser()
    group = parser.add_argument_group('inventory display')
    group.add_argument('--md', help='Show dates using MM-DD', action='store_true', default=False)
    args = parser0.parse_args()

    try:
        print title
        cls = data_class(args.command)

        extents = SpatialExtent.factory(cls, args.site, args.key, args.where, 
                                        args.tiles, args.pcov, args.ptile)
        for extent in extents:
            inv = DataInventory(cls, extent, TemporalExtent(args.dates, args.days), **vars(args))
            inv.pprint(md=args.md)            
           
    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Data inventory error: %s' % e

"""
def test():
    from gips.utils import data_sources
    import sys
    args = ['', '-s /etc/gips/test/NH.shp', '-v 3']
    for a in args:
        sys.argv.append(a)
    for src in data_sources():
        sys.argv[1] = src
        main()
"""

if __name__ == "__main__":
    main()
