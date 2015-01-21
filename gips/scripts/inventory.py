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
from gips.parsers import parser_add_inventory
from gips.inventory import DataInventory

def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    h = 'GIPS Data Inventory Utility v%s' % __version__
    parser0 = argparse.ArgumentParser(description=h, formatter_class=dhf)

    parser_add_inventory(parser0)

    args = parser0.parse_args()

    # Set options
    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetDefaultFormat(args.format)
    gippy.Options.SetChunkSize(args.chunksize) 
   
    # turn requested data into classname
    eval('from gips.data.%s import %s as dataclass' % args.data)
 
    try:
        print Colors.BOLD + 'GIPS Inventory v%s for %s' % (__version__, args.data) + Colors.OFF
        #inv = DataInventory(eval(args.data), site=args.site, tiles=args.tiles, dates=args.dates,
        #                    pcov=args.pcov, ptile=args.ptile, fetch=args.fetch, sensors=args.sensors, **kwargs) 
        inv = DataInventory(eval(args.data), **args)
        inv.pprint(md=args.md, compact=args.compact)
    except Exception, e:
        VerboseOut(traceback.format_exc(), 4)
        print 'Error with %s: %s' % (args.data, e)

if __name__ == "__main__":
    main()
