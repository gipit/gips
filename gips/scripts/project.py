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
from gips.parsers import GIPSParser
from gips.data.core import data_class
from gips.utils import Colors, VerboseOut, basename, mkdir


def main():
    title = Colors.BOLD + 'GIPS Data Project Utility v%s' % gipsversion + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser0.add_inventory_parser()
    parser0.add_process_parser()
    parser0.add_project_parser()
    parser0.add_warp_parser()
    parser0.add_data_sources()
    args = parser0.parse_args()

    try:
        print title
        cls = data_class(args.command)
        inv = cls.inventory(**vars(args))

        # create project directory
        suffix = '' if args.suffix is None else '_' + args.suffix
        if args.datadir is None:
            if args.res is None:
                args.res = cls.Asset._defaultresolution
            if args.res[0] == args.res[1]:
                resstr = str(args.res[0])
            else:
                resstr = '%sx%s' % (args.res[0], args.res[1])
            args.datadir = '%s_%s_%s%s' % (basename(args.site).replace('_', '-'), resstr, args.command, suffix)
        mkdir(args.datadir)

        # warp and mosaic
        for date in inv.dates:
            inv[date].process(overwrite=False)
            inv[date].mosaic(datadir=args.datadir, res=args.res, interpolation=args.interpolation,
                             crop=args.crop, overwrite=args.overwrite)

    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Data Project error: %s' % e


if __name__ == "__main__":
    main()
