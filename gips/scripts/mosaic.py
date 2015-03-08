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
from gips import __version__ as gipsversion
from gips.parsers import GIPSParser, parse_sites
from gips.data.core import data_class
from gips.utils import Colors, VerboseOut, mkdir, basename, open_vector


def main():
    title = Colors.BOLD + 'GIPS Mosaic Utility v%s' % gipsversion + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser0.add_inventory_parser()
    parser0.add_process_parser()
    parser0.add_project_parser()
    args = parser0.parse_args()

    try:
        print title
        cls = data_class(args.command)

        # create project directory SITENAME_DATATYPE
        suffix = '' if args.suffix is None else '_' + args.suffix

        for feature in open_vector(args.site, args.key, args.where):
            inv = cls.inventory(feature, **vars(args))
            datadir = os.path.join(args.outdir, '%s_%s%s' % (inv.spatial.sitename, args.command, suffix))
            mkdir(datadir)
            for date in inv.dates:
                # make sure back-end tiles are processed
                inv[date].process(overwrite=False)
                # mosaic the tiles
                inv[date].mosaic(datadir=datadir, overwrite=args.overwrite)

    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Mosaic error: %s' % e


if __name__ == "__main__":
    main()
