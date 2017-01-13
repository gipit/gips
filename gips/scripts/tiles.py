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
from gips import __version__
from gips.parsers import GIPSParser
from gips.core import SpatialExtent, TemporalExtent
from gips.utils import Colors, VerboseOut, mkdir, open_vector, import_data_class
from gips import utils
from gips.inventory import DataInventory
from gips.inventory import orm


def main():
    title = Colors.BOLD + 'GIPS Tiles (v%s)' % __version__ + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser0.add_inventory_parser()
    parser0.add_process_parser()
    parser0.add_project_parser()
    parser0.add_warp_parser()
    args = parser0.parse_args()

    cls = utils.gips_script_setup(args.command, args.stop_on_error)
    print title

    with utils.error_handler():
        # create output directory if needed
        # tld is "{}_tiles_{}_{}".format(DATATYPE, RESOLUTION, SUFFIX)
        if args.notld:
            tld = args.outdir
        else:
            tld = os.path.join(args.outdir, '%s_tiles' % args.command)
            if args.res is not None:
                tld = tld + '_%sx%s' % (args.res[0], args.res[1])
            if args.suffix != '':
                tld = tld + '_' + args.suffix
        mkdir(tld)

        extents = SpatialExtent.factory(cls, args.site, args.key, args.where, args.tiles, args.pcov, args.ptile)
        for extent in extents:
            inv = DataInventory(cls, extent, TemporalExtent(args.dates, args.days), **vars(args))
            for date in inv.dates:
                for tid in inv[date].tiles:
                    # make sure back-end tiles are processed
                    inv[date].tiles[tid].process(args.products, overwrite=False)
                    # warp the tiles & copy into place in the output dir
                    inv[date].tiles[tid].copy(tld, args.products, inv.spatial.site,
                                              args.res, args.interpolation, args.crop, args.overwrite, args.tree)

    utils.gips_exit() # produce a summary error report then quit with a proper exit status


if __name__ == "__main__":
    main()
