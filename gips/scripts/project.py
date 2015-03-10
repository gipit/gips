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
from gips.data.core import data_class
from gips.utils import Colors, VerboseOut, open_vector
from gips.inventory import DataInventory


def main():
    title = Colors.BOLD + 'GIPS Data Project (v%s)' % __version__ + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser0.add_inventory_parser(site_required=True)
    parser0.add_process_parser()
    parser0.add_project_parser()
    parser0.add_warp_parser()
    args = parser0.parse_args()

    try:
        print title
        cls = data_class(args.command)

        extents = SpatialExtent.factory(cls, args.site, args.key, args.where, args.tiles, args.pcov, args.ptile)
        
        # create tld: SITENAME--KEY_DATATYPE_SUFFIX
        if args.notld:
            tld = args.outdir
        else:
            key = '' if args.key == '' else '--' + args.key
            suffix = '' if args.suffix == '' else '_' + suffix
            res = '' if args.res is None else '_%sx%s' % (args.res[0], args.res[1])
            bname = extents[0].site.LayerName() + key + res + '_' + args.command + suffix
            tld = os.path.join(args.outdir, bname)

        for extent in extents:
            inv = DataInventory(cls, extent, TemporalExtent(args.dates, args.days), **vars(args))
            datadir = os.path.join(tld, extent.site.Value())
            if inv.numfiles > 0:
                inv.mosaic(datadir=datadir, tree=args.tree, overwrite=args.overwrite, 
                           res=args.res, interpolation=args.interpolation, crop=args.crop)

    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Data Project error: %s' % e


if __name__ == "__main__":
    main()
