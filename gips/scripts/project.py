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
    title = Colors.BOLD + 'GIPS Data Project Utility v%s' % gipsversion + Colors.OFF

    # argument parsing
    parser0 = GIPSParser(description=title)
    parser = inventory_parser()
    group = parser.add_argument_group('project options')
    group.add_argument('--res', nargs=2, help='Resolution of (warped) output rasters', default=None, type=float)
    group.add_argument('--crop', help='Crop down to minimum bounding box', default=False, action='store_true')
    group.add_argument('--nowarp', help='Do not warp (or crop)', default=False, action='store_true')
    group.add_argument('--interpolation', help='Interpolate using: 0-NN, 1-Bilinear, 2-Cubic', default=0, type=int)
    group.add_argument('--nomosaic', help='Do not mosaic (keep as tiles)', default=False, action='store_true')
    group.add_argument('--datadir', help='Directory to store output', default=None)
    group.add_argument('--suffix', help='Suffix on end of project directory', default=None)
    group.add_argument('--format', help='Format for output file', default="GTiff")
    group.add_argument('--chunksize', help='Chunk size in MB', type=float, default=512.0)
    parser0.add_data_sources(parents=[parser])
    args = parser0.parse_args()

    try:
        print title
        cls = data_class(args.command)
        inv = cls.inventory(**vars(args))
        del args.products
        inv.project(**vars(args))
    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Data Project error: %s' % e


if __name__ == "__main__":
    main()
