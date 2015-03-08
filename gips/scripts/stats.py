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

import gippy
from gips.parsers import GIPSParser
from gips.inventory import ProjectInventory
from gips.utils import Colors, VerboseOut, basename

__version__ = '0.1.0'

def main():
    title = Colors.BOLD + 'GIPS Image Statistics v%s' % __version__ + Colors.OFF

    parser0 = GIPSParser(datasources=False, description=title)
    parser0.add_default_parser()
    parser0.add_projdir_parser()
    group = parser0.add_argument_group('masking options')
    args = parser0.parse_args()

    # TODO - check that at least 1 of filemask or pmask is supplied

    try:
        print title
        inv = ProjectInventory(args.projdir, args.products)

        for date in inv.dates:
            VerboseOut('Calculating statistics for %s' % date)
            for p in inv.products(date):
                # don't mask any masks
                if p in available_masks:
                    continue
                meta = ''
                update = True if args.original else False
                img = inv[date].open(p, update=update)
        for band in img:
            stats = band.Statistics()

                if args.filemask is not None:
                    img.AddMask(mask_file[0])
                    meta = basename(args.filemask) + ' '
                for mask in available_masks:
                    img.AddMask(inv[date].open(mask)[0])
                    meta = meta + basename(inv[date][mask]) + ' '
                if meta != '':
                    if args.original:
                        VerboseOut('  %s' % (img.Basename()), 2)
                        img.Process()
                        img.SetMeta('MASKS', meta)
                    else:
                        fout = os.path.splitext(img.Filename())[0] + args.suffix + '.tif'
                        if not os.path.exists(fout) or overwrite:
                            VerboseOut('  %s -> %s' % (img.Basename(), basename(fout)), 2)
                            imgout = img.Process(fout)
                            imgout.SetMeta('MASKS', meta)
                            imgout = None
                img = None
        mask_file = None            

    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Masking error: %s' % e


if __name__ == "__main__":
    main()
