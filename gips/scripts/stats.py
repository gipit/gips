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
    title = Colors.BOLD + 'GIPS Image Statistics (v%s)' % __version__ + Colors.OFF

    parser0 = GIPSParser(datasources=False, description=title)
    parser0.add_default_parser()
    parser0.add_projdir_parser()
    group = parser0.add_argument_group('masking options')
    args = parser0.parse_args()

    # TODO - check that at least 1 of filemask or pmask is supplied

    try:
        print title
        header = ['min', 'max', 'mean', 'sd', 'skew', 'count']

        for projdir in args.projdir:
            VerboseOut('Stats for Project directory: %s' % projdir, 1)
            inv = ProjectInventory(projdir, args.products)
        
            files = {}
            for date in inv.dates:
                VerboseOut('Calculating statistics for %s' % date)
                for p in inv.products(date):
                    img = inv[date].open(p)
                    if p not in files.keys():
                        files[p] = open(os.path.join(projdir, p + '_stats.txt'), 'w')
                        # write header
                        files[p].write('date ')
                        if img.NumBands() == 1:
                            files[p].write(' '.join(header))
                        else:
                            for band in img:
                                h = [band.Description() + "-" + a for a in header]
                                files[p].write(' '.join(h) + ' ')
                        files[p].write('\n')
                    # print date and stats
                    files[p].write(date.strftime('%Y-%j'))
                    for band in img:
                        stats = band.Stats()
                        [files[p].write(' ' + str(s)) for s in stats]
                        files[p].write(' ')
                    files[p].write('\n')
                    img = None
            for f in files:
                files[f].close()

    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'Error: %s' % e


if __name__ == "__main__":
    main()
