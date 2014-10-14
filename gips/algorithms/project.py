#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    Copyright (C) 2014 Matthew A Hanson
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
from gips.algorithms.core import Algorithm
from gips.utils import basename


class Project(Algorithm):
    name = 'GIPS Project Utilities'
    __version__ = '1.0.0'

    def inventory(self, **kwargs):
        self.inv.pprint()

    def browse(self, quality=75, **kwargs):
        for date in self.inv.dates:
            for p in self.inv.product_list(date):
                try:
                    gippy.BrowseImage(self.inv[date].open(p), quality)
                except:
                    pass

    def stack(self, suffix='stack', **kwargs):
        """ Stack products (from single date) into single image file """
        for date in self.inv.dates:
            filenames = [self.inv[date].filenames[p] for p in self.inv.product_list(date)]
            img = gippy.GeoImage(filenames)
            bname = basename(filenames[0])
            bname = bname[0:bname.rfind('_', 0)]
            fout = os.path.join(self.inv.projdir, bname + '_' + suffix)
            imgout = img.Process(fout)
            imgout.CopyMeta(img)

    @classmethod
    def parser(cls, parser0):
        subparser, kwargs = cls.subparser(parser0, project=True)

        parser = subparser.add_parser('browse', help='Create browse imagery', **kwargs)
        parser.add_argument('-q', '--quality', help='JPG Quality', default=75)

        parser = subparser.add_parser('inventory', help='Print project inventory', **kwargs)

        parser = subparser.add_parser('stack', help='Stack given products into one image (per date)', **kwargs)
        parser.add_argument('--suffix', help='Suffix to give stacked output image', default='stack')
        #h = 'Create as GDAL virtual file format (VRT)'
        #parser.add_argument('--vrt', help=h, default=False, action='store_true')

        return parser0


def main():
    Project.main()
