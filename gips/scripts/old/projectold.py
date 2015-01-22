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

import os
import numpy
import pandas
import matplotlib.pyplot as plt
import gippy
from gippy.algorithms import BrowseImage
from gips.algorithms.core import Algorithm
from gips.utils import basename
from datetime import datetime


class Project(Algorithm):
    name = 'GIPS Project Utilities'
    __version__ = '1.0.0'

    def inventory(self, **kwargs):
        self.inv.pprint()

    def browse(self, quality=75, **kwargs):
        for date in self.inv.dates:
            for p in self.inv.products(date):
                try:
                    BrowseImage(self.inv[date].open(p), quality)
                except:
                    pass

    def stack(self, suffix='stack', **kwargs):
        """ Stack products (from single date) into single image file """
        for date in self.inv.dates:
            filenames = [self.inv[date].filenames[p] for p in self.inv.products(date)]
            img = gippy.GeoImage(filenames)
            bname = basename(filenames[0])
            bname = bname[0:bname.rfind('_', 0)]
            fout = os.path.join(self.inv.projdir, bname + '_' + suffix)
            imgout = img.Process(fout)
            imgout.CopyMeta(img)

    def plot(self, classes, **kwargs):
        """ Work in progress """
        from pdb import set_trace
        classimg = gippy.GeoImage(classes)
        dates = numpy.array([int(datetime.strftime(d, '%Y')) for d in self.inv.dates])
        dateset = set(dates)
        for p in self.inv.requested_products:
            for d in list(dateset):
                print 'Extracting pixels for %s' % d
                t1 = datetime.now()
                # find matching dates
                gimg = self.inv.get_timeseries(p, dates=numpy.array(self.inv.dates)[dates == d])
                arr = gimg.Extract(classimg[0])
                df = pandas.DataFrame(gimg.Extract(classimg[0])).T
                df[df == gimg[0].NoDataValue()] = numpy.nan
                print 'extract time %s' % (datetime.now() - t1)
                df.plot(kind='box')
                plt.show()
                set_trace()

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

        #parser = subparser.add_parser('plot', help='Create plots for products', **kwargs)
        #parser.add_argument('--classes', help='Classes/Pixels-Of-Interest image', default='', required=True)
        #parser.add_argument('--dates', help='Group by these dates', nargs='*', default='')

        return parser0


def main():
    Project.main()
