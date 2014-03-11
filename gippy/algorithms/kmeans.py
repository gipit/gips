#!/usr/bin/env python
################################################################################
#    GIPPY: Geospatial Image Processing library for Python
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

import os,argparse,datetime
import gippy

def add_options(subparser,parents=[]):
    parser = subparser.add_parser('kmeans',help='K-means classification aglorithm', 
        parents=parents, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_kmeans')
    group.add_argument('--classes', help='Number of classes in image', default=5, type=int)
    group.add_argument('--iter', help='Maximum number of iterations', default=5, type=int)
    group.add_argument('--thresh', help='Pixel Change Threshold (%%)', default=5, type=int)
    group.add_argument('--multclass', help='Multiple runs using from 3 to 12 classes', default=False, action='store_true')

def process(infile, *args, **kwargs):
    fbase, ext = os.path.splitext(os.path.basename(infile))
    if kwargs['multclass']:
        classes = range(3,13)
    else:
        classes = [kwargs['classes']]
    for c in classes:
        start = datetime.datetime.now()
        outfile = fbase + kwargs['suffix'] + str(c) + ext
        img = gippy.kmeans(gippy.GeoImage(infile), outfile, c, kwargs['iter'], kwargs['thresh'])
        print "%s -> %s: %s" % (fbase, os.path.basename(outfile), datetime.datetime.now()-start)