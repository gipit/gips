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
    parser = subparser.add_parser('permutations', help='Calculate all permutations of bands/files', 
        parents=parents, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('--file',help='Files to use in permutation', required=True)
    group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_permutation')

def process(infile, *args, **kwargs):
    start = datetime.datetime.now()
    fbase, ext = os.path.splitext(os.path.basename(infile))
    outfile = fbase + kwargs['suffix'] + ext
    img = gippy.Permutations(gippy.GeoImage(infile), gippy.GeoImage(kwargs['file']), outfile)
    print "%s -> %s: %s" % (fbase, os.path.basename(outfile), datetime.datetime.now()-start)