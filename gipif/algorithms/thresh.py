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
import numpy, gdal

def add_options(subparser,parents=[]):
    parser = subparser.add_parser('thresh', help='Threshold image', 
        parents=parents, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_thresh')
    group.add_argument('-t','--threshold', help='Threshold value (> threshold is true)', default=0.5, type=float)
    group.add_argument('-b', '--band', help='Band to threshold', default=1, type=int)

def process(infile, *args, **kwargs):
    fbase, ext = os.path.splitext(os.path.basename(infile))
    start = datetime.datetime.now()
    outfile = fbase + kwargs['suffix'] + str(kwargs['threshold']) + ext

    imgout = gippy.GeoImage(outfile, gippy.GeoImage(infile), gippy.GDT_Byte, 1)
    nodatavalue = imgout[0].NoDataValue()
    imgarr = numpy.zeros((imgout.YSize(), imgout.XSize()))
    imgout = None
    # Open input and threshold
    fh = gdal.Open(infile)
    fhb = fh.GetRasterBand(kwargs['band'])
    img = fhb.ReadAsArray()
    (i,j) = numpy.where(img > kwargs['threshold'])
    imgarr[i,j] = 1
    (i,j) = numpy.where(img == nodatavalue)
    imgarr[i,j] = 0
    # Write output
    fhout = gdal.Open(outfile,gdal.GA_Update)
    fhout.GetRasterBand(1).WriteArray(imgarr)
    fhout.GetRasterBand(1).SetNoDataValue(0)

    print "%s -> %s: %s" % (fbase, os.path.basename(outfile), datetime.datetime.now()-start)