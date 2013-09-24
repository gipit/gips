#!/usr/bin/env python

import os,argparse,datetime
import gippy
import numpy, gdal

def add_options(subparser,parents=[]):
    parser = subparser.add_parser('index2prob', help='Convert Indices to Probability', 
        parents=parents, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_prob')
    group.add_argument('--min', help='Minimum index (prob=0)', default=0.25, type=float)
    group.add_argument('--max', help='Maximum index (prob=1)', default=0.75, type=float)

def process(infile, *args, **kwargs):
    start = datetime.datetime.now()
    fbase, ext = os.path.splitext(os.path.basename(infile))
    outfile = fbase + kwargs['suffix'] + ext
    img = gippy.Index2Probability(gippy.GeoImage(infile), outfile, kwargs['min'], kwargs['max'])
    print "%s -> %s: %s" % (fbase, os.path.basename(outfile), datetime.datetime.now()-start)