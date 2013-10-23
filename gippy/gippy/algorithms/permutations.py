#!/usr/bin/env python

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