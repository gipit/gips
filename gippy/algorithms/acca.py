#!/usr/bin/env python

import os, argparse
import gippy

def add_options(subparser,parents=[]):
    parser = subparser.add_parser('acca',help='ACCA (Automatic Cloud Cover Assessment)',parents=parents,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_acca')
    #group.add_argument('--tolerance', help='Tolerance (1-5). Higher tolerance means fewer clouds', default=3, type=int)
    #group.add_argument('--dilate', help='Size of dilation filter', default=10, type=int)
    #group.add_argument('--shadow', help='Shadow threshold', default=0.02, type=float)

def process(f, *args, **kwargs):
    fbase,ext = os.path.splitext(os.path.basename(f))
    fout = fbase + kwargs['suffix']
    return gippy.ACCA(gippy.GeoImage(f), fout)