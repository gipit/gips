#!/usr/bin/env python

import argparse
import gippy

def add_options(subparser,parents=[]):
	parser = subparser.add_parser('acca',help='ACCA (Automatic Cloud Cove Assessment)',parents=parents,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#group = parser.add_argument_group('algorithm arguments')
	#group.add_argument('--tolerance', help='Tolerance (1-5). Higher tolerance means fewer clouds', default=3, type=int)
	#group.add_argument('--dilate', help='Size of dilation filter', default=10, type=int)
	#group.add_argument('--shadow', help='Shadow threshold', default=0.02, type=float)

def process(image, outfile, **kwargs):
	return gippy.ACCA(image, outfile)