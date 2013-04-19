#!/usr/bin/env python

import argparse
import gippy

def add_options(subparser,parents=[]):
	parser = subparser.add_parser('fmask',help='FMask cloud masking algorithm',parents=parents,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	group = parser.add_argument_group('algorithm arguments')
	group.add_argument('--cloud', help='Cloud probability (?)', default=50, type=int)
	group.add_argument('--shadow', help='Shadow threshold', default=0.02, type=float)

def process(image, outfile, cloud=50, shadow=0.2, verbose=0, **kwargs):
	img = gippy.Fmask(image, outfile)
	return img