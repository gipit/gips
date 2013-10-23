#!/usr/bin/env python

import argparse
import gippy

def add_options(subparser,parents=[]):
	parser = subparser.add_parser('fmask',help='FMask cloud masking algorithm',parents=parents,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	group = parser.add_argument_group('algorithm arguments')
	group.add_argument('--tolerance', help='Tolerance (1-5). Higher tolerance means fewer clouds', default=3, type=int)
	group.add_argument('--dilate', help='Size of dilation filter', default=10, type=int)
	#group.add_argument('--shadow', help='Shadow threshold', default=0.02, type=float)

def process(image, outfile, tolerance=3, dilate=10, verbose=0, **kwargs):
	gippy.Options.SetVerbose(verbose)
	img = gippy.Fmask(image, outfile, tolerance, dilate)
	filename = img.Filename()
	del img
	meta = gippy.landsat.readmeta(image.Filename())
	from gippy.algorithms.acloud import AddShadowMask
	AddShadowMask(filename, 90.0 - meta['solarzenith'], meta['solarazimuth'], 4000, 2)
	return gippy.GeoImage(filename)