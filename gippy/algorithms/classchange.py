#!/usr/bin/env python

def add_options(subparser,parents=[]):
	parser = subparser.add_parser('classmod',help='Classification Modifier', parents=[gparser], formatter_class=hformat)
	group = parser.add_argument_group('algorithm arguments')
	group.add_argument('image1', help='Input Image')
	group.add_argument('c2c', help='Change class from A to B when operator is true (A,B)')
	group.add_argument('image2', help='Operand Image')
	group.add_argument('op', help='Comparison operator (lt, gt, eq)')
	group.add_argument('th', help='Threshold used with comparison operator')

def process(image, outfile, cloud=cloud, shadow=shadow, verbose=verbose):
	return outfile