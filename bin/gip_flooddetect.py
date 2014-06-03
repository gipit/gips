#!/usr/bin/env python

import argparse
import numpy
import gippy
from gipif.inventory import project_inventory
from gipif.data.sar import SARData
from gipif.data.landsat import LandsatData
from gipif.utils import VerboseOut
from pdb import set_trace

__version__ = '0.1.0'

def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='GIPIF Flood Detect', formatter_class=dhf)

    #group = parser0.add_argument_group('inventory arguments')
    parser0.add_argument('datadir', help='GIPIF Project directory')
    parser0.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
    parser0.add_argument('-p', '--product', help='Product to operate on', required=True)

    group = parser0.add_argument_group('algorithm options')
    group.add_argument('-o', '--output', help='Output file name', default='flood.tif')
    group.add_argument('--water', help='Threshold for water', default=-12.0, type=float)
    group.add_argument('--land', help='Threshold for land', default=-10.0, type=float)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)

    # inventory project directory
    inv = project_inventory(args.datadir)

    VerboseOut('Flood Detect v%s' % __version__)
    VerboseOut('Processing %s dates' % len(inv))
    th0 = args.water
    th1 = args.land

    # Input image(s)
    filenames = [inv[date][args.product] for date in inv.keys()]
    img = gippy.GeoImage(filenames)
    nodata = img[0].NoDataValue()

    # Flood detect algorithm
    imgout = gippy.GeoImage(args.output, img, gippy.GDT_Byte, 3)
    imgout[0].SetDescription('hits')
    imgout[0].SetDescription('num_observations')
    imgout[2].SetDescription('norm_hits')
    imgout.SetNoData(0)

    last_water = (img[0] < th0).Read()
    last_water[numpy.where(last_water == nodata)] = 0
    #dried = (img[0] > th1).Read()
    hits = numpy.zeros(last_water.shape)
    numobs = numpy.zeros(last_water.shape)
    for b in range(1, img.NumBands()):
        VerboseOut('Process band %s' % str(b+1), 2)
        numobs = numobs + img[b].DataMask()
        # Current dry land mask
        dried = (img[b] > th1).Read()
        dried[numpy.where(dried == nodata)] = 0
        # Find transitional pixels (water->dry) and increment
        inds = numpy.where(numpy.logical_and(dried, last_water))
        hits[inds] = hits[inds] + 1
        #imgout[b+1].Write(hits)
        # update water mask
        water = (img[b] < th0).Read()
        water[numpy.where(water == nodata)] = 0
        last_water = last_water + water
        # Reset water mask if it was dry
        last_water[numpy.where(dried > 0)] = 0
    imgout[0].Write(hits)
    imgout[1].Write(numobs)
    imgout[2].Write(numpy.divide(hits, numobs) * 50)

    imgout = None
    img = None

if __name__ == "__main__":
    main()
