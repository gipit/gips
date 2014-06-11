#!/usr/bin/env python

import argparse
import numpy
import gippy
from gipif.inventory import project_inventory
from gipif.utils import VerboseOut

__version__ = '0.7.0'


def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='GIPIF Flood Detect', formatter_class=dhf)

    #group = parser0.add_argument_group('inventory arguments')
    parser0.add_argument('datadir', help='GIPIF Project directory', default='./')
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
    dates = sorted(inv.keys())
    filenames = [inv[date][args.product] for date in dates]
    img = gippy.GeoImage(filenames)
    nodata = img[0].NoDataValue()

    # Flood detect algorithm
    imgout = gippy.GeoImage(args.output, img, gippy.GDT_Byte, img.NumBands()+2)
    imgout[0].SetDescription('hits')
    imgout[1].SetDescription('num_observations')
    for b in range(0, img.NumBands()):
        imgout[b+2].SetDescription(str(dates[b]))
    imgout.SetNoData(0)
    VerboseOut('Processing %s' % dates[0], 2)
    if th0 > th1:
        th0 = -th0
        th1 = -th1
        for b in range(0, img.NumBands()):
            img[b] = img[b] * -1

    last_water = (img[0] < th0).Read()
    last_water[numpy.where(last_water == nodata)] = 0
    hits = numpy.zeros(last_water.shape)
    numobs = numpy.zeros(last_water.shape)
    days = last_water
    imgout[2].Write(days)
    for b in range(1, img.NumBands()):
        VerboseOut('Processing %s' % dates[b], 2)
        days = days + (last_water * (dates[b]-dates[0]).days)
        # Increment # of observations
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
        # add new water regions (1 day)
        nowater_inds = numpy.where(days == 0)
        days[nowater_inds] = water[nowater_inds]
        last_water = numpy.minimum(last_water + water, 1)
        # Reset water mask if it was dry
        dry_inds = numpy.where(dried > 0)
        last_water[dry_inds] = 0
        days[dry_inds] = 0
        imgout[b+2].Write(days)
    imgout[0].Write(hits)
    imgout[1].Write(numobs)

    imgout = None
    img = None

if __name__ == "__main__":
    main()
