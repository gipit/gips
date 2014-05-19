#!/usr/bin/env python

import argparse
import numpy
import gippy
from gipif.data.sar import SARData
from gipif.utils import VerboseOut
from pdb import set_trace

__version__ = '0.1.0'

def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='GIPIF Flood Detect', formatter_class=dhf)
    group = parser0.add_argument_group('inventory arguments')
    group.add_argument('-s', '--site', help='Vector file for region of interest', default=None)
    group.add_argument('-d', '--dates', help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    group.add_argument('--days', help='Include data within these days of year (doy1,doy2)', default=None)
    group.add_argument('--sensors', help='Sensors to include', nargs='*', default=None)
    group.add_argument('--%cov', dest='pcov', help='Threshold of %% coverage of tile over site', default=0, type=int)
    group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)

    group = parser0.add_argument_group('misc options')
    group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
    group.add_argument('--chunksize', help='Chunk size in MB', default=512.0)

    group = parser0.add_argument_group('algorithm options')
    group.add_argument('-o', '--output', help='Output file name', default='rice.tif')
    group.add_argument('--water', 'Sigma nought threshold for water', default=-12.0, type=float)
    group.add_argument('--land', 'Sigma nought threshold for land', default=-10.0, type=float)
    group.add_argument('--mean', help='Mean water signal', default=-12.0, type=float)
    group.add_argument('--sd', help='Water signal std dev', default=1.0, type=float)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetChunkSize(args.chunksize)

    inv = SARData.inventory(site=args.site, dates=args.dates, days=args.days,
                            products=['sign'], pcov=args.pcov, ptile=args.ptile, sensors=args.sensors)
    inv.print_inv()
    # TODO - get resolution from sensor
    inv.project(res=[100, 100])

    days = []
    for d in inv.dates:
        days.append((d-inv.dates[0]).days)
    imgs = inv.get_timeseries(product='sign')
    nodata = imgs[0].NoDataValue()

    VerboseOut('Flood Detect v%s' % __version__)
    VerboseOut('Processing %s dates' % len(inv.dates))

    # Flood detect algorithm
    imgout = gippy.GeoImage(args.output, imgs, gippy.GDT_Byte, 3)
    imgout.SetNoData(0)
    th0 = args.water
    th1 = args.land
    last_water = (imgs[0] < th0).Read()
    last_water[numpy.where(last_water == nodata)] = 0
    #dried = (imgs[0] > th1).Read()
    hits = numpy.zeros(last_water.shape)
    numobs = numpy.zeros(last_water.shape)
    for b in range(1, imgs.NumBands()):
        VerboseOut('Process band %s' % b, 2)
        numobs = numobs + imgs[b].DataMask()
        # Current dry land mask
        dried = (imgs[b] > th1).Read()
        dried[numpy.where(dried == nodata)] = 0
        # Find transitional pixels (water->dry) and increment
        inds = numpy.where(numpy.logical_and(dried, last_water))
        hits[inds] = hits[inds] + 1
        #imgout[b+1].Write(hits)
        # update water mask
        water = (imgs[b] < th0).Read()
        water[numpy.where(water == nodata)] = 0
        last_water = last_water + water
        # Reset water mask if it was dry
        last_water[numpy.where(dried > 0)] = 0
    imgout[0].Write(hits)
    imgout[1].Write(numobs)
    imgout[2].Write(numpy.divide(hits, numobs) * 50)

    imgout = None
    imgs = None

if __name__ == "__main__":
    main()
