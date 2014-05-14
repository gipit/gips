#!/usr/bin/env python

import argparse
import numpy
import gippy
from gipif.data.sar import SARData
from pdb import set_trace


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

    print 'dates = ', len(inv.dates)

    # Flood detect algorithm
    imgout = gippy.GeoImage(args.output, imgs, gippy.GDT_Byte, 1)
    th0 = args.mean + args.sd
    th1 = args.mean + 3 * args.sd
    last_water = (imgs[0] > th0).Read()
    #dried = (imgs[0] > th1).Read()
    hits = numpy.zeros(last_water.shape)
    for b in range(1, imgs.NumBands()):
        print b
        water = (imgs[b] > th0).Read()
        dried = (imgs[b] > th1).Read()
        print 'before where'
        inds = numpy.where(numpy.logical_and(dried, last_water))
        print 'after where'
        hits[inds] = hits[inds] + 1
        last_water = water
    imgout[0].Write(hits)

    imgout = None
    imgs = None

    set_trace()


if __name__ == "__main__":
    main()
