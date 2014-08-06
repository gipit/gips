#!/usr/bin/env python

import argparse
import gippy
from gipif.data.sar import SARData

from pdb import set_trace


def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='AGS Rice Detect', formatter_class=dhf)
    group = parser0.add_argument_group('inventory arguments')
    group.add_argument('-s', '--site', help='Vector file for region of interest', default=None)
    group.add_argument('-d', '--dates', help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    group.add_argument('--days', help='Include data within these days of year (doy1,doy2)', default=None)
    group.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
    group.add_argument('--%cov', dest='pcov', help='Threshold of %% coverage of tile over site', default=0, type=int)
    group.add_argument('--%tile', dest='ptile', help='Threshold of %% tile used', default=0, type=int)

    group = parser0.add_argument_group('algorithm options')
    group.add_argument('-o', '--output', help='Output file name', default='rice.tif')
    group.add_argument('--th0', help='Low threshold', default=-14.0, type=float)
    group.add_argument('--th1', help='High threshold', default=-8.0, type=float)
    group.add_argument('--dth0', help='Low threshold', default=90, type=int)
    group.add_argument('--dth1', help='High threshold', default=130, type=int)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetChunkSize(64.0)

    inv = SARData.inventory(site=args.site, dates=args.dates, days=args.days,
                            products=['sign'], pcov=args.pcov, ptile=args.ptile, sensors=['AWB1'])

    inv.printcalendar()
    inv.project(res=[100, 100])

    days = []
    for d in inv.dates:
        days.append((d-inv.dates[0]).days)

    images = inv.get_timeseries(product='sign')
    img = gippy.RiceDetect(images, args.output, days, args.th0, args.th1, args.dth0, args.dth1)

if __name__ == "__main__":
    main()
