#!/usr/bin/env python

import os, sys, argparse
import gdal
import numpy
import gippy
import datetime

if __name__ == "__main__":
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(prog=os.path.split(__file__)[1], formatter_class=dhf,
                                      description='Data analysis')
    #subparser = parser0.add_subparsers(dest='command')

    parser0.add_argument('files', nargs='*', help='Imagery files to analyze')
    parser0.add_argument('-t','--threshold', help='Difference threshold between LSWI and NDVI', default=0.5, type=float)
    parser0.add_argument('-s','--suffix', help='Append suffix to filename for output', default='thresh')
    #parser0.add_argument('-v','--vector',help='Shapefile of site to be extracted from imagery files', required=True)
    #parser0.add_argument('-w','--warp', help='Reproject raster files to match vector', action='store_true', default=False)
    #parser0.add_argument('-o', '--overwrite', help='Overwrite existing files', action='store_true', default=False)

    args = parser0.parse_args()

    for f in args.files:
        start = datetime.datetime.now()
        fbase,ext = os.path.splitext(os.path.basename(f))
        # Create output file
        fout = fbase + args.suffix
        imgout = gippy.GeoImage(fout, gippy.GeoImage(f), gippy.GDT_Byte, 1)
        imgarr = numpy.zeros((imgout.YSize(), imgout.XSize()))
        fout = imgout.Filename()
        imgout = None
        # Open input and threshold
        fh = gdal.Open(f)
        img1 = fh.GetRasterBand(1).ReadAsArray()
        img2 = fh.GetRasterBand(3).ReadAsArray()
        img = img1 - img2
        (i,j) = numpy.where(img < args.threshold)
        imgarr[i,j] = 1
        (i,j) = numpy.where(img1 == fh.GetRasterBand(1).GetNoDataValue())
        imgarr[i,j] = 0
        # Write output
        fhout = gdal.Open(fout,gdal.GA_Update)
        fhout.GetRasterBand(1).WriteArray(imgarr)
        fhout.GetRasterBand(1).SetNoDataValue(0)

        # update to get no data value from source....for that matter do this whole thing with API not cmdline utils
        #if args.overwrite: cmd = cmd + ' -overwrite'
        #if os.path.exists(fout) and args.overwrite: os.remove(fout)
        #out = commands.getstatusoutput(cmd)
        print '%s: %s' % (os.path.basename(f), datetime.datetime.now()-start)

