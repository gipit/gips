#!/usr/bin/env python

import os, sys, argparse
import gdal
import numpy
import gippy
import agspy.utils.zonal as zonal

if __name__ == "__main__":
    prog = os.path.split(__file__)[1]
    parser0 = argparse.ArgumentParser(prog=prog, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Class stats')
    #subparser = parser0.add_subparsers(dest='command')

    parser0.add_argument('files', nargs='*', help='Class maps to analyze')
    parser0.add_argument('-v','--vector',help='Count pixels in each class for provided shapefile training data')
    parser0.add_argument('-f','--field', help='The field in vector file to use as training set', default='FID')
    #parser0.add_argument('-n','--nodata', help='This is a nodata value, discard from total pixels', default=0, type=int)

    args = parser0.parse_args()

    nodata = 0

    for f in args.files:
        fh = gdal.Open(f)
        img = fh.GetRasterBand(1).ReadAsArray()
        sz = img.shape
        sz = sz[0]*sz[1] 
        numnodata = len(numpy.where(img == nodata)[0])

        print os.path.basename(f) + ':'
        print 'NoData (0): %4.2f%% of scene' % (float(numnodata)/sz*100)
        sz = sz - numnodata

        numclass = int(numpy.max(img))

        head = '{0:^5} {1:^6}'.format('Class','Scene%')
        if args.vector:
            #dbftable = zonal.get_dbftable(args.vector)
            #attrmap = dbftable.create_map('FID',args.field)
            rvector = zonal.rasterize_vector(args.vector, f, field=args.field)
            fh_vec = gdal.Open(rvector)
            imgv = fh_vec.GetRasterBand(1).ReadAsArray()
            numclass_vec = numpy.max(imgv)
            totalpix_vec = []
            for iclass in range(1,numclass_vec+1): 
                head = head + ' {0:^6}'.format('TS%s'%iclass)
                totalpix_vec.append( len(numpy.where(imgv == iclass)[0]) )

        print head
        for i in range(1,numclass+1):
            inds = numpy.where(img == i)
            num = float(len(inds[0]))
            line = '{i:^5} {p:^6.2f}'.format(i=i,p=num/sz*100)
            if args.vector:
                img_vec_sub = imgv[inds]
                for iTS in range(1,numclass_vec+1):
                    vinds = numpy.where(img_vec_sub == iTS)
                    numpix = float(len(vinds[0]))
                    line = line + ' {ts:^6.1f}'.format(ts=numpix/totalpix_vec[iTS-1]*100)
            if num > 0:
                print line


        fh = None
        fhv = None

            


