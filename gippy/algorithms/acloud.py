#!/usr/bin/env python

import os
import argparse
import datetime
import gippy

# Defaults
_cloud = 4000
_dilate = 10
_minred = 0.2
_maxtemp = 14
_maxndvi = 0.2


# for shadowmask
from osgeo import gdal
import numpy
from scipy import weave
from scipy.weave import converters

#Issues:
#1. use CLOUD_COVER from MTL file, create a histogram from the data (Therm) and estimate the threshold.
#2. think up an automated way to estimate cloud height. For unknown or multiple cloud heights, use a 
#   larger dilation parameter 

# absolute filename bad, should at least come from landsatlib file
vangname = '/titan/data/landsat/documents/StandardLandsatViewAngle.tif'

def AddShadowMask(inputname, se_degrees, sa_degrees, cloudheight, bandnum=1):

    #required input is 1. a full cloudmask, 2. a full map of view angles in radians,
    #3. an output cloud shadowmask, 4. sunelevation in radians, 5. solarazimuth in radians,
    #6. cloudheight in meters, #7. nx, #8. ny. 
    code = """
        int i, j, sj, si;
        float shiftdistance, shiftdistance_x, shiftdistance_y;
        float sinoforbitangle, cosoforbitangle;
        float distance, yshift, xshift;
        sinoforbitangle = 0.1582;
        cosoforbitangle = 0.9874;
        distance = cloudheight/tan(sunelevation);
        yshift = cos(solarazimuth) * distance;
        xshift = -1.0 * sin(solarazimuth) * distance;
        for (j=0;j<ny;j++) {
            for (i=0;i<nx;i++) {
                if (inmask(j,i) == 1) {
                    //calculate shift for actual cloud location
                    shiftdistance = cloudheight/tan(viewangle(j,i));
                    shiftdistance_x = shiftdistance * cosoforbitangle;
                    shiftdistance_y = shiftdistance * sinoforbitangle;
                    //those on the west side will have positive angles/shiftdistances
                    //those on the east side will have negative angles/shiftdistances
                    //so ADD the shiftdistances (for x i'm sure, not sure for y)
                    sj = j+(shiftdistance_y+yshift)/30.0;
                    si = i+(shiftdistance_x+xshift)/30.0;
                    sj = j+(yshift)/30.0;
                    si = i+(xshift)/30.0;
                    if (sj > 0 && si > 0 && sj < ny && si < nx) {
                        outmask(sj,si) = 1;
                    }
                }
            }
        }
    #"""

    gdal.UseExceptions()

    #open full cloudmask
    
    fo_in=gdal.Open(inputname, gdal.GA_Update)
    tband = fo_in.GetRasterBand(1)
    inmask=tband.ReadAsArray()
    proj = fo_in.GetProjection()
    gt = fo_in.GetGeoTransform()
    nx = fo_in.RasterXSize
    ny = fo_in.RasterYSize
    xcenter = int(numpy.floor(float(nx)/2.0))
    ycenter = int(numpy.floor(float(ny)/2.0))

    #open and cut viewangle layer (in radians)
    v_in=gdal.Open(vangname)
    vangle=v_in.ReadAsArray().astype(float)
    vanx=v_in.RasterXSize
    vany=v_in.RasterYSize 
    v_in = None
    MidPoint=int(numpy.floor(float(vanx)/2.0))
    addvalx = 0
    if numpy.mod(vanx,2) != 0:
        addvalx = 1 
    addvaly = 0
    if numpy.mod(vany,2) != 0:
        addvaly = 1 
    ystart=MidPoint-ycenter
    yend=MidPoint+ycenter+addvaly
    ysize=yend-ystart
    xstart=MidPoint-xcenter
    xend=MidPoint+xcenter+addvalx
    xsize=xend-xstart
    viewangle=numpy.zeros((ny,nx),float)
    viewangle[:,:]=vangle[ystart:yend,xstart:xend]

    #create cloud shadow mask
    #convert sunelevation and solarazimuth to radians
    #estimate cloud height
    outmask=numpy.zeros((ny,nx))
    sunelevation=se_degrees*numpy.pi/180.0
    solarazimuth=sa_degrees*numpy.pi/180.0
#    cloudheight=4000.0
#    if len(sys.argv) > 4:
#        cloudheight=float(sys.argv[4])

    weave.inline(code, ['inmask', 'outmask', 'viewangle', 'sunelevation', 'solarazimuth', 'cloudheight', 'nx', 'ny'], \
       support_code='#include<cmath>', type_converters=converters.blitz, verbose=0)

    (i,j)=numpy.where(outmask == 1)
    outmask[i,j] = 2
    (i,j)=numpy.where(inmask == 1)
    outmask[i,j] = 1

    #outputname = inputname.split('.tif')[0]+'_wshad.tif'
    #FORMAT = 'GTiff'
    #DATATYPE = gdal.GDT_Byte
    #OPTIONS = []
    #XSIZE = nx
    #YSIZE = ny
    #NBANDS = 1
    #driver = gdal.GetDriverByName(FORMAT)
    #tfh = driver.Create(outputname, XSIZE, YSIZE, NBANDS, DATATYPE, OPTIONS)
    #tfh.SetProjection(proj)
    #tfh.SetGeoTransform(gt)

    tband = fo_in.GetRasterBand(bandnum)
    tband.WriteArray(outmask)
 
    fo_in = None
    #tfh = None

# Shagen's Algorithm processing given list of files
def add_options(subparser,parents=[]):
    parser = subparser.add_parser('acloud',help='AutoCloud Mask Algorithm',parents=parents,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('--cloud', help='Cloud height in meters', default=_cloud, type=int)
    group.add_argument('--minred', help='Minimum threshold for red band reflectivity', default=_minred, type=float)
    group.add_argument('--maxtemp', help='Max temperature threshold in Celsius (thermal band)', default=_maxtemp, type=float)
    group.add_argument('--maxndvi', help='Max threshold for NDVI', default=_maxndvi, type=float)
    group.add_argument('--dilate', help='Number of pixels for morphological dilation', default=_dilate, type=int)

def process(image, outfile, cloud=_cloud, minred=_minred, maxtemp=_maxtemp, maxndvi=_maxndvi, dilate=_dilate, **kwargs):
    img = gippy.AutoCloud(image, outfile, cloud, minred, maxtemp, maxndvi, dilate)
    img.ClearNoData()
    # metadata should be put into image (img.Meta(metaname))
    # replace with full C++ implementation
    # also assuming landsat here
    filename = img.Filename()
    del img
    from agspy.data.landsatlib import readmtl
    meta = readmtl(image.Filename())
    AddShadowMask(filename, 90.0 - meta['solarzenith'], meta['solarazimuth'], _cloud)

    #import pdb
    #pdb.set_trace()
    return gippy.GeoImage(filename)