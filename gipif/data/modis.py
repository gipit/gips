#!/usr/bin/env python
################################################################################
#    GIPPY: Geospatial Image Processing library for Python
#
#    Copyright (C) 2014 Matthew A Hanson
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>
################################################################################

import sys
import os
import re
import time
import datetime
import urllib
from osgeo import gdal
from collections import OrderedDict

import math
import numpy as np

import gippy
from gippy import GeoImage
from gippy.data.core import Repository, Asset, Data
from gippy.data.inventory import DataInventory
from gippy.utils import File2List, List2File, VerboseOut

from pdb import set_trace


def binmask(arr, bit):
    """ Return boolean array indicating which elements as binary have a 1 in
        a specified bit position. Input is Numpy array.
    """
    return arr & (1<<(bit-1)) == (1<<(bit-1))


class ModisRepository(Repository):

    _rootpath = '/titan/data/modis'
    _stagedir = os.path.join(_rootpath, 'stage')

    _tiles_vector = 'modis_sinusoidal_grid_world.shp'
    #_tilesdir = 'tiles.dev'

    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")
        h = str(int(feature.GetField(fldindex_h))).zfill(2)
        v = str(int(feature.GetField(fldindex_v))).zfill(2)
        tile = "h%sv%s" % (h, v)
        return tile


class ModisAsset(Asset):
    Repository = ModisRepository

    _sensors = {
        'MOD': {'description': 'Terra'},
        'MYD': {'description': 'Aqua'},
        'MCD': {'description': 'Combined'}
    }

    _assets = {
        'MOD10A1': {
            'pattern': 'MOD10A1*hdf',
            'url': 'ftp://n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005',
            'startdate': datetime.date(2000, 2, 24),
            'latency': -3
        },
        'MOD11A1': {
            'pattern': 'MOD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005',
            'startdate': datetime.date(2000, 3, 5),
            'latency': -1
        },
        'MYD11A1': {
            'pattern': 'MYD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.005',
            'startdate': datetime.date(2002, 7, 8),
            'latency': -1
        }
    }

    # Not used anywhere
    _launchdate = {
        'MOD': datetime.date(1999, 12, 1),
        'MYD': datetime.date(2002, 5, 1),
        'MCD': None
    }

    # Should this be specified on a per asset basis?
    _defaultresolution = [926.625433138333392, -926.625433139166944]

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(ModisAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        self.asset = bname[0:7]
        self.tile = bname[17:23]
        year = bname[9:13]
        doy = bname[13:16]

        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()
        self.sensor = bname[:3]

        datafiles = self.datafiles()

        # I don't understand what this is used for
        # because the asset could be many different things
        # self.products = {'sds1': datafiles[0]}

    def datafiles(self):
        indexfile = self.filename + '.index'

        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:
            gdalfile = gdal.Open(self.filename)
            subdatasets = gdalfile.GetSubDatasets()
            datafiles = [s[0] for s in subdatasets]
            List2File(datafiles, indexfile)
        return datafiles

    @classmethod
    def fetch(cls, asset, tile, date):
        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)

        httploc = cls._assets[asset]['url']

        year, month, day = date.timetuple()[:3]
        doy = date.timetuple()[7]

        pattern = ''.join(['(', asset, '.A', str(year), str(doy).zfill(3), '.', tile, '.005.\d{13}.hdf)'])
        mainurl = ''.join([httploc, '/', str(year), '.', '%02d' % month, '.', '%02d' % day])

        VerboseOut('%s: mainurl %s, pattern %s' % (asset, mainurl, pattern), 4)

        try:
            VerboseOut('opening listing', 1)
            listing = urllib.urlopen(mainurl).readlines()
        except Exception, e:
            listing = None
            VerboseOut('unable to access %s' % mainurl, 1)
            return 2

        cpattern = re.compile(pattern)
        name = None
        success = False

        for item in listing:
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                name = cpattern.findall(item)[0]
                VerboseOut('found %s in %s' % (name, item.strip()), 4)
                url = ''.join([mainurl, '/', name])
                VerboseOut('the url is %s' % url, 4)
                try:
                    urllib.urlretrieve(url, os.path.join(cls.Repository.spath(), name))
                    VerboseOut('retrieved %s' % name, 4)
                    success = True
                except Exception, e:
                    VerboseOut('unable to retrieve %s from %s' % (name, url), 4)

        if not success:
            VerboseOut('did not find a match for %s in listing of %s' % (pattern, mainurl), 4)
            return 1
        else:
            return 0


class ModisData(Data):
    """ A tile of data (all assets and products) """
    name = 'Modis'
    Asset = ModisAsset
    _pattern = '*.tif'
 
    _products = {
        'temp': {
            'description': 'Surface temperature data',
            'assets': ['MYD11A1', 'MOD11A1']
        },
        'snow': {
            'description': 'Snow and ice cover data',
            'assets': ['MOD10A1']
        }
    }
    


    def process(self, products):

        # modis process --temp -t h12v08 -d 2010-001 -v4 --overwrite

        # print "self._products"
        # print self._products

        # print "products"
        # print products

        platformnames = {'MOD':'Terra', 'MYD':'Aqua'}

        for key, val in products.items():
        
            # print "key, val"
            # print key, val

            outfname = os.path.join(self.path, self.basename + '_' + key)        
            VerboseOut("outfname: %s" % outfname, 4)

            # SNOW/ICE COVER PRODUCT
            if val[0] == 'snow':

                assets = self._products['snow']['assets']

                # at some point can also include MYD10A1
                allsds = []
                missingassets = []
                for asset in assets:
                    try:
                        sds = self.assets[asset].datafiles()
                    except:
                        missingassets.append(asset)
                    else:
                        allsds.extend(sds)
                if len(missingassets) > 0:
                    print "Skipping missing data:", self.date, self.id, missingassets
                    continue

                # for i,sds in enumerate(allsds):
                #     print i,sds

                snowsds = [allsds[0], allsds[3]]
                img = gippy.GeoImage(snowsds) # read only

                # get the data values for both bands
                cover = img[0].Read()
                frac = img[1].Read()

                # check out frac
                wbad1 = np.where((frac==200)|(frac==201)|(frac==211)|(frac==250)|(frac==254)|(frac==255))
                wsurface1 = np.where((frac==225)|(frac==237)|(frac==239))
                wvalid1 = np.where((frac>=0)&(frac<=100))
                wmostly1 = np.where((frac>50)&(frac<=100))

                nbad1 = len(wbad1[0])
                nsurface1 = len(wsurface1[0])
                nvalid1 = len(wvalid1[0])
                assert nbad1 + nsurface1 + nvalid1 == frac.size, "frac contains invalid values"

                # check out cover
                wbad2 = np.where((cover==0)|(cover==1)|(cover==11)|(cover==50)|(cover==254)|(cover==255))
                wsurface2 = np.where((cover==25)|(cover==37)|(cover==39))
                wvalid2 = np.where((cover==100)|(cover==200))

                nbad2 = len(wbad2[0])
                nsurface2 = len(wsurface2[0])
                nvalid2 = len(wvalid2[0])
                assert nbad2 + nsurface2 + nvalid2 == cover.size, "cover contains invalid values"

                # assign output data here because frac will get overwritten
                coverout = np.zeros_like(cover, dtype=np.uint8)
                fracout = np.zeros_like(frac, dtype=np.uint8)

                fracout[wvalid1] = frac[wvalid1]
                fracout[wsurface1] = 0
                fracout[wbad1] = 127

                coverout[wvalid2] = 100
                coverout[wsurface2] = 0
                coverout[wbad2] = 127

                # now complete the checks and create metadata
                frac[wbad1] = 0
                frac[wsurface1] = 1
                frac[wvalid1] = 2
                cover[wbad2] = 0
                cover[wsurface2] = 1
                cover[wvalid2] = 2

                meta = {}
                meta['fracmissingcoverclear'] = np.sum((frac==0)&(cover==1))/float(cover.size)
                meta['fracmissingcoversnow'] = np.sum((frac==0)&(cover==2))/float(cover.size)
                meta['fracclearcovermissing'] = np.sum((frac==1)&(cover==0))/float(cover.size)
                meta['fracclearcoversnow'] = np.sum((frac==1)&(cover==2))/float(cover.size)
                meta['fracsnowcovermissing'] = np.sum((frac==2)&(cover==0))/float(cover.size)
                meta['fracsnowcoverclear'] = np.sum((frac==2)&(cover==1))/float(cover.size)
                meta['fracmostlybutclear'] = np.sum(cover[wmostly1]==1)/float(cover.size)
                meta['totsnowfrac'] = np.sum(fracout[wvalid1])
                meta['totsnowcover'] = np.sum(coverout[wvalid2])

                # create output gippy image
                imgout = gippy.GeoImage(outfname, img, gippy.GDT_Byte, 2)

                imgout.SetNoData(127)
                imgout.SetOffset(0.0)
                imgout.SetGain(1.0)

                imgout[0].Write(coverout)
                imgout[1].Write(fracout)

                imgout.SetColor('Snow Cover', 1)
                imgout.SetColor('Fractional Snow Cover', 2)

                # print imgout.BandNames()

                for k, v in meta.items():
                    imgout.SetMeta(k, str(v))

                # print imgout.GetMetaGroup('')

            # TEMPERATURE PRODUCT
            if val[0] == 'temp':

                assets = self._products['temp']['assets']

                # print "assets"
                # print assets

                # names = {
                #     '0': 'LST_Day_1km',
                #     '1': 'QC_Day',
                #     '2': 'Day_view_time',
                #     '3': 'Day_view_angl',
                #     '4': 'LST_Night_1km',
                #     '5': 'QC_Night',
                #     '6': 'Night_view_time',
                #     '7': 'Night_view_angl',
                #     '8': 'Emis_31',
                #     '9': 'Emis_32',
                #     '10': 'Clear_day_cov',
                #     '11': 'Clear_night_cov',
                # }

                allsds = []
                missingassets = []

                for asset in assets:
                    try:
                        sds = self.assets[asset].datafiles()
                    except:
                        missingassets.append(asset)
                    else:
                        allsds.extend(sds)

                if len(missingassets) > 0:
                    print "Skipping missing data:", self.date, self.id, missingassets
                    continue

                # for i,sds in enumerate(allsds):
                #     print i,sds

                tempsds = [allsds[0], allsds[4], allsds[12], allsds[16]]                    
                qcsds = [allsds[1], allsds[5], allsds[13], allsds[17]]
                hoursds = [allsds[2], allsds[6], allsds[14], allsds[18]]

                tempbands = GeoImage(tempsds) # read only
                qcbands = GeoImage(qcsds) # read only
                hourbands = GeoImage(hoursds) # read only

                metanames = {}

                # there are four temperature bands
                for iband in range(4):
 
                    tempnodatamask = tempbands[iband].NoDataMask()
                    tempnodatavalue = tempbands[iband].NoDataValue()
                    hournodatavalue = hourbands[iband].NoDataValue()

                    print "tempnodatavalue", tempnodatavalue
                    print "hournodatavalue", hournodatavalue

                    temp = tempbands[iband].Read()

                    print "temp.min, temp.max", temp.min(), temp.max()

                    print "len(nodatamask)", (tempnodatamask[tempnodatamask == 1]).size
                    print "len(temp == 0)", np.sum(temp == 0)
                    print "len(temp != 0)", np.sum(temp != 0)

                    qc = qcbands[iband].Read()

                    # first two bits are 10 or 11
                    newmaskbad = binmask(qc, 2)

                    # first two bits are 00 or 01
                    newmaskgood = ~binmask(qc, 2)

                    # first two bits are 00
                    newmaskbest = ~binmask(qc, 1) & ~binmask(qc, 2)

                    print "newmaskbest.dtype"
                    print newmaskbest.dtype

                    if iband == 0:
                        bestmask = newmaskbest.astype('uint16')
                    else:
                        bestmask += (math.pow(2, iband)*newmaskbest).astype('uint16')

                    numbad = np.sum(newmaskbad)
                    fracbad = np.sum(newmaskbad)/float(newmaskbad.size)
                    print "numbad, fracbad", numbad, fracbad
                    print "numgood (by difference)", qc.size - numbad

                    numgood = np.sum(newmaskgood)
                    fracgood = np.sum(newmaskgood)/float(newmaskgood.size)
                    print "numgood, fracgood", numgood, fracgood

                    numbest = np.sum(newmaskbest)
                    fracbest = np.sum(newmaskbest)/float(newmaskbest.size)
                    print "numbest, fracbest", numbest, fracbest

                    # overpass time
                    hour = hourbands[iband].Read()
                    hour = hour[hour != hournodatavalue]

                    try:
                        hourmin = hour.min()
                        hourmean = hour.mean()
                        hourmax = hour.max()
                        print "hour.min(), hour.mean(), hour.max()", hour.min(), hour.mean(), hour.max()
                    except:
                        hourmean = 0

                    basename = tempbands[iband].Basename()
                    platform = basename[:3]
                    platform = platformnames[platform]

                    dayornight = basename.split()[1]
                    dayornight = dayornight.replace('time', '')
                    assert dayornight in ('day', 'night')
                    metaname = "MEANOVERPASSTIME_%s_%s" % (platform, dayornight)
                    metaname = metaname.upper()

                    print "metaname", metaname
                    metanames[metaname] = str(hourmean)

                    print


                print "tempbands.GetMetaGroup('')"
                print tempbands.GetMetaGroup('')
                # print "tempbands.GetMeta()"
                # print tempbands.GetMeta()

                print "tempbands.Info()"
                print tempbands.Info()
                print "tempbands.BandNames()"
                print tempbands.BandNames()
                print "tempbands.Basename()"
                print tempbands.Basename()
                print "tempbands.DataType()"
                print tempbands.DataType()
                print "tempbands.NumBands()"
                print tempbands.NumBands()
                print "tempbands.NoDataMask()"
                print tempbands.NoDataMask()

                print
                print "tempbands[0].Description()"
                print tempbands[0].Description()
                print "tempbands[0].Basename()"
                print tempbands[0].Basename()
                print "tempbands[0].NoDataValue()"
                print tempbands[0].NoDataValue()
                print "tempbands[0].NoData()"
                print tempbands[0].NoData()
                print "tempbands[0].Info()"
                print tempbands[0].Info()
                print "tempbands[0].Units()"
                print tempbands[0].Units()
                print "tempbands[0].Gain()"
                print tempbands[0].Gain()
                print "tempbands[0].Offset()"
                print tempbands[0].Offset()
                print "tempbands[0].Filename()"
                print tempbands[0].Filename()
                print "tempbands[0].Format()"
                print tempbands[0].Format()

                print "writing", outfname
                imgout = gippy.GeoImage(outfname, tempbands, gippy.GDT_UInt16, 5)
                imgout.SetNoData(0)
                imgout.SetGain(0.02)

                for i in range(4):
                    imgout[i].Process(tempbands[i])

                    print "imgout[%d].Gain()" % i, imgout[i].Gain()
                    print "imgout[%d].Offset()" % i, imgout[i].Offset()

                    print imgout[i].Stats()

                print "bestmask.min(), bestmask.max()", bestmask.min(), bestmask.max()

                imgout[4].SetGain(1.0)
                imgout[4].Write(bestmask)

                print "imgout[4].Gain()", imgout[4].Gain()
                print "imgout[4].Offset()", imgout[4].Offset()

                print imgout[4].Stats()

                for k, v in metanames.items():
                    imgout.SetMeta(k, v)

                print imgout.GetMetaGroup('')


def main():
    DataInventory.main(ModisData)
