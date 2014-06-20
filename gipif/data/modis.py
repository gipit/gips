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
from gipif.core import Repository, Asset, Data
from gipif.inventory import DataInventory
from gipif.utils import File2List, List2File, VerboseOut
import gipif.settings as settings


def binmask(arr, bit):
    """ Return boolean array indicating which elements as binary have a 1 in
        a specified bit position. Input is Numpy array.
    """
    return arr & (1 << (bit-1)) == (1 << (bit-1))


class ModisRepository(Repository):
    repo = settings.REPOS['modis']
    _rootpath = repo.get('rootpath', Repository._rootpath)
    _tiles_vector = repo.get('tiles_vector', Repository._tiles_vector)
    _tile_attribute = repo.get('tile_attribute', Repository._tile_attribute)

    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")
        h = str(int(feature.GetField(fldindex_h))).zfill(2)
        v = str(int(feature.GetField(fldindex_v))).zfill(2)
        tile = "h%sv%s" % (h, v)
        return tile

    # @classmethod
    # def find_dates(cls, tile):
    #     """ Get list of dates available in repository for a tile """
    #     tdir = cls.path(tile=tile)
    #     if os.path.exists(tdir):
    #         return [datetime.strptime(os.path.basename(d), cls._datedir).date() for d in os.listdir(tdir)]
    #     else:
    #         return []


class ModisAsset(Asset):
    Repository = ModisRepository

    _sensors = {
        'MOD': {'description': 'Terra'},
        'MYD': {'description': 'Aqua'},
        'MCD': {'description': 'Combined'}
    }

    _asset_layers = {
        'MOD11': {
            '0': 'LST_Day_1km',
            '1': 'QC_Day',
            '2': 'Day_view_time',
            '3': 'Day_view_angl',
            '4': 'LST_Night_1km',
            '5': 'QC_Night',
            '6': 'Night_view_time',
            '7': 'Night_view_angl',
            '8': 'Emis_31',
            '9': 'Emis_32',
            '10': 'Clear_day_cov',
            '11': 'Clear_night_cov',
        }
    }

    _assets = {
        'MCD43A4': {
            'pattern': 'MCD43A4*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOTA/MCD43A4.005',
            'startdate': datetime.date(2000, 2, 18),
            'latency': -15
        },
        'MCD43A2': {
            'pattern': 'MCD43A2*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOTA/MCD43A2.005',
            'startdate': datetime.date(2000, 2, 18),
            'latency': -15
        },
        'MOD10A1': {
            'pattern': 'MOD10A1*hdf',
            'url': 'ftp://n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005',
            'startdate': datetime.date(2000, 2, 24),
            'latency': -3
        },
        'MYD10A1': {
            'pattern': 'MYD10A1*hdf',
            'url': 'ftp://n5eil01u.ecs.nsidc.org/SAN/MOSA/MYD10A1.005',
            'startdate': datetime.date(2002, 7, 4),
            'latency': -3
        },
        'MOD11A1': {
            'pattern': 'MOD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005',
            'startdate': datetime.date(2000, 3, 5),
            'latency': -1,
            'layers': _asset_layers['MOD11']
        },
        'MYD11A1': {
            'pattern': 'MYD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.005',
            'startdate': datetime.date(2002, 7, 8),
            'latency': -1,
            'layers': _asset_layers['MOD11']
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


    def datafiles(self):

        indexfile = self.filename + '.index'

        # TODO: get File2List to handle missing file and empty file in the same way
        try:
            datafiles = File2List(indexfile)
        except:
            datafiles = []

        if not datafiles:
            gdalfile = gdal.Open(self.filename)
            subdatasets = gdalfile.GetSubDatasets()
            datafiles = [s[0] for s in subdatasets]
            List2File(datafiles, indexfile)

        return datafiles

    @classmethod
    def _filepattern(cls, asset, tile, date):
        year, month, day = date.timetuple()[:3]
        doy = date.timetuple()[7]
        pattern = ''.join(['(', asset, '.A', str(year), str(doy).zfill(3), '.', tile, '.005.\d{13}.hdf)'])
        return pattern


    @classmethod
    def _remote_subdirs(cls, asset, tile, date):
        year, month, day = date.timetuple()[:3]
        httploc = cls._assets[asset]['url']
        mainurl = ''.join([httploc, '/', str(year), '.', '%02d' % month, '.', '%02d' % day])
        return mainurl


    @classmethod
    def fetch(cls, asset, tile, date):
        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)

        if date.date() < cls._assets[asset]['startdate']:
            print "date is too early"
            return 3

        if date > datetime.datetime.now() - datetime.timedelta(cls._assets[asset]['latency']):
            print "date is too recent"
            return 3

        pattern = cls._filepattern(asset, tile, date)
        mainurl = cls._remote_subdirs(asset, tile, date)

        VerboseOut('%s: mainurl %s, pattern %s' % (asset, mainurl, pattern), 4)

        try:
            VerboseOut('opening listing', 1)
            listing = urllib.urlopen(mainurl).readlines()
        except Exception, e:
            listing = None
            VerboseOut('unable to access %s' % mainurl, 1)
            # sys.exit()
            return 2

        cpattern = re.compile(pattern)
        name = None
        success = False

        outdir = cls.Repository.spath()
        # outdir = os.path.join(cls.Repository.spath(), asset)
        # if not os.path.exists(outdir):
        #     os.mkdir(outdir)

        for item in listing:
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                name = cpattern.findall(item)[0]
                VerboseOut('found %s in %s' % (name, item.strip()), 4)
                url = ''.join([mainurl, '/', name])
                VerboseOut('the url is %s' % url, 4)
                outpath = os.path.join(outdir, name)
                try:                        
                    urllib.urlretrieve(url, outpath)
                except Exception, e:
                    VerboseOut('unable to retrieve %s from %s' % (name, url), 4)
                else:
                    VerboseOut('retrieved %s' % name, 4)
                    success = True

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
            'assets': ['MOD11A1', 'MYD11A1']
        },
        'snow': {
            'description': 'Snow and ice cover data',
            'assets': ['MOD10A1', 'MYD10A1']
        },
        'indices': {
            'description': 'Land indices',
            'assets': ['MCD43A4', 'MCD43A2']
        }
    }
    
    def process(self, products, **kwargs):

        VerboseOut(kwargs, 3)
        VerboseOut(products, 3)

        platformnames = {'MOD':'Terra', 'MYD':'Aqua'}

        for key, val in products.items():

            productname = '_'.join(self.basename.split('_')[:-1] + ['MOD'])
            outfname = os.path.join(self.path, productname + '_' + key)        

            VerboseOut("self.path: %s" % self.path, 4)
            VerboseOut("productname: %s" % productname, 4)
            VerboseOut("outfname: %s" % outfname, 4)

            ########################
            # LAND VEGETATION INDICES PRODUCT
            if val[0] == "indices":
                VERSION = "1.0"
                assets = self._products['indices']['assets']

                allsds = []
                missingassets = []
                availassets = []
                assetids = []

                for asset in assets:
                    try:
                        sds = self.assets[asset].datafiles()
                    except Exception,e:
                        missingassets.append(asset)
                    else:
                        assetids.append(assets.index(asset))
                        availassets.append(asset)
                        allsds.extend(sds)

                if missingassets:
                    VerboseOut('There are missing assets: %s,%s,%s' % (str(self.date), str(self.id), str(missingassets)), 4)
                    continue
                    # raise Exception, "There are missing assets"
                    # print message, then continue

                print "assetids", assetids
                print "availassets", availassets
                print "missingassets", missingassets
                for i,sds in enumerate(allsds):
                    print "i, sds", i,sds

                # there should be 11 SDSs, 7 bands and 4 QC layers

                reflsds = [allsds[i] for i in range(7)]
                qcsds = [allsds[i] for i in range(7,11)]

                refl = gippy.GeoImage(reflsds) 
                qc = gippy.GeoImage(qcsds) 

                missing = 32767 

                redimg = refl[0].Read()
                nirimg = refl[1].Read()
                bluimg = refl[2].Read()
                grnimg = refl[3].Read()
                mirimg = refl[5].Read()

                qcimg = qc[0].Read()
                qcimg = qcimg.astype(np.int16)
                qcimg[qcimg==255] = missing

                redimg[redimg<0.0] = 0.0
                nirimg[nirimg<0.0] = 0.0
                bluimg[bluimg<0.0] = 0.0
                grnimg[grnimg<0.0] = 0.0
                mirimg[mirimg<0.0] = 0.0

                meta = {}
                meta['AVAILABLE_ASSETS'] = "['MCD43A4', 'MCD43A2']"
                meta['VERSION'] = VERSION

                ndvi = missing + np.zeros_like(redimg)
                wg = np.where((redimg != missing)&(nirimg != missing)&(redimg+nirimg != 0.0))
                ndvi[wg] = (nirimg[wg] - redimg[wg])/(nirimg[wg] + redimg[wg])
                print "ndvi" 
                print len(wg[0])
                print ndvi.min(), ndvi.max()
                print ndvi[wg].min(), ndvi[wg].max()

                lswi = missing + np.zeros_like(redimg)
                wg = np.where((nirimg != missing)&(mirimg != missing)&(nirimg+mirimg != 0.0))
                lswi[wg] = (nirimg[wg] - mirimg[wg])/(nirimg[wg] + mirimg[wg])
                print "lswi"
                print len(wg[0])
                print lswi.min(), lswi.max()
                print lswi[wg].min(), lswi[wg].max()

                vari = missing + np.zeros_like(redimg)
                wg = np.where((grnimg != missing)&(redimg != missing)&(bluimg != missing)&(grnimg+redimg-bluimg != 0.0))
                vari[wg] = (grnimg[wg] - redimg[wg]) / (grnimg[wg] + redimg[wg] - bluimg[wg])
                print "vari"
                print len(wg[0])
                print vari.min(), vari.max()
                print vari[wg].min(), vari[wg].max()

                brgt = missing + np.zeros_like(redimg)
                wg = np.where((nirimg != missing)&(redimg != missing)&(bluimg != missing)&(grnimg != missing))
                brgt[wg] = 0.3*bluimg[wg] + 0.3*redimg[wg] + 0.1*nirimg[wg] + 0.3*grnimg[wg]
                print "brgt"
                print len(wg[0])
                print brgt.min(), brgt.max()
                print brgt[wg].min(), brgt[wg].max()


                # create output gippy image
                imgout = gippy.GeoImage(outfname, refl, gippy.GDT_Int16, 5)

                imgout.SetNoData(missing)
                imgout.SetOffset(0.0)
                imgout.SetGain(0.0001)

                imgout[0].Write(ndvi)
                imgout[1].Write(lswi)
                imgout[2].Write(vari)
                imgout[3].Write(brgt)

                imgout[4].SetGain(1.0)
                imgout[4].Write(qcimg)

                imgout.SetColor('NDVI', 1)
                imgout.SetColor('LSWI', 2)
                imgout.SetColor('VARI', 3)
                imgout.SetColor('BRGT', 4)
                imgout.SetColor('Best quality', 5)

                for k, v in meta.items():
                    imgout.SetMeta(k, str(v))


            ########################
            # SNOW/ICE COVER PRODUCT
            ########################

            if val[0] == "snow":
                VERSION = "1.0"
                assets = self._products['snow']['assets']

                allsds = []
                missingassets = []
                availassets = []
                assetids = []

                for asset in assets:
                    try:
                        sds = self.assets[asset].datafiles()
                    except Exception,e:
                        missingassets.append(asset)
                    else:
                        assetids.append(assets.index(asset))
                        availassets.append(asset)
                        allsds.extend(sds)

                # print "assetids", assetids
                # print "availassets", availassets
                # print "missingassets", missingassets
                # for i,sds in enumerate(allsds):
                #     print "i, sds", i,sds

                if not availassets:
                    VerboseOut('Both assets are missing: %s,%s,%s' % (str(self.date), str(self.id), str(missingassets)), 4)
                    continue
                    # raise Exception, "Both assets are missing"


                if not missingassets:
                    availbands = [0,1]
                    snowsds = [allsds[0], allsds[3], allsds[4], allsds[7]]
                elif missingassets[0] == 'MYD10A1':
                    availbands = [0]
                    snowsds = [allsds[0], allsds[3]]
                elif missingassets[0] == 'MOD10A1':
                    availbands = [1]
                    snowsds = [allsds[0], allsds[3]]
                else:
                    raise

                img = gippy.GeoImage(snowsds) # read only

                meta = {}
                meta['AVAILABLE_ASSETS'] = str(availassets)
                meta['VERSION'] = VERSION

                # there are two snow bands
                for iband, band in enumerate(availbands):

                    # get the data values for both bands
                    cover = img[2*iband].Read()
                    frac = img[2*iband+1].Read()

                    # check out frac
                    wbad1 = np.where((frac==200)|(frac==201)|(frac==211)|(frac==250)|(frac==254)|(frac==255))
                    wsurface1 = np.where((frac==225)|(frac==237)|(frac==239))
                    wvalid1 = np.where((frac>=0)&(frac<=100))

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

                    # assign output data here
                    coverout = np.zeros_like(cover, dtype=np.uint8)
                    fracout = np.zeros_like(frac, dtype=np.uint8)

                    fracout[wvalid1] = frac[wvalid1]
                    fracout[wsurface1] = 0
                    fracout[wbad1] = 127
                    coverout[wvalid2] = 100
                    coverout[wsurface2] = 0
                    coverout[wbad2] = 127

                    if len(availbands)==2:
                        if iband==0:
                            fracout1 = np.copy(fracout)
                            coverout1 = np.copy(coverout)
                        else:
                            # both the current and previous are valid
                            w = np.where((fracout != 127)&(fracout1 != 127))
                            fracout[w] = np.mean(np.array([fracout[w], fracout1[w]]), axis=0).astype('uint8')
                            # print "averaged this many values", len(w[0])

                            # the current is not valid but previous is valid
                            w = np.where((fracout == 127)&(fracout1 != 127))
                            fracout[w] = fracout1[w]
                            # print "replaced this many values", len(w[0])

                            # both the current and previous are valid
                            w = np.where((coverout != 127)&(coverout1 != 127))
                            coverout[w] = np.mean(np.array([coverout[w], coverout1[w]]), axis=0).astype('uint8')
                            # print "averaged this many values", len(w[0])

                            # the current is not valid but previous is valid
                            w = np.where((coverout == 127)&(coverout1 != 127))
                            coverout[w] = coverout1[w]
                            # print "replaced this many values", len(w[0])


                fracmissingcoverclear = np.sum((fracout==127)&(coverout==0))
                fracmissingcoversnow = np.sum((fracout==127)&(coverout==100))
                fracclearcovermissing = np.sum((fracout==0)&(coverout==127))
                fracclearcoversnow = np.sum((fracout==0)&(coverout==100))
                fracsnowcovermissing = np.sum((fracout>0)&(fracout<=100)&(coverout==127))
                fracsnowcoverclear = np.sum((fracout>0)&(fracout<=100)&(coverout==0))
                fracmostlycoverclear = np.sum((fracout>50)&(fracout<=100)&(coverout==0))
                totsnowfrac = int(0.01*np.sum(fracout[fracout<=100]))
                totsnowcover = int(0.01*np.sum(coverout[coverout<=100]))
                numvalidfrac = np.sum(fracout!=127)
                numvalidcover = np.sum(coverout!=127)

                if totsnowcover == 0 or totsnowfrac == 0:
                    print "no snow or ice: skipping", str(self.date), str(self.id), str(missingassets)


                meta['FRACMISSINGCOVERCLEAR'] = fracmissingcoverclear
                meta['FRACMISSINGCOVERSNOW'] = fracmissingcoversnow
                meta['FRACCLEARCOVERMISSING'] = fracclearcovermissing
                meta['FRACCLEARCOVERSNOW'] = fracclearcoversnow
                meta['FRACSNOWCOVERMISSING'] = fracsnowcovermissing
                meta['FRACSNOWCOVERCLEAR'] = fracsnowcoverclear
                meta['FRACMOSTLYCOVERCLEAR'] = np.sum((fracout>50)&(fracout<=100)&(coverout==0))
                meta['TOTSNOWFRAC'] = totsnowfrac
                meta['TOTSNOWCOVER'] = totsnowcover
                meta['NUMVALIDFRAC'] = numvalidfrac
                meta['NUMVALIDCOVER'] = numvalidcover

                # create output gippy image
                imgout = gippy.GeoImage(outfname, img, gippy.GDT_Byte, 2)

                imgout.SetNoData(127)
                imgout.SetOffset(0.0)
                imgout.SetGain(1.0)

                imgout[0].Write(coverout)
                imgout[1].Write(fracout)

                imgout.SetColor('Snow Cover', 1)
                imgout.SetColor('Fractional Snow Cover', 2)

                for k, v in meta.items():
                    imgout.SetMeta(k, str(v))


            #####################
            # TEMPERATURE PRODUCT

            # TODO:
            # use a missing value that is not zero, e.g. 32767
            # because 0 is a valid value in the QC mask
            # OR:
            # use 1,2 in the QC mask instead of 0,1

            if val[0] == "temp":
                VERSION = "1.0"
                assets = self._products['temp']['assets']

                allsds = []
                missingassets = []
                availassets = []
                assetids = []

                for asset in assets:
                    try:
                        sds = self.assets[asset].datafiles()
                    except Exception,e:
                        missingassets.append(asset)
                    else:
                        assetids.append(assets.index(asset))
                        availassets.append(asset)
                        allsds.extend(sds)

                # print "assetids", assetids
                # print "availassets", availassets
                # print "missingassets", missingassets
                # for i,sds in enumerate(allsds):
                #     print "i, sds", i,sds

                if not availassets:
                    VerboseOut('Both assets are missing: %s,%s,%s' % (str(self.date), str(self.id), str(missingassets)), 4)
                    continue

                if not missingassets:
                    availbands = [0,1,2,3]
                    tempsds = [allsds[0], allsds[4], allsds[12], allsds[16]]                    
                    qcsds = [allsds[1], allsds[5], allsds[13], allsds[17]]
                    hoursds = [allsds[2], allsds[6], allsds[14], allsds[18]]
                elif missingassets[0] == 'MYD11A1':
                    availbands = [0,1]
                    tempsds = [allsds[0], allsds[4]]                    
                    qcsds = [allsds[1], allsds[5]]
                    hoursds = [allsds[2], allsds[6]]
                elif missingassets[0] == 'MOD11A1':
                    availbands = [2,3]
                    tempsds = [allsds[0], allsds[4]]                    
                    qcsds = [allsds[1], allsds[5]]
                    hoursds = [allsds[2], allsds[6]]
                else:
                    raise

                tempbands = gippy.GeoImage(tempsds) # read only
                qcbands = gippy.GeoImage(qcsds) # read only
                hourbands = gippy.GeoImage(hoursds) # read only

                metanames = {}
                metanames['AVAILABLE_ASSETS'] = str(availassets)
                metanames['VERSION'] = VERSION

                # there are four temperature bands
                for iband, band in enumerate(availbands):

                    # get meta name template info
                    basename = tempbands[iband].Basename()
                    platform = basename[:3]
                    platform = platformnames[platform]
                    dayornight = basename.split()[1]
                    dayornight = dayornight.replace('time', '')
                    assert dayornight in ('day', 'night')

                    qc = qcbands[iband].Read()

                    # first two bits are 10 or 11
                    newmaskbad = binmask(qc, 2)
                    # first two bits are 00 or 01
                    newmaskgood = ~binmask(qc, 2)
                    # first two bits are 00
                    newmaskbest = ~binmask(qc, 1) & ~binmask(qc, 2)

                    if iband == 0:
                        bestmask = newmaskbest.astype('uint16')
                    else:
                        bestmask += (math.pow(2, iband)*newmaskbest).astype('uint16')

                    numbad = np.sum(newmaskbad)
                    fracbad = np.sum(newmaskbad)/float(newmaskbad.size)

                    numgood = np.sum(newmaskgood)
                    fracgood = np.sum(newmaskgood)/float(newmaskgood.size)
                    assert numgood == qc.size - numbad

                    numbest = np.sum(newmaskbest)
                    fracbest = np.sum(newmaskbest)/float(newmaskbest.size)

                    metaname = "NUMBAD_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    # print "metaname", metaname
                    metanames[metaname] = str(numbad)

                    metaname = "NUMGOOD_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    # print "metaname", metaname
                    metanames[metaname] = str(numgood)

                    metaname = "NUMBEST_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    # print "metaname", metaname
                    metanames[metaname] = str(numbest)

                    # overpass time
                    hournodatavalue = hourbands[iband].NoDataValue()
                    hour = hourbands[iband].Read()
                    hour = hour[hour != hournodatavalue]
                    try:
                        hourmin = hour.min()
                        hourmean = hour.mean()
                        hourmax = hour.max()
                        # print "hour.min(), hour.mean(), hour.max()", hour.min(), hour.mean(), hour.max()
                    except:
                        hourmean = 0

                    metaname = "MEANOVERPASSTIME_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    metanames[metaname] = str(hourmean)


                VerboseOut('writing %s' % outfname, 4)
                imgout = gippy.GeoImage(outfname, tempbands, gippy.GDT_UInt16, 5)
                imgout.SetNoData(0)
                imgout.SetGain(0.02)

                imgout.SetColor('Temperature Daytime Terra', 1)
                imgout.SetColor('Temperature Nighttime Terra', 2)
                imgout.SetColor('Temperature Daytime Aqua', 3)
                imgout.SetColor('Temperature Nighttime Aqua', 4)
                imgout.SetColor('Temperature Best Quality', 5)

                for iband, band in enumerate(availbands):
                    tempbands[iband].Process(imgout[band])
                    # print "imgout[%d].Gain()" % band, imgout[band].Gain()
                    # print "imgout[%d].Offset()" % band, imgout[band].Offset()
                    # print imgout[band].Stats()

                imgout[4].SetGain(1.0)
                imgout[4].Write(bestmask)
                # print "imgout[4].Gain()", imgout[4].Gain()
                # print "imgout[4].Offset()", imgout[4].Offset()
                # print imgout[4].Stats()

                for k, v in metanames.items():
                    imgout.SetMeta(k, v)

            # add product to inventory
            self.products[val[0]] = imgout.Filename()

def main():
    DataInventory.main(ModisData)
