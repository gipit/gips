#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  matt.a.hanson@gmail.com 
#
#    Copyright (C) 2014 Applied Geosolutions
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

import os
import re
import datetime
import urllib
import urllib2

import math
import numpy as np

import gippy
# TODO: Use this:
# from gippy.algorithms import Indices
from gips.data.core import Repository, Asset, Data
from gips.utils import VerboseOut


PROJ = """PROJCS["WELD_CONUS",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_84",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["longitude_of_center",-96.0],PARAMETER["Standard_Parallel_1",29.5],PARAMETER["Standard_Parallel_2",45.5],PARAMETER["latitude_of_center",23.0],UNIT["Meter",1.0]]"""


def binmask(arr, bit):
    """ Return boolean array indicating which elements as binary have a 1 in
        a specified bit position. Input is Numpy array.
    """
    return arr & (1 << (bit - 1)) == (1 << (bit - 1))


class weldRepository(Repository):
    name = 'WELD'
    description = 'WELD Landsat'

    @classmethod
    def feature2tile(cls, feature):
        """ convert tile field attributes to tile identifier """
        fldindex_h = feature.GetFieldIndex("h")
        fldindex_v = feature.GetFieldIndex("v")
        h = str(int(feature.GetField(fldindex_h))).zfill(2)
        v = str(int(feature.GetField(fldindex_v))).zfill(2)
        return "h%sv%s" % (h, v)


class weldAsset(Asset):
    Repository = weldRepository

    _sensors = {
        'WELD': {'description': 'WELD Landsat'},
    }

    _assets = {
        'CONUS': {
            'pattern': 'CONUS.week??.????.h??v??.doy???to???.v1.5.hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/WELD/WELDUSWK.001',
            'startdate': datetime.date(2003, 1, 1),
            'latency': None
        },
    }

    _defaultresolution = [30.0, 30.0]

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(weldAsset, self).__init__(filename)
        bname = os.path.basename(filename)
        parts = bname.split('.')
        self.asset = parts[0]
        self.tile = parts[3]
        year = parts[2]
        doy1 = parts[4][3:6]
        doy2 = parts[4][8:11]
        week = parts[1][4:6]
        self.sensor = "WELD"
        doy = str(datetime.timedelta(7*(int(week) - 1) + 1).days).zfill(3)
        self.date = datetime.datetime.strptime(year + doy, "%Y%j").date()

    @classmethod
    def fetch(cls, asset, tile, date):
        year, month, day = date.timetuple()[:3]
        mainurl = '%s/%s.%02d.%02d' % (cls._assets[asset]['url'], str(year), month, day)
        try:
            listing = urllib.urlopen(mainurl).readlines()
        except Exception:
            raise Exception("Unable to access %s" % mainurl)
        pattern = '(%s.week\d{2}.%s.%s.doy\d{3}to\d{3}.v1.5.hdf)' % (asset, str(year), tile)
        cpattern = re.compile(pattern)
        success = False
        for item in listing:
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                name = cpattern.findall(item)[0]
                url = ''.join([mainurl, '/', name])
                outpath = os.path.join(cls.Repository.path('stage'), name)
                if os.path.exists(outpath):
                    continue
                try:
                    connection = urllib2.urlopen(url)
                    output = open(outpath, 'wb')
                    output.write(connection.read())
                    output.close()
                except Exception:
                    # TODO - implement pre-check to only attempt on valid dates
                    # then uncomment this
                    #raise Exception('Unable to retrieve %s from %s' % (name, url))
                    pass
                else:
                    VerboseOut('Retrieved %s' % name, 2)
                    success = True
        if not success:
            # TODO - implement pre-check to only attempt on valid dates then uncomment below
            #raise Exception('Unable to find remote match for %s at %s' % (pattern, mainurl))
            # VerboseOut('Unable to find remote match for %s at %s' % (pattern, mainurl), 4)
            pass


class weldData(Data):
    """ A tile of data (all assets and products) """
    name = 'WELD'
    version = '0.1.0'
    Asset = weldAsset
    _pattern = '*.tif'
    _productgroups = {
        "indices": ['ndsi', 'ndvi', 'brgt'],
    }
    _products = {
        'ndsi': {
            'description': 'Snow index',
            'assets': ['CONUS'],
        },
        'ndvi': {
            'description': 'Vegetation index',
            'assets': ['CONUS'],
        },
        'brgt': {
            'description': 'Brightness index',
            'assets': ['CONUS'],
        },
        'snow': {
            'description': 'Snow cover',
            'assets': ['CONUS'],
        },
    }

    def process(self, *args, **kwargs):
        """ Process all products """
        products = super(weldData, self).process(*args, **kwargs)
        if len(products) == 0:
            return
        bname = os.path.join(self.path, self.basename)
        for key, val in products.requested.items():
            start = datetime.datetime.now()
            sensor = 'WELD'
            # Check for asset availability
            assets = self._products[val[0]]['assets']
            missingassets = []
            availassets = []
            allsds = []
            for asset in assets:
                try:
                    sds = self.assets[asset].datafiles()
                except Exception:
                    missingassets.append(asset)
                else:
                    availassets.append(asset)
                    allsds.extend(sds)
            if not availassets:
                # some products aren't available for every day but this is trying every day
                VerboseOut('There are no available assets (%s) on %s for tile %s'
                           % (str(missingassets), str(self.date), str(self.id), ), 5)
                continue
            meta = self.meta_dict()
            meta['AVAILABLE_ASSETS'] = ' '.join(availassets)

            # SNOW ICE INDEX PRODUCT
            if val[0] == "ndsi":
                VERSION = "1.0"
                meta['VERSION'] = VERSION
                sensor = "WELD"
                fname = "%s_%s_%s" % (bname, sensor, key)
                refl = gippy.GeoImage(allsds)
                # band 2
                grnimg = refl[1].Read()
                missing = refl[1].NoDataValue()
                # band 4
                nirimg = refl[3].Read()
                assert refl[3].NoDataValue() == missing
                # band 5
                swrimg = refl[4].Read()
                assert refl[4].NoDataValue() == missing
                cldimg = refl[11].Read()
                # accaimg = refl[13].Read()
                ndsi = missing + np.zeros_like(grnimg)
                wg = np.where((grnimg != missing) & (swrimg != missing) & (grnimg + swrimg != 0.0) & (cldimg == 0))
                ng = len(wg[0])
                print "ng", ng
                if ng == 0:
                    continue
                ndsi[wg] = (grnimg[wg] - swrimg[wg]) / (grnimg[wg] + swrimg[wg])
                print ndsi.min(), ndsi.max()
                print ndsi[wg].min(), ndsi[wg].max()
                print "writing", fname
                imgout = gippy.GeoImage(fname, refl, gippy.GDT_Float32, 1)
                imgout.SetNoData(float(missing))
                imgout.SetOffset(0.0)
                imgout.SetGain(1.0)
                imgout.SetProjection(PROJ)
                imgout[0].Write(ndsi)
                imgout.SetBandName('NDSI', 1)

            # SNOW ICE COVER PRODUCT
            if val[0] == "snow":
                VERSION = "1.0"
                meta['VERSION'] = VERSION
                sensor = "WELD"
                fname = "%s_%s_%s" % (bname, sensor, key)
                refl = gippy.GeoImage(allsds)
                # band 2
                grnimg = refl[1].Read()
                missing = refl[1].NoDataValue()
                # band 4
                nirimg = refl[3].Read()
                assert refl[3].NoDataValue() == missing
                # band 5
                swrimg = refl[4].Read()
                assert refl[4].NoDataValue() == missing
                cldimg = refl[11].Read()
                accaimg = refl[13].Read()
                snow = 127 + np.zeros_like(grnimg)
                ndsi = missing + np.zeros_like(grnimg)
                wg = np.where((grnimg != missing) & (swrimg != missing) & (grnimg + swrimg != 0.0) & (cldimg == 0))
                ng = len(wg[0])
                print "ng", ng
                if ng == 0:
                    continue
                ndsi[wg] = (grnimg[wg] - swrimg[wg]) / (grnimg[wg] + swrimg[wg])
                ws = np.where((ndsi != missing) & (ndsi > 0.4) & (nirimg > 0.11) & (swrimg > 0.1))
                wc = np.where((ndsi != missing) & (ndsi > 0.4) & (nirimg <= 0.11) & (swrimg <= 0.1))
                ns = len(ws[0])
                nc = len(wc[0])
                print ng, ns, nc
                if (ns > 0):
                    snow[ws] = 1
                if (nc > 0):
                    snow[wc] = 0
                print "writing", fname
                imgout = gippy.GeoImage(fname, refl, gippy.GDT_Byte, 1)
                imgout.SetNoData(127)
                imgout.SetOffset(0.)
                imgout.SetGain(1.)
                imgout.SetProjection(PROJ)
                imgout[0].Write(snow)
                imgout.SetBandName('SNOW', 1)

            # VEGETATION INDEX PRODUCT
            if val[0] == "ndvi":
                VERSION = "1.0"
                meta['VERSION'] = VERSION
                sensor = "WELD"
                fname = "%s_%s_%s" % (bname, sensor, key)
                refl = gippy.GeoImage(allsds)
                # band 3
                redimg = refl[2].Read()
                missing = refl[2].NoDataValue()
                # band 4
                nirimg = refl[3].Read()
                assert refl[3].NoDataValue() == missing
                cldimg = refl[11].Read()
                accaimg = refl[13].Read()
                ndvi = missing + np.zeros_like(redimg)
                wg = np.where((redimg != missing) & (nirimg != missing) & (redimg + nirimg != 0.0) & (cldimg == 0))
                ng = len(wg[0])
                print "ng", ng
                if ng == 0:
                    continue
                ndvi[wg] = (nirimg[wg] - redimg[wg]) / (nirimg[wg] + redimg[wg])
                print ndvi.min(), ndvi.max()
                print ndvi[wg].min(), ndvi[wg].max()
                print "writing", fname
                imgout = gippy.GeoImage(fname, refl, gippy.GDT_Float32, 1)
                imgout.SetNoData(float(missing))
                imgout.SetOffset(0.0)
                imgout.SetGain(1.0)
                imgout.SetProjection(PROJ)
                imgout[0].Write(ndvi)
                imgout.SetBandName('NDVI', 1)

            # BRIGHTNESS PRODUCT
            if val[0] == "brgt":
                VERSION = "1.0"
                meta['VERSION'] = VERSION
                sensor = "WELD"
                fname = "%s_%s_%s" % (bname, sensor, key)
                refl = gippy.GeoImage(allsds)
                # band 2
                grnimg = refl[1].Read()
                missing = refl[1].NoDataValue()
                # band 3
                redimg = refl[2].Read()
                assert refl[2].NoDataValue() == missing
                # band 4
                nirimg = refl[3].Read()
                assert refl[3].NoDataValue() == missing
                cldimg = refl[11].Read()
                brgt = missing + np.zeros_like(redimg)
                wg = np.where((grnimg != missing) & (redimg != missing) & (nirimg != missing) & (cldimg == 0))
                ng = len(wg[0])
                print "ng", ng
                if ng == 0:
                    continue
                brgt[wg] = (grnimg[wg] + redimg[wg] + nirimg[wg])/3.
                print brgt.min(), brgt.max()
                print brgt[wg].min(), brgt[wg].max()
                print "writing", fname
                imgout = gippy.GeoImage(fname, refl, gippy.GDT_Float32, 1)
                imgout.SetNoData(float(missing))
                imgout.SetOffset(0.0)
                imgout.SetGain(1.0)
                imgout.SetProjection(PROJ)
                imgout[0].Write(brgt)
                imgout.SetBandName('BRGT', 1)

            # set metadata
            meta = {k: str(v) for k, v in meta.iteritems()}
            imgout.SetMeta(meta)

            # add product to inventory
            self.AddFile(sensor, key, imgout.Filename())
            VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname), datetime.datetime.now() - start), 1)
