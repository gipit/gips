#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    Copyright (C) 2014 Bobby Braswell
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
from osgeo import gdal
import math
import numpy as np

import gippy
from gips.core import Repository, Asset, Data
from gips.inventory import DataInventory
from gips.utils import File2List, List2File, VerboseOut
import gips.settings as settings
from pdb import set_trace


def binmask(arr, bit):
    """ Return boolean array indicating which elements as binary have a 1 in
        a specified bit position. Input is Numpy array.
    """
    return arr & (1 << (bit - 1)) == (1 << (bit - 1))


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
        return "h%sv%s" % (h, v)

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
        'MOD09Q1': {
            'pattern': 'MOD09Q1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLT/MOD09Q1.005/',
            'startdate': datetime.date(2000, 2, 18),
            'latency': -7,
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
            'latency': -1
        },
        'MYD11A1': {
            'pattern': 'MYD11A1*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.005',
            'startdate': datetime.date(2002, 7, 8),
            'latency': -1
        },
        'MOD11A2': {
            'pattern': 'MOD11A2*hdf',
            'url': 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.005',
            'startdate': datetime.date(2000, 3, 5),
            'latency': -7
        },
    }

    # Should this be specified on a per asset basis? (in which case retrieve from asset)
    _defaultresolution = [926.625433138333392, -926.625433139166944]

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(ModisAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        self.asset = bname[0:7]
        self.tile = bname[17:23]
        year = bname[9:13]
        doy = bname[13:16]

        self.date = datetime.datetime.strptime(year + doy, "%Y%j").date()
        self.sensor = bname[:3]

    @classmethod
    def fetch(cls, asset, tile, date):
        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)

        if date.date() < cls._assets[asset]['startdate']:
            raise Exception("date earlier than %s" % cls._assets[asset]['startdate'])
        latest_date = datetime.datetime.now() - datetime.timedelta(cls._assets[asset]['latency'])
        if date > latest_date:
            raise Exception("date is too recent (after %s)" % latest_date)

        year, month, day = date.timetuple()[:3]
        mainurl = '%s/%s.%02d.%02d' % (cls._assets[asset]['url'], str(year), month, day)
        try:
            listing = urllib.urlopen(mainurl).readlines()
        except Exception:
            raise Exception("Unable to access %s" % mainurl)

        pattern = '(%s.A%s%s.%s.005.\d{13}.hdf)' % (asset, str(year), str(date.timetuple()[7]).zfill(3), tile)
        cpattern = re.compile(pattern)
        success = False

        for item in listing:
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                name = cpattern.findall(item)[0]
                url = ''.join([mainurl, '/', name])
                outpath = os.path.join(cls.Repository.spath(), name)
                try:
                    urllib.urlretrieve(url, outpath)
                except Exception:
                    raise Exception('Unable to retrieve %s from %s' % (name, url))
                else:
                    VerboseOut('Retrieved %s' % name, 2)
                    success = True

        if not success:
            raise Exception('Unable to find remote match for %s at %s' % (pattern, mainurl))


class ModisData(Data):
    """ A tile of data (all assets and products) """
    name = 'Modis'
    Asset = ModisAsset
    _pattern = '*.tif'
    _products = {
        # MCD Products
        'indices': {
            'description': 'Land indices based',
            'assets': ['MCD43A4', 'MCD43A2'],
            'group': "Nadir BRDF-Adjusted 16-day"
        },
        # Daily
        'snow': {
            'description': 'Snow and ice cover data',
            'assets': ['MOD10A1', 'MYD10A1'],
            'group': 'Terra/Aqua Daily'
        },
        'temp': {
            'description': 'Surface temperature data',
            'assets': ['MOD11A1', 'MYD11A1'],
            'group': 'Terra/Aqua Daily'
        },
        # Misc
        'ndvi8': {
            'description': 'Normalized Difference Vegetation Index: 250m',
            'assets': ['MOD09Q1'],
            'group': 'Terra 8-day'
        },
        'temp8': {
            'description': 'Surface temperature: 1km',
            'assets': ['MOD11A2'],
            'group': 'Terra 8-day'
        },

    }

    def process(self, products, **kwargs):
        """ Process all products """

        bname = '_'.join(self.basename.split('_')[:-1] + ['MOD'])

        for key, val in products.items():
            start = datetime.datetime.now()
            fname = os.path.join(self.path, bname + '_' + key)

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
                VerboseOut('There are no available assets (%s) on %s for tile %s'
                           % (str(missingassets), str(self.date), str(self.id), ), 3)
                continue

            meta = {}
            meta['AVAILABLE_ASSETS'] = ' '.join(availassets)

            # LAND VEGETATION INDICES PRODUCT
            if val[0] == "indices":
                VERSION = "1.0"
                meta['VERSION'] = VERSION

                # there should be 11 SDSs, 7 bands and 4 QC layers
                reflsds = [allsds[i] for i in range(7)]
                qcsds = [allsds[i] for i in range(7, 11)]

                refl = gippy.GeoImage(reflsds)
                qc = gippy.GeoImage(qcsds)

                missing = 32767

                redimg = refl[0].Read()
                nirimg = refl[1].Read()
                bluimg = refl[2].Read()
                grnimg = refl[3].Read()
                mirimg = refl[5].Read()
                swir2img = refl[6].Read()

                qcimg = qc[0].Read()
                qcimg = qcimg.astype(np.int16)
                qcimg[qcimg == 255] = missing

                redimg[redimg < 0.0] = 0.0
                nirimg[nirimg < 0.0] = 0.0
                bluimg[bluimg < 0.0] = 0.0
                grnimg[grnimg < 0.0] = 0.0
                mirimg[mirimg < 0.0] = 0.0
                swir2img[swir2img < 0.0] = 0.0

                ndvi = missing + np.zeros_like(redimg)
                wg = np.where((redimg != missing) & (nirimg != missing) & (redimg + nirimg != 0.0))
                ndvi[wg] = (nirimg[wg] - redimg[wg]) / (nirimg[wg] + redimg[wg])
                print "ndvi"
                print len(wg[0])
                print ndvi.min(), ndvi.max()
                print ndvi[wg].min(), ndvi[wg].max()

                lswi = missing + np.zeros_like(redimg)
                wg = np.where((nirimg != missing) & (mirimg != missing) & (nirimg + mirimg != 0.0))
                lswi[wg] = (nirimg[wg] - mirimg[wg]) / (nirimg[wg] + mirimg[wg])
                print "lswi"
                print len(wg[0])
                print lswi.min(), lswi.max()
                print lswi[wg].min(), lswi[wg].max()

                vari = missing + np.zeros_like(redimg)
                wg = np.where((grnimg != missing) & (redimg != missing) & (bluimg != missing) & (grnimg + redimg - bluimg != 0.0))
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

                satvi = missing + np.zeros_like(redimg)
                wg = np.where((redimg != missing)&(mirimg != missing)&(swir2img != missing)&(((mirimg + redimg + 0.5)*swir2img) != 0.0))
                satvi[wg] = (((mirimg[wg] - redimg[wg])/(mirimg[wg] + redimg[wg] + 0.5))*1.5) - (swir2img[wg] / 2.0)
                print "satvi"
                print len(wg[0])
                print satvi.min(), satvi.max()
                print satvi[wg].min(), satvi[wg].max()

                # create output gippy image
                imgout = gippy.GeoImage(fname, refl, gippy.GDT_Int16, 6)

                imgout.SetNoData(missing)
                imgout.SetOffset(0.0)
                imgout.SetGain(0.0001)

                imgout[0].Write(ndvi)
                imgout[1].Write(lswi)
                imgout[2].Write(vari)
                imgout[3].Write(brgt)
                imgout[4].Write(satvi)

                imgout[5].SetGain(1.0)
                imgout[5].Write(qcimg)

                imgout.SetBandName('NDVI', 1)
                imgout.SetBandName('LSWI', 2)
                imgout.SetBandName('VARI', 3)
                imgout.SetBandName('BRGT', 4)
                imgout.SetBandName('SATVI', 5)
                imgout.SetBandName('Best quality', 6)

                imgout.SetMeta(meta)

            # SNOW/ICE COVER PRODUCT
            if val[0] == "snow":
                VERSION = "1.0"
                meta['VERSION'] = VERSION

                if not missingassets:
                    availbands = [0, 1]
                    snowsds = [allsds[0], allsds[3], allsds[4], allsds[7]]
                elif missingassets[0] == 'MYD10A1':
                    availbands = [0]
                    snowsds = [allsds[0], allsds[3]]
                elif missingassets[0] == 'MOD10A1':
                    availbands = [1]
                    snowsds = [allsds[0], allsds[3]]
                else:
                    raise

                img = gippy.GeoImage(snowsds)

                # there are two snow bands
                for iband, band in enumerate(availbands):

                    # get the data values for both bands
                    cover = img[2 * iband].Read()
                    frac = img[2 * iband + 1].Read()

                    # check out frac
                    wbad1 = np.where((frac == 200) | (frac == 201) | (frac == 211) |
                                     (frac == 250) | (frac == 254) | (frac == 255))
                    wsurface1 = np.where((frac == 225) | (frac == 237) | (frac == 239))
                    wvalid1 = np.where((frac >= 0) & (frac <= 100))

                    nbad1 = len(wbad1[0])
                    nsurface1 = len(wsurface1[0])
                    nvalid1 = len(wvalid1[0])
                    assert nbad1 + nsurface1 + nvalid1 == frac.size, "frac contains invalid values"

                    # check out cover
                    wbad2 = np.where((cover == 0) | (cover == 1) | (cover == 11) |
                                     (cover == 50) | (cover == 254) | (cover == 255))
                    wsurface2 = np.where((cover == 25) | (cover == 37) | (cover == 39))
                    wvalid2 = np.where((cover == 100) | (cover == 200))

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

                    if len(availbands) == 2:
                        if iband == 0:
                            fracout1 = np.copy(fracout)
                            coverout1 = np.copy(coverout)
                        else:
                            # both the current and previous are valid
                            w = np.where((fracout != 127) & (fracout1 != 127))
                            fracout[w] = np.mean(np.array([fracout[w], fracout1[w]]), axis=0).astype('uint8')

                            # the current is not valid but previous is valid
                            w = np.where((fracout == 127) & (fracout1 != 127))
                            fracout[w] = fracout1[w]

                            # both the current and previous are valid
                            w = np.where((coverout != 127) & (coverout1 != 127))
                            coverout[w] = np.mean(np.array([coverout[w], coverout1[w]]), axis=0).astype('uint8')

                            # the current is not valid but previous is valid
                            w = np.where((coverout == 127) & (coverout1 != 127))
                            coverout[w] = coverout1[w]

                fracmissingcoverclear = np.sum((fracout == 127) & (coverout == 0))
                fracmissingcoversnow = np.sum((fracout == 127) & (coverout == 100))
                fracclearcovermissing = np.sum((fracout == 0) & (coverout == 127))
                fracclearcoversnow = np.sum((fracout == 0) & (coverout == 100))
                fracsnowcovermissing = np.sum((fracout > 0) & (fracout <= 100) & (coverout == 127))
                fracsnowcoverclear = np.sum((fracout > 0) & (fracout <= 100) & (coverout == 0))
                #fracmostlycoverclear = np.sum((fracout > 50) & (fracout <= 100) & (coverout == 0))
                totsnowfrac = int(0.01 * np.sum(fracout[fracout <= 100]))
                totsnowcover = int(0.01 * np.sum(coverout[coverout <= 100]))
                numvalidfrac = np.sum(fracout != 127)
                numvalidcover = np.sum(coverout != 127)

                if totsnowcover == 0 or totsnowfrac == 0:
                    print "no snow or ice: skipping", str(self.date), str(self.id), str(missingassets)

                meta['FRACMISSINGCOVERCLEAR'] = fracmissingcoverclear
                meta['FRACMISSINGCOVERSNOW'] = fracmissingcoversnow
                meta['FRACCLEARCOVERMISSING'] = fracclearcovermissing
                meta['FRACCLEARCOVERSNOW'] = fracclearcoversnow
                meta['FRACSNOWCOVERMISSING'] = fracsnowcovermissing
                meta['FRACSNOWCOVERCLEAR'] = fracsnowcoverclear
                meta['FRACMOSTLYCOVERCLEAR'] = np.sum((fracout > 50) & (fracout <= 100) & (coverout == 0))
                meta['TOTSNOWFRAC'] = totsnowfrac
                meta['TOTSNOWCOVER'] = totsnowcover
                meta['NUMVALIDFRAC'] = numvalidfrac
                meta['NUMVALIDCOVER'] = numvalidcover

                # create output gippy image
                imgout = gippy.GeoImage(fname, img, gippy.GDT_Byte, 2)
                imgout.SetNoData(127)
                imgout.SetOffset(0.0)
                imgout.SetGain(1.0)
                imgout.SetBandName('Snow Cover', 1)
                imgout.SetBandName('Fractional Snow Cover', 2)
                imgout.SetMeta(meta)

                imgout[0].Write(coverout)
                imgout[1].Write(fracout)

            # TEMPERATURE PRODUCT (DAILY)
            # TODO: # use a missing value that is not zero, e.g. 32767 because 0 is a valid value in the QC mask
            #  or use 1,2 in the QC mask instead of 0,1
            if val[0] == "temp":
                VERSION = "1.1"
                meta['VERSION'] = VERSION

                if not missingassets:
                    availbands = [0, 1, 2, 3]
                    tempsds = [allsds[0], allsds[4], allsds[12], allsds[16]]
                    qcsds = [allsds[1], allsds[5], allsds[13], allsds[17]]
                    hoursds = [allsds[2], allsds[6], allsds[14], allsds[18]]
                elif missingassets[0] == 'MYD11A1':
                    availbands = [0, 1]
                    tempsds = [allsds[0], allsds[4]]
                    qcsds = [allsds[1], allsds[5]]
                    hoursds = [allsds[2], allsds[6]]
                elif missingassets[0] == 'MOD11A1':
                    availbands = [2, 3]
                    tempsds = [allsds[0], allsds[4]]
                    qcsds = [allsds[1], allsds[5]]
                    hoursds = [allsds[2], allsds[6]]
                else:
                    raise

                tempbands = gippy.GeoImage(tempsds)
                qcbands = gippy.GeoImage(qcsds)
                hourbands = gippy.GeoImage(hoursds)

                imgout = gippy.GeoImage(fname, tempbands, gippy.GDT_UInt16, 5)
                imgout.SetNoData(65535)
                imgout.SetGain(0.02)

                imgout.SetBandName('Temperature Daytime Terra', 1)
                imgout.SetBandName('Temperature Nighttime Terra', 2)
                imgout.SetBandName('Temperature Daytime Aqua', 3)
                imgout.SetBandName('Temperature Nighttime Aqua', 4)
                imgout.SetBandName('Temperature Best Quality', 5)

                # there are four temperature bands
                for iband, band in enumerate(availbands):
                    # get meta name template info
                    basename = tempbands[iband].Basename()
                    print "basename", basename
                    platform = basename[:3]
                    platform = self.Asset._sensors[platform]

                    dayornight = basename.split()[2]
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
                        bestmask = np.zeros_like(qc, dtype='uint16')

                    bestmask += (math.pow(2, band) * newmaskbest).astype('uint16')

                    numbad = np.sum(newmaskbad)
                    #fracbad = np.sum(newmaskbad) / float(newmaskbad.size)

                    numgood = np.sum(newmaskgood)
                    #fracgood = np.sum(newmaskgood) / float(newmaskgood.size)
                    assert numgood == qc.size - numbad

                    numbest = np.sum(newmaskbest)
                    #fracbest = np.sum(newmaskbest) / float(newmaskbest.size)

                    metaname = "NUMBAD_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    # print "metaname", metaname
                    meta[metaname] = str(numbad)

                    metaname = "NUMGOOD_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    # print "metaname", metaname
                    meta[metaname] = str(numgood)

                    metaname = "NUMBEST_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    # print "metaname", metaname
                    meta[metaname] = str(numbest)

                    # overpass time
                    hournodatavalue = hourbands[iband].NoDataValue()
                    hour = hourbands[iband].Read()
                    hour = hour[hour != hournodatavalue]
                    try:
                        #hourmin = hour.min()
                        hourmean = hour.mean()
                        #hourmax = hour.max()
                        # print "hour.min(), hour.mean(), hour.max()", hour.min(), hour.mean(), hour.max()
                    except:
                        hourmean = 0

                    metaname = "MEANOVERPASSTIME_%s_%s" % (dayornight, platform)
                    metaname = metaname.upper()
                    meta[metaname] = str(hourmean)

                    tempbands[iband].Process(imgout[band])

                imgout[4].SetGain(1.0)
                imgout[4].Write(bestmask)

                imgout.SetMeta(meta)

            # NDVI (8-day) - Terra only
            if val[0] == "ndvi8":
                VERSION = "1.0"
                meta['VERSION'] = VERSION

                refl = gippy.GeoImage(allsds)
                refl.SetBandName("RED", 1)
                refl.SetBandName("NIR", 2)
                refl.SetNoData(32767)

                fouts = dict(gippy.Indices(refl, {'ndvi': fname}, meta))
                imgout = gippy.GeoImage(fouts['ndvi'])

            # TEMPERATURE PRODUCT (8-day) - Terray only
            if val[0] == "temp8":
                VERSION = "1.0"
                meta['VERSION'] = VERSION

                set_trace()
                os.symlink(allsds[0], fname + '.tif')
                set_trace()
                imgout = gippy.GeoImage(fname)

            # add product to inventory
            self.products[val[0]] = imgout.Filename()
            VerboseOut(' -> %s: processed in %s' % (os.path.basename(fname), datetime.now() - start), 1)


def main():
    DataInventory.main(ModisData)
