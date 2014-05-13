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

import os
import datetime
import gdal
import ftplib
import numpy
import glob

import gippy
from gipif.core import Repository, Asset, Data
from gipif.inventory import DataInventory
from gipif.utils import File2List, List2File, VerboseOut
import gipif.settings as settings


class AODRepository(Repository):
    repo = settings.REPOS['AOD']
    _rootpath = repo.get('rootpath', Repository._rootpath)
    _tiles_vector = repo.get('tiles_vector', Repository._tiles_vector)
    _tile_attribute = repo.get('tile_attribute', Repository._tile_attribute)

    _datedir = '%Y%j'

    @classmethod
    def path(cls, tile='', date=''):
        path = os.path.join(cls._rootpath, cls._tilesdir)
        if date != '':
            path = os.path.join(path, str(date.strftime('%Y')), str(date.strftime('%j')))
        return path

    @classmethod
    def find_tiles(cls):
        return ['']

    @classmethod
    def find_dates(cls, tile=''):
        """ Get list of dates available in repository for a tile """
        #tdir = cls.path()
        #if os.path.exists(tdir):
        dates = []
        for year in os.listdir(cls.path()):
            days = os.listdir(os.path.join(cls.path(), year))
            for day in days:
                dates.append(datetime.datetime.strptime(year+day, '%Y%j').date())
        return dates

    @classmethod
    def vector2tiles(cls, *args, **kwargs):
        """ There are no tiles """
        return {'': (1, 1)}


class AODAsset(Asset):
    Repository = AODRepository

    # ???? Not specific to MODIS
    _sensors = {
        'MOD': {'description': 'MODIS Terra'},
        #'MYD': {'description': 'MODIS Aqua'},
    }
    _assets = {
        'MOD08': {
            'pattern': 'MOD08_D3*hdf',
            'url': 'ladsweb.nascom.nasa.gov/allData/51/MOD08_D3'
        },
        #'MYD08': {
        #    'pattern': 'MYD08_D3*hdf',
        #    'url': 'ladsweb.nascom.nasa.gov/allData/51/MYD08_D3'
        #}
    }

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(AODAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        self.asset = bname[0:5]
        self.tile = ''
        year = bname[10:14]
        doy = bname[14:17]
        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()
        self.sensor = bname[:3]
        #datafiles = self.datafiles()
        prefix = 'HDF4_EOS:EOS_GRID:"'
        self.products = {'aod': prefix + filename + '":mod08:Optical_Depth_Land_And_Ocean_Mean'}

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
    def archive(cls, path='.', recursive=False, keep=False):
        assets = super(AODAsset, cls).archive(path, recursive, keep)
        # this creates new LTA files every archiving
        #dates = [a.date for a in assets]
        #for date in set(dates):
        #    AODData.process_aerolta_daily(date.strftime('%j'))


class AODData(Data):
    name = 'Globally Gridded Atmospheric Data'
    Asset = AODAsset

    _products = {
        'aod': {
            'description': 'Aerosol Optical Depth',
            # the list of asset types associated with this product
            'assets': ['MOD08'],  # , 'MYD08'],
        },
        'ltad': {
            'description': 'Average daily AOD',
            # the list of asset types associated with this product
            'assets': ['MOD08'],  # , 'MYD08'],
            'composite': True,
        },
        'lta': {
            'description': 'Average AOD',
            'assets': ['MOD08'],
            'composite': True
        }
    }

    #def process(self, products):
    #    start = datetime.datetime.now()
        #bname = os.path.basename(self.assets[''].filename)
    #    for product in products:
            #if product == 'aerolta':
            #    self.process_aerolta()
    #        VerboseOut(' -> %s: processed %s in %s' % (fout, product, datetime.datetime.now()-start))

    @classmethod
    def process_composites(cls, inventory, products, **kwargs):
        for product in products:
            cpath = cls.Asset.Repository.cpath('ltad')
            path = os.path.join(cpath, 'ltad')
            # Calculate AOT long-term multi-year averages (lta) for given day
            if product == 'ltad':
                for day in range(inventory.start_day, inventory.end_day+1):
                    dates = [d for d in inventory.dates if int(d.strftime('%j')) == day]
                    filenames = [inventory[d].tiles[''].products['aod'] for d in dates]
                    fout = path + '%s.tif' % str(day).zfill(3)
                    imgout = cls.process_mean(filenames, fout)
            # Calculate single average per pixel (all days and years)
            if product == 'lta':
                filenames = glob.glob(path+'*.tif')
                if len(filenames) > 0:
                    fout = os.path.join(cpath, 'lta.tif')
                    imgout = cls.process_mean(filenames, fout)
                else:
                    raise Exception('No daily LTA files exist!')
            imgout = None

    @classmethod
    def process_mean(cls, filenames, fout):
        """ Calculates mean of all filenames, and per pixel variances """
        start = datetime.datetime.now()
        if len(filenames) > 0:
            img = gippy.GeoImage(filenames)
            imgout = gippy.GeoImage(fout, img, gippy.GDT_Float32, 2)
            imgout.SetNoData(-32768)
            img.Mean(imgout[0])
            meanimg = imgout[0].Read()
            for band in range(0, img.NumBands()):
                data = img[band].Read()
                mask = img[band].DataMask()
                var = numpy.multiply(numpy.power(data-meanimg, 2), mask)
                if band == 0:
                    totalvar = var
                    counts = mask
                else:
                    totalvar = totalvar + var
                    counts = counts + mask
            inds = numpy.where(counts == 0)
            totalvar[inds] = -32768
            inds = numpy.where(counts != 0)
            totalvar[inds] = numpy.divide(totalvar[inds], counts[inds])
            imgout[1].Write(totalvar)
            t = datetime.datetime.now()-start
            VerboseOut('%s: mean + variance for %s files processed in %s' % (os.path.basename(fout), len(filenames), t))
        return imgout

    def _read_point(self, filename, roi):
        """ Read single point from mean/var file and return if valid, or mean/var of 3x3 neighborhood """
        img = gippy.GeoImage(filename)
        vals = img[0].Read(roi).squeeze()
        variances = img[1].Read(roi)
        if numpy.isnan(vals[1, 1]):
            val = numpy.nanmean(vals)
            var = numpy.nanmean(variances)
        else:
            val = vals[1, 1]
            var = variances[1, 1]
        img = None
        return val, var

    def get_point(self, lat, lon, product='aod'):
        pixx = int(numpy.round(float(lon) + 179.5))
        pixy = int(numpy.round(89.5 - float(lat)))
        roi = gippy.iRect(pixx-1, pixy-1, 3, 3)
        #img = self.open(product=product)
        img = gippy.GeoImage(self.products[product], False)
        nodata = img[0].NoDataValue()
        vals = img[0].Read(roi).squeeze()
        # TODO - do this automagically in swig wrapper
        vals[numpy.where(vals == nodata)] = numpy.nan

        aod = vals[1, 1]
        source = 'actual'
        if numpy.isnan(aod):
            aod = numpy.nanmean(vals)
            source = 'actual spatial average'

        day = self.date.strftime('%j')

        # Calculate best estimate from multiple sources
        if numpy.isnan(aod):
            aod = 0.0
            norm = 0.0
            cnt = 0

            source = 'best estimate'
            # LTA-Daily
            filename = os.path.join(self.Repository.cpath('ltad'), 'ltad%s.tif' % str(day).zfill(3))
            val, var = self._read_point(filename, roi)
            if not numpy.isnan(val):
                aod = val/var
                norm = var
                cnt = cnt + 1
            VerboseOut('AOD: LTA-Daily = %s, %s' % (val, var), 4)

            # LTA
            val, var = self._read_point(self.Repository.cpath('lta.tif'), roi)
            if not numpy.isnan(val):
                aod = aod + val/var
                norm = norm + var
                cnt = cnt + 1
            VerboseOut('AOD: LTA = %s, %s' % (val, var), 4)

            # TODO - adjacent days

            # Final AOD estimate
            aod = aod/norm
            VerboseOut('AOD = %s, var=%s' % (aod, norm/cnt))

        if numpy.isnan(aod):
            aod = 0.17
            source = 'default'

        VerboseOut('AOD (%s) = %s' % (source, aod), 3)
        return aod


def main():
    DataInventory.main(AODData)
