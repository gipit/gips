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

import gippy
from gippy.data.inventory import DataInventory
from gippy.data.core import Repository, Asset, Data
from gippy.utils import File2List, List2File, VerboseOut
import gippy.settings as settings

import traceback
from pdb import set_trace


class AtmosRepository(Repository):
    _rootpath = '/titan/data/atmos'
    _cpath = 'composites'
    _datedir = '%Y%j'

    @classmethod
    def path(cls, tile='', date=''):
        path = os.path.join(cls._rootpath, cls._tilesdir)
        if date != '':
            path = os.path.join(path, str(date.strftime('%Y')), str(date.strftime('%j')))
        return path

    @classmethod
    def cpath(cls):
        return os.path.join(cls._rootpath, cls._cpath)

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


class AtmosAsset(Asset):
    Repository = AtmosRepository

    # ???? Not specific to MODIS
    _sensors = {
        'MOD': {'description': 'MODIS Terra'},
        'MYD': {'description': 'MODIS Aqua'},
    }
    _assets = {
        'MOD08': {
            'pattern': 'MOD08_D3*hdf',
            'url': 'ladsweb.nascom.nasa.gov/allData/51/MOD08_D3'
        },
        'MYD08': {
            'pattern': 'MYD08_D3*hdf',
            'url': 'ladsweb.nascom.nasa.gov/allData/51/MYD08_D3'
        }
    }

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(AtmosAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        self.asset = bname[0:5]
        self.tile = ''
        year = bname[10:14]
        doy = bname[14:17]
        self.date = datetime.datetime.strptime(year+doy, "%Y%j").date()
        self.sensor = bname[:3]
        #datafiles = self.datafiles()
        prefix = 'HDF4_EOS:EOS_GRID:"'
        self.products = {'aero': prefix + filename + '":mod08:Optical_Depth_Land_And_Ocean_Mean'}

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
        url = cls._assets[asset]['url']
        ftpurl = url.split('/')[0]
        ftpdir = url[len(ftpurl):]

        try:
            ftp = ftplib.FTP(ftpurl)
            ftp.login('anonymous', settings.EMAIL)
            ftp.cwd(os.path.join(ftpdir, date.strftime('%Y'), date.strftime('%j')))
            ftp.set_pasv(True)

            filenames = []
            ftp.retrlines('LIST', filenames.append)

            for f in ftp.nlst('*.hdf'):
                VerboseOut("Downloading %s" % f, 3)
                ftp.retrbinary('RETR %s' % f, open(os.path.join(cls.Repository.spath(), f), "wb").write)

            ftp.close()
        except:
            VerboseOut(traceback.format_exc(), 4)
            raise Exception("Error downloading")


class AtmosData(Data):
    name = 'Globally Gridded Atmospheric Data'
    Asset = AtmosAsset

    _products = {
        'aero': {
            'description': 'Aerosols',
            # the list of asset types associated with this product
            'assets': ['MOD08'],  # , 'MYD08'],
        },
        'aerolta': {
            'description': 'Aerosols, average daily',
            # the list of asset types associated with this product
            'assets': ['MOD08'],  # , 'MYD08'],
        },
    }

    #@classmethod
    def process_aerolta(cls):
        """ Calculate long-term averages of aerosol optical thickness """
        fouts = []
        for day in range(1, 366):
            start = datetime.datetime.now()
            inv = cls.inventory(products=['aero'], days="%s,%s" % (day, day))
            fnames = [inv[d].tiles[''].products['aero'] for d in inv.dates]
            if len(fnames) > 0:
                img = gippy.GeoImage(fnames)
                fout = os.path.join(cls.Repository.cpath(), 'aerolta', 'aerolta_%s.tif' % day)
                aerolta = img.Mean(fout)
                fouts.append(fout)
                t = datetime.datetime.now()-start
                VerboseOut('Long-term average aerosol optical depth: processed day %s in %s' % (day, t))
        start = datetime.datetime.now()
        aero = gippy.GeoImage(fouts)
        fout = os.path.join(cls.Repository.cpath(), 'aerolta', 'aerolta.tif')
        aerolta.Mean(fout)
        t = datetime.datetime.now()-t
        VerboseOut('Long-term average aerosol optical depth: processed in %s' % (t))
        # calculate long term average

    def process(self, products):
        start = datetime.datetime.now()
        #bname = os.path.basename(self.assets[''].filename)
        for product in products:
            if product == 'aerolta':
                self.process_aerolta()
            VerboseOut(' -> %s: processed %s in %s' % (fout, product, datetime.datetime.now()-start))

    def get_point(self, lat, lon, product=''):
        pixx = int(numpy.round(float(lon) + 179.5))
        pixy = int(numpy.round(89.5 - float(lat)))
        roi = gippy.iRect(pixx, pixy, 0, 0)
        img = self.open(product=product)
        val = img[0].Read(roi).squeeze()
        day = self.date.strftime('%j')
        print 'val', val
        if val == img[0].NoData():
            fname = os.path.join(self.Repository.cpath(), 'aerolta', 'aerolta_%s' % day)
            img = gippy.GeoImage(fname)
            val = img[0].Read(roi).squeeze()
            print 'val', val
        if val == img[0].NoData():
            print 'still nodata aerosols'
            val = 0.17
        return val


def main():
    DataInventory.main(AtmosData)
