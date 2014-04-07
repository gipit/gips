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

from gippy.data.inventory import DataInventory
from gippy.data.core import Repository, Asset, Data
from gippy.utils import File2List, List2File, VerboseOut
import gippy.settings as settings

import traceback
from pdb import set_trace


class AtmosRepository(Repository):
    _rootpath = '/titan/data/atmos'
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


class AtmosAsset(Asset):
    Repository = ModisAtmosRepository

    # ???? Not specific to MODIS
    _sensors = {
        'MOD': {'description': 'Terra'},
        'MYD': {'description': 'Aqua'},
        'MCD': {'description': 'Combined'}
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
        super(ModisAtmosAsset, self).__init__(filename)

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
    Asset = ModisAtmosAsset

    _products = {
        'aero': {
            'description': 'Aerosols',
            # the list of asset types associated with this product
            'assets': ['MOD08'],  # , 'MYD08'],
        },
    }


def main():
    DataInventory.main(AtmosData)
