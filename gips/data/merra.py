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
import datetime
import time
from pydap.client import open_url
import numpy

import gippy
from gips.data.core import Repository, Asset, Data
from gips.utils import VerboseOut, basename

requirements = ['pydap']


class MerraRepository(Repository):
    name = 'Merra'
    description = 'Modern Era Retrospective-Analysis for Research and Applications (weather and climate)'

    @classmethod
    def tile_bounds(cls, tile):
        """ Get the bounds of the tile (in same units as tiles vector) """
        tilesvector = cls.vector()
        for fid in tilesvector.get_fids():
            feature = tilesvector.get_feature(fid)
            if feature['tileid'] == tile:
                break
        bounds = eval(feature['bounds'])
        return bounds


class MerraAsset(Asset):
    Repository = MerraRepository

    _sensors = {
        'MERRA': {
            'description': 'Modern Era Retrospective-analysis for Resarch and Application',
        }
    }

    _bandnames = ['%02d30GMT' % i for i in range(24)]

    _assets = {
        'TS': {
            'description': 'Surface skin temperature',
            'pattern': 'MERRA_TS_*.tif',
            'url': 'http://goldsmr2.sci.gsfc.nasa.gov:80/opendap/MERRA/MAT1NXSLV.5.2.0',
            'source': 'MERRA%s.prod.assim.tavg1_2d_slv_Nx.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 11),
            'latency': 60,
            'bandnames': _bandnames
        },
        'T2M': {
            'description': 'Temperature at 2 m above the displacement height',
            'pattern': 'MERRA_T2M_*.tif',
            'url': 'http://goldsmr2.sci.gsfc.nasa.gov:80/opendap/MERRA/MAT1NXSLV.5.2.0',
            'source': 'MERRA%s.prod.assim.tavg1_2d_slv_Nx.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 11),
            'latency': 60,
            'bandnames': _bandnames
        },
        'T10M': {
            'description': 'Temperature at 10 m above the displacement height',
            'pattern': 'MERRA_T10M_*.tif',
            'url': 'http://goldsmr2.sci.gsfc.nasa.gov:80/opendap/MERRA/MAT1NXSLV.5.2.0',
            'source': 'MERRA%s.prod.assim.tavg1_2d_slv_Nx.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 11),
            'latency': 60,
            'bandnames': _bandnames
        },
        #'PREC': {
        #    'description': 'Total Precipitation (kg m-2 s-1)',
        #    'pattern': 'MST1NXMLD_PREC_*.tif',
        #    'url': 'http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MST1NXMLD.5.2.0',
        #    'source': 'MERRA%s.prod.simul.tavg1_2d_mld_Nx.%04d%02d%02d.hdf',
        #    'startdate': datetime.date(1980, 1, 1),
        #    'latency': 60
        #},
        'PROFILE': {
            'description': 'Atmospheric Profile',
            'pattern': 'MAI6NVANA_PROFILE_*.tif',
            'url': 'http://goldsmr3.sci.gsfc.nasa.gov:80/opendap/MERRA/MAI6NVANA.5.2.0',
            'source': 'MERRA%s.prod.assim.inst6_3d_ana_Nv.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 1),
            'latency': 60,
            'bandnames': ['0000GMT', '0600GMT', '1200GMT', '1800GMT']
        },
        'PROFILEP': {
            'description': 'Atmospheric Profile',
            'pattern': 'MAI6NVANA_PROFILE_*.tif',
            'url': 'http://goldsmr3.sci.gsfc.nasa.gov:80/opendap/MERRA/MAI6NPANA.5.2.0',
            'source': 'MERRA%s.prod.assim.inst6_3d_ana_Np.%04d%02d%02d.hdf',
            'startdate': datetime.date(1979, 1, 1),
            'latency': 60,
            'bandnames': ['0000GMT', '0600GMT', '1200GMT', '1800GMT']
        }
    }

    _origin = (-180., -90.)
    _defaultresolution = (0.666666666666667, 0.50)

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(MerraAsset, self).__init__(filename)
        parts = basename(filename).split('_')
        self.sensor = 'MERRA'
        self.asset = parts[1]
        self.tile = parts[2]
        # e.g., ['MERRA', 'TS', 'h06v05', '2010001']
        self.date = datetime.datetime.strptime(parts[3], '%Y%j').date()
        self.products[self.asset] = filename

    @classmethod
    def opendap_fetch(cls, asset, date):
        """ Get array proxy from OpenDap for this asset and date """
        url = cls._assets[asset].get('url', '')
        if url == '':
            raise Exception("%s: URL not defined for asset %s" % (cls.__name__, asset))
        success = False
        for ver in ['100', '200', '300', '301']:
            f = cls._assets[asset]['source'] % (ver, date.year, date.month, date.day)
            loc = "%s/%04d/%02d/%s" % (url, date.year, date.month, f)
            try:
                dataset = open_url(loc)
            except:
                continue
            else:
                success = True
                break
        if not success:
            raise Exception('Data unavailable (%s)' % loc)
        return dataset

    @classmethod
    def lonlat2xy(cls, lon, lat):
        """ Convert from lon-lat to x-y in array """
        x = int(round((lon - cls._origin[0]) / cls._defaultresolution[0]))
        y = int(round((lat - cls._origin[1]) / cls._defaultresolution[1]))
        return (x, y)

    @classmethod
    def fetch(cls, asset, tile, date):
        """ Get this asset for this tile and date (using OpenDap service) """
        #super(MerraAsset, cls).fetch(asset, tile, date)

        dataset = cls.opendap_fetch(asset, date)

        # Find the bounds of the tile requested
        bounds = cls.Repository.tile_bounds(tile)
        # TODO - get origin from shapefile
        ORIGIN = (-180., -90.)
        dx = bounds[2] - bounds[0]
        dy = bounds[3] - bounds[1]
        xsize = int(round(dx / cls._defaultresolution[0]))
        ysize = int(round(dy / cls._defaultresolution[1]))

        # verify bounds
        #h = int(tile[1:3])
        #v = int(tile[4:6])
        #x0_alt = ORIGIN[0] + xsize * (h - 1)
        #y0_alt = ORIGIN[1] + ysize * (18 - v)
        #assert x0 == x0_alt
        #assert y0 == y0_alt
        ix0 = int(round((bounds[0] - ORIGIN[0]) / cls._defaultresolution[0]))
        iy0 = int(round((bounds[1] - ORIGIN[1]) / cls._defaultresolution[1]))
        ix1 = ix0 + xsize
        iy1 = iy0 + ysize

        VerboseOut('Retrieving data for bounds (%s, %s) - (%s, %s)' % (bounds[0], bounds[1], bounds[2], bounds[3]), 3)

        # TODO - what keys to get?
        shape = dataset['PS'].shape

        print dataset.keys()
        print shape

        data = dataset[asset][:, iy0:iy1, ix0:ix1].astype('float32')
        # What's this for?
        data = data[:, ::-1, :]

        # Save tile data
        description = cls._assets[asset]['description']
        meta = {'ASSET': asset, 'TILE': tile, 'DATE': str(date.date()), 'DESCRIPTION': description}
        doy = date.strftime('%j')
        fout = os.path.join(cls.Repository.spath(), "MERRA_%s_%s_%4d%s.tif" % (asset, tile, date.year, doy))
        # TODO - use GIPPY to write
        from agspy.utils import raster
        proj = raster.create_proj(4326)
        geo = (bounds[0], cls._defaultresolution[0], 0.0, bounds[3], 0.0, -cls._defaultresolution[1])
        raster.write_raster(fout, data, proj, geo, meta, bandnames=cls._assets[asset]['bandnames'])


class MerraData(Data):
    """ A tile of data (all assets and products) """
    name = 'Merra'
    version = '0.9.0'
    Asset = MerraAsset

    _products = {
        'temp': {
            'description': 'Air temperature data',
            'assets': ['TS', 'T2M', 'T10M']
        },
        'T2M': {
            'description': 'Temperature at 2 m above the displacement height',
            'assets': ['T2M']
        },
        'T10M': {
            'description': 'Temperature at 10 m above the displacement height',
            'assets': ['T10M']
        },
        'TS': {
            'description': 'Surface temperature',
            'assets': ['TS']
        }

        #'daily_weather': {
        #    'description': 'Climate forcing data, e.g. for DNDC',
        #    'assets': ['T2M', 'PRECTOT']
        #},
        #'profile': {
        #    'description': 'Atmospheric Profile',
        #    'assets': ['PROFILE'],
        #}

    }

    @classmethod
    def profile(cls, lon, lat, dtime):
        """ Retrieve atmospheric profile directly from merra data via OpenDap """
        dataset = cls.Asset.opendap_fetch('PROFILE', dtime)
        (x, y) = cls.Asset.lonlat2xy(lon, lat)

        # TODO - I know these are hours (0, 6, 12, 18), but it's still an assumption
        times = [datetime.datetime.combine(dtime.date(), datetime.time(int(d / 60.0))) for d in dataset['TIME'][:]]
        unixtime = time.mktime(dtime.timetuple())
        timediff = numpy.array([unixtime - time.mktime(t.timetuple()) for t in times])
        timeind = numpy.abs(timediff).argmin()

        p = dataset['PS'][timeind, y, x].squeeze()
        pthick = dataset['DELP'][timeind, :, y, x].squeeze()[::-1]
        pressure = []
        for pt in pthick:
            pressure.append(p)
            p = p - pt
        pressure = numpy.array(pressure)
        inds = numpy.where(pressure > 0)

        data = {
            # Pa -> mbar
            'pressure': numpy.array(pressure)[inds] / 100.0,
            # Kelvin -> Celsius
            'temp': dataset['T'][timeind, :, y, x].squeeze()[::-1][inds] - 273.15,
            # kg/kg -> g/kg (Mass mixing ratio)
            'humidity': dataset['QV'][timeind, :, y, x].squeeze()[::-1][inds] * 1000,
            'ozone': dataset['O3'][timeind, :, y, x].squeeze()[::-1][inds] * 1000,
        }

        return data

    def getlonlat(self):
        """ return the center coordinates of the MERRA tile """
        hcoord = int(self.id[1:3])
        vcoord = int(self.id[4:])
        lon = -180. + (hcoord - 1) * 12. + 6.
        lat = 90. - vcoord * 10. - 5.
        return lon, lat

    def gmtoffset(self):
        """ return the approximate difference between local time and GMT """
        lon = self.getlonlat()[0]
        houroffset = lon * (12. / 180.)
        return houroffset

    def process(self, *args, **kwargs):
        products = super(MerraData, self).process(*args, **kwargs)
        if len(products) == 0:
            return

        self.basename = self.basename + '_' + self.sensor_set[0]
        for key, val in products.requested.items():
            try:
                assets = self.asset_filenames(val[0])
            except:
                # Required assets unavailable, continue to next product
                continue
            fout = os.path.join(self.path, self.basename + '_' + key)

            ####################################################################
            """
            # This doesn't currently save anything
            if val[0] == "daily_weather":
                t10 = img[0].Read()
                prectot = img[1].Read()

                print t10.shape
                print t10.min(), t10.mean(), t10.max()

                print prectot.shape
                print prectot.min(), prectot.mean(), prectot.max()
            """
            ####################################################################
            if val[0] == "temp":
                img = gippy.GeoImage(assets[0])

                imgout = gippy.GeoImage(fout, img, img.DataType(), 4)

                # Aqua AM, Terra AM, Aqua PM, Terra PM
                localtimes = [1.5, 10.5, 13.5, 22.5]
                strtimes = ['0130LT', '1030LT', '1330LT', '2230LT']
                hroffset = self.gmtoffset()

                # TODO: loop across the scene in latitude
                # calculate local time for each latitude column
                print 'localtimes', localtimes
                for itime, localtime in enumerate(localtimes):
                    print itime
                    picktime = localtime - hroffset
                    pickhour = int(picktime)

                    if pickhour < 0:
                        # next day local time
                        pickday = +1
                    elif pickhour > 24:
                        # previous day local time
                        pickday = -1
                    else:
                        # same day local time
                        pickday = 0

                    pickidx = pickhour % 24
                    print "localtime", localtime
                    print "picktime", picktime
                    print "pickhour", pickhour
                    print "pickday", pickday
                    print "pickidx", pickidx

                    img[pickidx].Process(imgout[itime])

                    obsdate = self.date + datetime.timedelta(pickday)
                    descr = " ".join([strtimes[itime], obsdate.isoformat()])
                    imgout.SetBandName(descr, itime + 1)
            ####################################################################
            elif val[0] == 'profile':
                pass
