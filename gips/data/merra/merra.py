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
import signal

import gippy
from gips.data.core import Repository, Asset, Data
from gips.utils import VerboseOut, basename, open_vector

# TODO: use gippy instead
from gips.data.merra import raster


from pdb import set_trace

requirements = ['pydap']

MISSING = 9.999999870e+14

class Timeout():
    """Timeout class using ALARM signal."""

    class Timeout(Exception):
        pass

    def __init__(self, sec):
        self.sec = sec

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.raise_timeout)
        signal.alarm(self.sec)

    def __exit__(self, *args):
        signal.alarm(0)    # disable alarm

    def raise_timeout(self, *args):
        raise Timeout.Timeout()


class merraRepository(Repository):
    name = 'merra'
    description = 'Modern Era Retrospective-Analysis for Research and Applications (weather and climate)'
    _tile_attribute = 'tileid'

    @classmethod
    def tile_bounds(cls, tile):
        """ Get the bounds of the tile (in same units as tiles vector) """

        # vector = open_vector(cls.get_setting('tiles'))
        # extent = vector.where('tileid==%s' % tile).extent()
        # return [extent.x0(), extent.y0(), extent.x1(), extent.y1()]

        vector = open_vector(cls.get_setting('tiles'))
        features = vector.where('tileid', tile)
        if len(features) != 1:
            raise Exception('there should be a single tile with id %s' % tile)

            # set_trace()

        extent = features[0].Extent()
        return [extent.x0(), extent.y0(), extent.x1(), extent.y1()]



class merraAsset(Asset):
    Repository = merraRepository

    _sensors = {
        'merra': {
            'description': 'Modern Era Retrospective-analysis for Research and Applications',
        }
    }

    _bandnames = ['%02d30GMT' % i for i in range(24)]

    _assets = {
        # MERRA2 SLV
        'TS': {
            'description': 'Surface skin temperature',
            'pattern': 'MERRA_TS_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_slv_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        'T2M': {
            'description': 'Temperature at 2 m above the displacement height',
            'pattern': 'MERRA_T2M_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_slv_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        'T10M': {
            'description': 'Temperature at 10 m above the displacement height',
            'pattern': 'MERRA_T10M_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_slv_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        'PS': {
            'description': 'Time averaged surface pressure (Pa)',
            'pattern': 'MERRA_PS_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_slv_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        'QV2M': {
            'description': 'Specific humidity at 2 m above the displacement height (kg kg-1)',
            'pattern': 'MERRA_QV2M_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_slv_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        # MERRA2 FLX
        'PRECTOT': {
            'description': 'Total Precipitation (kg m-2 s-1)',
            'pattern': 'MERRA_PRECTOT_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_flx_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        'SPEED': {
            'description': '3-dimensional wind speed for surface fluxes (m s-1)',
            'pattern': 'MERRA_SPEED_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_flx_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
        # MERRA2 RAD
        'SWGDN': {
            'description': 'Surface incident shortwave flux (W m-2)',
            'pattern': 'MERRA_SWGDN_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov/opendap/MERRA2/M2T1NXRAD.5.12.4',
            'source': 'MERRA2_%s.tavg1_2d_rad_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
            'bandnames': _bandnames
        },
         'PROFILE': {
             'description': 'Atmospheric Profile',
             'pattern': 'MAI6NVANA_PROFILE_*.tif',
             'url': 'http://goldsmr5.sci.gsfc.nasa.gov/opendap/MERRA2/M2I6NVANA.5.12.4',
             'source': 'MERRA2_%s.inst6_3d_ana_Nv.%04d%02d%02d.nc4',
             'startdate': datetime.date(1980, 1, 1),
             'latency': 60,
             'bandnames': ['0000GMT', '0600GMT', '1200GMT', '1800GMT']
         },
        # 'PROFILEP': {
        #     'description': 'Atmospheric Profile',
        #     'pattern': 'MAI6NVANA_PROFILE_*.tif',
        #     'url': 'http://goldsmr3.sci.gsfc.nasa.gov:80/opendap/MERRA/MAI6NPANA.5.2.0',
        #     'source': 'MERRA%s.prod.assim.inst6_3d_ana_Np.%04d%02d%02d.hdf',
        #     'startdate': datetime.date(1980, 1, 1),
        #     'latency': 60,
        #     'bandnames': ['0000GMT', '0600GMT', '1200GMT', '1800GMT']
        # },
        # MERRA2 CONST
        'FRLAND': {
            'description': 'Land Fraction',
            'pattern': 'MERRA_FRLAND_*.tif',
            'url': 'http://goldsmr4.sci.gsfc.nasa.gov:80/opendap/MERRA2_MONTHLY/M2C0NXASM.5.12.4',
            'source': 'MERRA2_%s.const_2d_asm_Nx.%04d%02d%02d.nc4',
            'startdate': datetime.date(1980, 1, 1),
            'latency': 0,
            'bandnames': ['FRLAND']
        }

    }

    _origin = (-180., -90.)
    _defaultresolution = (0.625, 0.50)
    # _defaultresolution = (0.666666666666667, 0.50)

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(merraAsset, self).__init__(filename)
        parts = basename(filename).split('_')
        self.sensor = 'merra'
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

        for ver in ['100', '200', '300', '301', '400']:

            if asset != "FRLAND":
                f = cls._assets[asset]['source'] % (ver, date.year, date.month, date.day)
                loc = "%s/%04d/%02d/%s" % (url, date.year, date.month, f)
            else:
                f = cls._assets[asset]['source'] % (ver, 0, 0, 0)
                loc = "%s/1980/%s" % (url, f)
            try:
                with Timeout(30):
                    dataset = open_url(loc)
            except Timeout.Timeout:
                print "Timeout"
            except Exception,e:
                pass
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

        if cls._assets[asset]['latency'] is None:
            assert date == datetime.datetime(1980, 1, 1)

        # if date > datetime.datetime.now():
        #     print "These data are not available for future dates."
        #     return None

        if date > (datetime.datetime.now() - datetime.timedelta(cls._assets[asset]['latency'])):
            print "These data are not available for specified dates."
            return None

        try:
            dataset = cls.opendap_fetch(asset, date)
        except Exception:
            print "Fetch: data not available", asset, tile, date
            return

        # Find the bounds of the tile requested
        bounds = cls.Repository.tile_bounds(tile)
        # TODO - get origin from shapefile
        ORIGIN = (-180., -90.)
        dx = bounds[2] - bounds[0]
        dy = bounds[3] - bounds[1]
        xsize = int(round(dx / cls._defaultresolution[0]))
        ysize = int(round(dy / cls._defaultresolution[1]))

        ix0 = int(round((bounds[0] - ORIGIN[0]) / cls._defaultresolution[0]))
        iy0 = int(round((bounds[1] - ORIGIN[1]) / cls._defaultresolution[1]))
        ix1 = ix0 + xsize
        iy1 = iy0 + ysize

        VerboseOut('Retrieving data for bounds (%s, %s) - (%s, %s)' % (bounds[0], bounds[1], bounds[2], bounds[3]), 3)

        data = dataset[asset][:, iy0:iy1, ix0:ix1].astype('float32')
        data = data[:, ::-1, :]

        # Save tile data
        description = cls._assets[asset]['description']
        meta = {'ASSET': asset, 'TILE': tile, 'DATE': str(date.date()), 'DESCRIPTION': description}
        doy = date.strftime('%j')
        fout = os.path.join(cls.Repository.path('stage'), "MERRA_%s_%s_%4d%s.tif" % (asset, tile, date.year, doy))

        # TODO: use gippy instead
        proj = raster.create_proj(4326)
        geo = (bounds[0], cls._defaultresolution[0], 0.0, bounds[3], 0.0, -cls._defaultresolution[1])
        print "writing", fout
        raster.write_raster(fout, data, proj, geo, meta, bandnames=cls._assets[asset]['bandnames'])


class merraData(Data):
    """ A tile of data (all assets and products) """
    name = 'merra'
    version = '0.9.0'
    Asset = merraAsset

    _products = {
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
        },
        'PRECTOT': {
            'description': 'Precipitation rate [kg m-2 s-1]',
            'assets': ['PRECTOT']
        },
        'SPEED': {
            'description': '3-dimensional wind speed for surface fluxes (m s-1)',
            'assets': ['SPEED']
        },
        'FRLAND': {
            'description': 'Land Fraction (static map)',
            'assets': ['FRLAND']
        },

        'tave': {
            'description': 'Ave daily air temperature data',
            'assets': ['T2M']
        },
        'tmin': {
            'description': 'Min daily air temperature data',
            'assets': ['T2M']
        },
        'tmax': {
            'description': 'Max daily air temperature data',
            'assets': ['T2M']
        },
        'prcp': {
            'description': 'Daily total precipitation (mm day-1)',
            'assets': ['PRECTOT']
        },
        'wind': {
            'description': 'Daily mean wind speed (m s-1)',
            'assets': ['SPEED']
        },
        'shum': {
            'description': 'Relative humidity (kg kg-1)',
            'assets': ['QV2M']
        },
        'srad': {
            'description': 'Incident solar radiation (W m-2)',
            'assets': ['SWGDN']
        },
        'patm': {
            'description': 'Surface atmospheric pressure (mb)',
            'assets': ['PS']
        },
        'rhum': {
            'description': 'Relative humidity (%)',
            'assets': ['QV2M', 'PS', 'T2M']
        },

        '_temps': {
            'description': 'Air temperature data',
            'assets': ['TS', 'T2M', 'T10M']
        },


        #'daily_weather': {
        #    'description': 'Climate forcing data, e.g. for DNDC',
        #    'assets': ['T2M', 'PRECTOT']
        #},
        #'profile': {
        #    'description': 'Atmospheric Profile',
        #    'assets': ['PROFILE'],
        #}

    }


    # @classmethod
    # def process_composites(cls, inventory, products, **kwargs):
    #     for product in products:
    #         cpath = os.path.join(cls.Asset.Repository.path('composites'), 'ltad')
    #         path = os.path.join(cpath, 'ltad')
    #         # Calculate AOT long-term multi-year averages (lta) for given day
    #         if product == 'ltad':
    #             for day in range(inventory.start_day, inventory.end_day + 1):
    #                 dates = [d for d in inventory.dates if int(d.strftime('%j')) == day]
    #                 filenames = [inventory[d].tiles[''].products['aod'] for d in dates]
    #                 fout = path + '%s.tif' % str(day).zfill(3)
    #                 cls.process_mean(filenames, fout)
    #         # Calculate single average per pixel (all days and years)
    #         if product == 'lta':
    #             filenames = glob.glob(path + '*.tif')
    #             if len(filenames) > 0:
    #                 fout = os.path.join(cls.Asset.Repository.path('composites'), 'lta.tif')
    #                 cls.process_mean(filenames, fout)
    #             else:
    #                 raise Exception('No daily LTA files exist!')


    @classmethod
    def profile(cls, lon, lat, dtime):
        """ Retrieve atmospheric profile directly from merra data via OpenDap """
        dataset = cls.Asset.opendap_fetch('PROFILE', dtime)
        (x, y) = cls.Asset.lonlat2xy(lon, lat)

        # TODO - I know these are hours (0, 6, 12, 18), but it's still an assumption
        times = [datetime.datetime.combine(dtime.date(), datetime.time(int(d / 60.0))) for d in dataset['time'][:]]
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
        products = super(merraData, self).process(*args, **kwargs)
        if len(products) == 0:
            return
        for key, val in products.requested.items():
            try:
                assets = self.asset_filenames(val[0])
            except:
                # Required assets unavailable, continue to next product
                continue
            fout = os.path.join(self.path, self.basename + "_merra_" + key)

            ####################################################################

            if val[0] == "tave":
                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                thourly = numpy.ma.MaskedArray(img.Read().squeeze())
                thourly.mask = (thourly == MISSING)

                temp = thourly.mean(axis=0)
                temp = temp - 273.16

                temp[temp.mask] = MISSING
                imgout[0].Write(numpy.array(temp))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('C')
                imgout.SetNoData(MISSING)

            ####################################################################

            if val[0] == "tmax":
                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                thourly = numpy.ma.MaskedArray(img.Read().squeeze())
                thourly.mask = (thourly == MISSING)

                temp = thourly.max(axis=0)
                temp = temp - 273.16

                temp[temp.mask] = MISSING
                imgout[0].Write(numpy.array(temp))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('C')
                imgout.SetNoData(MISSING)

            ####################################################################

            if val[0] == "tmin":
                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                thourly = numpy.ma.MaskedArray(img.Read().squeeze())
                thourly.mask = (thourly == MISSING)

                temp = thourly.min(axis=0)
                temp = temp - 273.16

                temp[temp.mask] = MISSING
                imgout[0].Write(numpy.array(temp))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('C')
                imgout.SetNoData(MISSING)

            ####################################################################

            if val[0] == "prcp":

                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                phourly = numpy.ma.MaskedArray(img.Read().squeeze())
                phourly.mask = (phourly == MISSING)

                prcp = phourly.mean(axis=0)
                prcp = prcp * 36. * 24. * 1000

                prcp[prcp.mask] = MISSING
                imgout[0].Write(numpy.array(prcp))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('mm d-1')
                imgout.SetNoData(MISSING)


            ####################################################################

            if val[0] == "wind":

                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                hourly = numpy.ma.MaskedArray(img.Read().squeeze())
                hourly.mask = (hourly == MISSING)
                daily = hourly.mean(axis=0)
                # daily = daily * 1.0
                daily[daily.mask] = MISSING
                imgout[0].Write(numpy.array(daily))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('m s-1')
                imgout.SetNoData(MISSING)

            ####################################################################

            if val[0] == "srad":

                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                hourly = numpy.ma.MaskedArray(img.Read().squeeze())
                hourly.mask = (hourly == MISSING)
                daily = hourly.mean(axis=0)
                # daily = daily * 1.0
                daily[daily.mask] = MISSING
                imgout[0].Write(numpy.array(daily))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('W m-2')
                imgout.SetNoData(MISSING)

            ####################################################################

            if val[0] == "shum":

                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                hourly = numpy.ma.MaskedArray(img.Read().squeeze())
                hourly.mask = (hourly == MISSING)
                daily = hourly.mean(axis=0)
                # daily = daily * 1.0
                daily[daily.mask] = MISSING
                imgout[0].Write(numpy.array(daily))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('kg kg-1')
                imgout.SetNoData(MISSING)

            ####################################################################

            if val[0] == "patm":

                img = gippy.GeoImage(assets[0])
                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                hourly = numpy.ma.MaskedArray(img.Read().squeeze())
                hourly.mask = (hourly == MISSING)
                daily = hourly.mean(axis=0)
                daily = daily / 100.
                daily[daily.mask] = MISSING
                imgout[0].Write(numpy.array(daily))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('mb')
                imgout.SetNoData(MISSING)


            ####################################################################

            # ##' @param qair specific humidity, dimensionless (e.g. kg/kg) ratio of water mass / total air mass
            # ##' @param temp degrees C
            # ##' @param press pressure in mb
            # ##' @return rh relative humidity, ratio of actual water mixing ratio to saturation mixing ratio
            # ##' @export
            # ##' @author David LeBauer
            # qair2rh <- function(qair, temp, press = 1013.25){
            #     es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
            #     e <- qair * press / (0.378 * qair + 0.622)
            #     rh <- e / es
            #     rh[rh > 1] <- 1
            #     rh[rh < 0] <- 0
            #     return(rh)
            # }


            if val[0] == "rhum":
                # based on 'QV2M', 'PS', 'T2M

                img = gippy.GeoImage(assets[0])
                qv2m = numpy.ma.MaskedArray(img.Read().squeeze()) # kg kg-1
                qv2m.mask = (qv2m == MISSING)
                img = gippy.GeoImage(assets[1])
                ps = numpy.ma.MaskedArray(img.Read().squeeze()) # Pa
                ps.mask = (ps == MISSING)
                img = gippy.GeoImage(assets[2])
                t2m = numpy.ma.MaskedArray(img.Read().squeeze()) # K
                t2m.mask = (t2m == MISSING)

                temp = t2m - 273.15
                press = ps/100.
                qair = qv2m

                es = 6.112*numpy.exp((17.67*temp)/(temp + 243.5))

                e = qair*press/(0.378*qair + 0.622)
                rh = 100. * (e/es)
                rh[rh > 100.] = 100.
                rh[rh < 0.] = 0.

                daily = rh.mean(axis=0)
                daily[daily.mask] = MISSING

                imgout = gippy.GeoImage(fout, img, img.DataType(), 1)
                imgout[0].Write(numpy.array(daily))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('%')
                imgout.SetNoData(MISSING)


            ####################################################################

            if val[0] == "temp_modis":

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
