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

from __future__ import print_function

import os
import re
import sys
import datetime
import time

import tempfile
import urllib
import requests
import signal

import numpy as np
from netCDF4 import Dataset

import gippy
from gips.data.core import Repository, Asset, Data
from gips.utils import basename, open_vector
from gips import utils


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

    # NASA assets require special authentication
    _manager_url = "https://urs.earthdata.nasa.gov"

    @classmethod
    def tile_bounds(cls, tile):
        """ Get the bounds of the tile (in same units as tiles vector) """
        vector = open_vector(cls.get_setting('tiles'))
        features = vector.where('tileid', tile)
        if len(features) != 1:
            raise Exception('there should be a single tile with id %s' % tile)
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
    # used in _assets[asset_type]['pattern'], which is used by data/core.py to search the filesystem
    _asset_re_pattern = '^MERRA2_\d\d\d\.%s\.\d{4}\d{2}\d{2}\.nc4$'
    # used in _assets[asset_type]['re_pattern'], which is used exclusively by query_service
    _asset_re_format_pattern = "MERRA2_\d\d\d\.{name}\.%04d%02d%02d\.nc4"

    _assets = {
        # MERRA2 SLV
        ## TS: Surface skin temperature (K)
        ## T2M: Temperature at 2 m above the displacement height (K)
        ## T10M: Temperature at 10 m above the displacement height (K)
        ## PS: Time averaged surface pressure (Pa)
        ## QV2M: Specific humidity at 2 m above the displacement height (kg kg-1)
        'SLV': {
            'shortname': 'M2T1NXSLV',
            'description': '2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Single-Level Diagnostics V5.12.4',
            'url': 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4',
            'pattern': _asset_re_pattern % 'tavg1_2d_slv_Nx',
            're_pattern': _asset_re_format_pattern.format(name='tavg1_2d_slv_Nx'),
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
        },
        # MERRA2 FLX
        ## PRECTOT: Total Precipitation (kg m-2 s-1)'
        ## SPEED: 3-dimensional wind speed for surface fluxes (m s-1)'
        'FLX': {
            'shortname': 'M2T1NXFLX',
            'description': '2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux Diagnostics V5.12.4',
            'url': 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4',
            'pattern': _asset_re_pattern % 'tavg1_2d_flx_Nx',
            're_pattern': _asset_re_format_pattern.format(name='tavg1_2d_flx_Nx'),
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
        },
        # MERRA2 RAD
        ## SWGDN: Surface incident shortwave flux (W m-2)
        'RAD': {
            'shortname': 'M2T1NXRAD',
            'description': '2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics V5.12.4',
            'url': 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXRAD.5.12.4',
            'pattern': _asset_re_pattern % 'tavg1_2d_rad_Nx',
            're_pattern': _asset_re_format_pattern.format(name='tavg1_2d_rad_Nx'),
            'startdate': datetime.date(1980, 1, 1),
            'latency': 60,
        },
        # MERRA2 ASM
        ## FRLAND: Fraction of land (fraction)
        'ASM': {
            'shortname': 'M2C0NXASM',
            'description': '2d, constants V5.12.4',
            'url': 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980',
            'pattern': _asset_re_pattern % 'const_2d_asm_Nx',
            're_pattern': _asset_re_format_pattern.format(name='const_2d_asm_Nx'),
            'startdate': datetime.date(1980, 1, 1),
            'latency': None,
        }

        #'PROFILE': {
        #     'description': 'Atmospheric Profile',
        #     'pattern': 'MAI6NVANA_PROFILE_*.tif',
        #     'url': 'http://goldsmr5.sci.gsfc.nasa.gov/opendap/MERRA2/M2I6NVANA.5.12.4',
        #     'source': 'MERRA2_%s.inst6_3d_ana_Nv.%04d%02d%02d.nc4',
        #     'startdate': datetime.date(1980, 1, 1),
        #     'latency': 60,
        #     'bandnames': ['0000GMT', '0600GMT', '1200GMT', '1800GMT']
        #},
        # 'PROFILEP': {
        #     'description': 'Atmospheric Profile',
        #     'pattern': 'MAI6NVANA_PROFILE_*.tif',
        #     'url': 'http://goldsmr3.sci.gsfc.nasa.gov:80/opendap/MERRA/MAI6NPANA.5.2.0',
        #     'source': 'MERRA%s.prod.assim.inst6_3d_ana_Np.%04d%02d%02d.hdf',
        #     'startdate': datetime.date(1980, 1, 1),
        #     'latency': 60,
        #     'bandnames': ['0000GMT', '0600GMT', '1200GMT', '1800GMT']
        # },
    }

    def __init__(self, filename):
        """ Inspect a single file and get some metadata """
        super(merraAsset, self).__init__(filename)
        self.sensor = 'merra'
        self.tile = 'h01v01'
        parts = basename(filename).split('.')
        self.asset = parts[1].split('_')[2].upper()
        self._version = int(parts[0].split('_')[1])
        if self.asset == "ASM":
            # assign date of static data sets
            self.date = self._assets[self.asset]['startdate']
        else:
            self.date = datetime.datetime.strptime(parts[2], '%Y%m%d').date()

    @classmethod
    def query_service(cls, asset, tile, date):
        year, month, day = date.timetuple()[:3]
        if asset != "ASM":
            mainurl = "%s/%04d/%02d" % (cls._assets[asset]['url'], year, month)
            pattern = cls._assets[asset]['re_pattern'] % (year, month, day)
        else:
            # asset ASM is for constants which all go into 1980-01-01
            mainurl = cls._assets[asset]['url']
            pattern = cls._assets[asset]['re_pattern'] % (0, 0, 0)
        cpattern = re.compile(pattern)
        with utils.error_handler("Error downloading"):
            # obtain the list of files
            response = cls.Repository.managed_request(mainurl, verbosity=2)
            if response is None:
                return []
        available = []
        for item in response.readlines():
            # inspect the page and extract the full name of the needed file
            if cpattern.search(item):
                if 'xml' in item:
                    continue
                basename = cpattern.findall(item)[0]
                url = '/'.join([mainurl, basename])
                available.append({'basename': basename, 'url': url})
        if len(available) == 0:
            msg = 'Unable to find a remote match for {} at {}'
            utils.verbose_out(msg.format(pattern, mainurl), 4)
        return available

    @classmethod
    def fetch(cls, asset, tile, date):
        """Standard Asset.fetch implementation for downloading assets."""
        if asset == "ASM" and date.date() != cls._assets[asset]['startdate']:
            #TODO: which should it be? if message then remove comment
            #raise Exception, "constants are available for %s only" % cls._assets[asset]['startdate']
            utils.verbose_out('constants are available for %s only' % cls._assets[asset]['startdate'])
            return

        available_assets = cls.query_service(asset, tile, date)
        retrieved_filenames = []

        for asset_info in available_assets:
            basename = asset_info['basename']
            url = asset_info['url']
            outpath = os.path.join(cls.Repository.path('stage'), basename)

            with utils.error_handler("Asset fetch error", continuable=True):
                # obtain the data
                response = cls.Repository.managed_request(url)
                if response is None:
                    continue
                tmp_outpath = tempfile.mkstemp(
                    suffix='.nc4', prefix='downloading',
                    dir=cls.Repository.path('stage')
                )[1]
                with open(tmp_outpath, 'w') as fd:
                    fd.write(response.read())

                # verify that it is a netcdf file
                try:
                    ncroot = Dataset(tmp_outpath)
                    os.rename(tmp_outpath, outpath)
                except Exception as e:
                    text = ''
                    if e.message.endswith('Unknown file format'):
                        token = 'Authorize NASA GESDISC DATA ARCHIVE'
                        html = open(tmp_outpath, 'r').read(100000)
                        if token in html:
                            text = ('\n\nYou need to {t} \nfor your NASA '
                                    'EarthData login.\n').format(t=token)
                    raise Exception(e.message + text)
            retrieved_filenames.append(outpath)
            utils.verbose_out('Retrieved %s' % basename, 2)

        return retrieved_filenames


class merraData(Data):
    """ A tile of data (all assets and products) """
    name = 'merra'
    version = '1.0.0'
    Asset = merraAsset

    _geotransform = (-180.3125, 0.625, 0.0, 90.25, 0.0, -0.5)
    _projection = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

    _products = {

        'tave': {
            'description': 'Ave daily air temperature data (K)',
            'assets': ['SLV'],
            'layers': ['T2M'],
            'bands': ['tave']
        },
        'tmin': {
            'description': 'Min daily air temperature data',
            'assets': ['SLV'],
            'layers': ['T2M'],
            'bands': ['tmin']
        },
        'tmax': {
            'description': 'Max daily air temperature data',
            'assets': ['SLV'],
            'layers': ['T2M'],
            'bands': ['tmax']
        },
        'patm': {
            'description': 'Surface atmospheric pressure (mb)',
            'assets': ['SLV'],
            'layers': ['PS'],
            'bands': ['patm']
        },
        'shum': {
            'description': 'Relative humidity (kg kg-1)',
            'assets': ['SLV'],
            'layers': ['QV2M'],
            'bands': ['shum']
        },
        'rhum': {
            'description': 'Relative humidity (%)',
            'assets': ['SLV'],
            'layers': ['QV2M', 'PS', 'T2M'],
            'bands': ['rhum']
        },
        'prcp': {
            'description': 'Daily total precipitation (mm day-1)',
            'assets': ['FLX'],
            'layers': ['PRECTOT'],
            'bands': ['prcp']
        },
        'wind': {
            'description': 'Daily mean wind speed (m s-1)',
            'assets': ['FLX'],
            'layers': ['SPEED'],
            'bands': ['wind']
        },
        'srad': {
            'description': 'Incident solar radiation (W m-2)',
            'assets': ['RAD'],
            'layers': ['SWGDN'],
            'bands': ['srad']
        },
        'frland': {
            'description': 'Fraction of land (fraction)',
            'assets': ['ASM'],
            'layers': ['FRLAND'],
            'bands': ['frland']
        }
        # TODO:
        # 'temp_modis': {
        #    'description': 'Fraction of land (fraction)',
        #    'assets': ['SLV'],
        #    'layers': ['T2M'],
        #    'bands': ['temp_modis']
        #},
        #'_temps': {
        #    'description': 'Air temperature data',
        #    'assets': ['TS', 'T2M', 'T10M']
        #},
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

    # TODO: move this into Landsat or atmos module to support the wtemp pro
    #@classmethod
    #def lonlat2xy(cls, lon, lat):
    #    """ Convert from lon-lat to x-y in array """
    #    x = int(round((lon - cls._origin[0]) / cls._defaultresolution[0]))
    #    y = int(round((lat - cls._origin[1]) / cls._defaultresolution[1]))
    #    return (x, y)
    #
    #@classmethod
    #def profile(cls, lon, lat, dtime):
    #    """ Retrieve atmospheric profile directly from merra data via OpenDap """
    #    dataset = cls.Asset.opendap_fetch('PROFILE', dtime)
    #    (x, y) = cls.Asset.lonlat2xy(lon, lat)
    #    # TODO - I know these are hours (0, 6, 12, 18), but it's still an assumption
    #    times = [datetime.datetime.combine(dtime.date(), datetime.time(int(d / 60.0))) for d in dataset['time'][:]]
    #    unixtime = time.mktime(dtime.timetuple())
    #    timediff = numpy.array([unixtime - time.mktime(t.timetuple()) for t in times])
    #    timeind = numpy.abs(timediff).argmin()
    #    p = dataset['PS'][timeind, y, x].squeeze()
    #    pthick = dataset['DELP'][timeind, :, y, x].squeeze()[::-1]
    #    pressure = []
    #    for pt in pthick:
    #        pressure.append(p)
    #        p = p - pt
    #    pressure = numpy.array(pressure)
    #    inds = numpy.where(pressure > 0)
    #    data = {
    #        # Pa -> mbar
    #        'pressure': numpy.array(pressure)[inds] / 100.0,
    #        # Kelvin -> Celsius
    #        'temp': dataset['T'][timeind, :, y, x].squeeze()[::-1][inds] - 273.15,
    #        # kg/kg -> g/kg (Mass mixing ratio)
    #        'humidity': dataset['QV'][timeind, :, y, x].squeeze()[::-1][inds] * 1000,
    #        'ozone': dataset['O3'][timeind, :, y, x].squeeze()[::-1][inds] * 1000,
    #    }
    #    return data


    def getlonlat(self):
        """ return the center coordinates of the MERRA tile
            used only by product temp_modis
        """
        hcoord = int(self.id[1:3])
        vcoord = int(self.id[4:])
        lon = -180. + (hcoord - 1) * 12. + 6.
        lat = 90. - vcoord * 10. - 5.
        return lon, lat

    def gmtoffset(self):
        """ return the approximate difference between local time and GMT
            used only by product temp_modis
        """
        lon = self.getlonlat()[0]
        houroffset = lon * (12. / 180.)
        return houroffset

    def write_reduced(self, prod, fun, fout, meta, units):
        """ apply a function to reduce to a daily value """
        assetname = self._products[prod]['assets'][0]
        layername = self._products[prod]['layers'][0]
        bandnames = self._products[prod]['bands']
        assetfile = self.assets[assetname].filename
        ncroot = Dataset(assetfile)
        var = ncroot.variables[layername]
        missing = float(var.missing_value)
        scale = var.scale_factor
        assert scale == 1.0, "Handle non-unity scale functions"
        hourly = np.ma.MaskedArray(var[:])
        hourly.mask = (hourly == missing)
        nb, ny, nx = hourly.shape
        # apply reduce rule
        daily = fun(hourly)
        daily[daily.mask] = missing
        utils.verbose_out('writing %s' % fout, 4)
        imgout = gippy.GeoImage(fout, nx, ny, 1, gippy.GDT_Float32)
        imgout[0].Write(np.array(np.flipud(daily)).astype('float32'))
        imgout.SetBandName(prod, 1)
        imgout.SetUnits(units)
        imgout.SetNoData(missing)
        imgout.SetProjection(self._projection)
        imgout.SetAffine(np.array(self._geotransform))


    @Data.proc_temp_dir_manager
    def process(self, *args, **kwargs):
        """Produce requested products."""
        products = super(merraData, self).process(*args, **kwargs)
        if len(products) == 0:
            return
        bname = os.path.join(self.path, self.basename)
        sensor = "merra"
        for key, val in products.requested.items():
            fout = self.temp_product_filename(sensor, key)
            meta = {}
            VERSION = "1.0"
            meta['VERSION'] = VERSION

            if val[0] == "tave":
                fun = lambda x: x.mean(axis=0) - 273.15
                units = "C"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "tmin":
                fun = lambda x: x.min(axis=0) - 273.15
                units = "C"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "tmax":
                fun = lambda x: x.max(axis=0) - 273.15
                units = "C"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "prcp":
                # conversion from (kg m-2 s-1) to (mm d-1)
                fun = lambda x: x.mean(axis=0)*36.*24.*100.
                units = "mm d-1"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "srad":
                fun = lambda x: x.mean(axis=0)
                units = "W m-2"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "wind":
                fun = lambda x: x.mean(axis=0)
                units = "m s-1"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "shum":
                fun = lambda x: x.mean(axis=0)
                units = "kg kg-1"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "patm":
                fun = lambda x: x.mean(axis=0)/100.
                units = "mb"
                self.write_reduced(val[0], fun, fout, meta, units)

            elif val[0] == "rhum":
                # uses layers QV2M PS T2M from asset SLV
                prod = val[0]
                bandnames = self._products[prod]['bands']
                assetname = self._products[prod]['assets'][0]
                assetfile = self.assets[assetname].filename
                ncroot = Dataset(assetfile)
                qv2m = ncroot.variables['QV2M']
                ps = ncroot.variables['PS']
                t2m = ncroot.variables['T2M']
                missing = float(qv2m.missing_value)
                assert missing == float(ps.missing_value)
                assert missing == float(t2m.missing_value)
                assert qv2m.scale_factor == 1.
                assert ps.scale_factor == 1.
                assert t2m.scale_factor == 1.
                qv2m = np.ma.MaskedArray(qv2m[:])
                ps = np.ma.MaskedArray(ps[:])
                t2m = np.ma.MaskedArray(t2m[:])
                nb, ny, nx = qv2m.shape
                qv2m.mask = (qv2m == missing)
                ps.mask = (ps == missing)
                t2m.mask = (t2m == missing)
                temp = t2m - 273.15
                press = ps/100.
                qair = qv2m
                es = 6.112*np.exp((17.67*temp)/(temp + 243.5))
                e = qair*press/(0.378*qair + 0.622)
                rh = 100. * (e/es)
                rh[rh > 100.] = 100.
                rh[rh < 0.] = 0.
                rhday = rh.mean(axis=0)
                rhday[rhday.mask] = missing
                utils.verbose_out('writing %s' % fout, 4)
                imgout = gippy.GeoImage(fout, nx, ny, 1, gippy.GDT_Float32)
                imgout[0].Write(np.array(np.flipud(rhday)).astype('float32'))
                imgout.SetBandName(val[0], 1)
                imgout.SetUnits('%')
                imgout.SetNoData(missing)
                imgout.SetProjection(self._projection)
                imgout.SetAffine(np.array(self._geotransform))

            elif val[0] == "frland":
                startdate = merraAsset._assets[self._products[val[0]]['assets'][0]]['startdate']
                if self.date != startdate:
                    utils.verbose_out('constants are available for %s only' % startdate)
                    continue
                prod = val[0]
                bandnames = self._products[prod]['bands']
                assetname = self._products[prod]['assets'][0]
                assetfile = self.assets[assetname].filename
                ncroot = Dataset(assetfile)
                var = ncroot.variables['FRLAND']
                missing = float(var.missing_value)
                scale = var.scale_factor
                assert scale == 1.0, "Handle non-unity scale functions"
                frland = np.ma.MaskedArray(var[:])
                frland.mask = (frland == missing)
                nb, ny, nx = frland.shape
                frland = frland.squeeze()
                if frland.mask.sum() > 0:
                    frland[frland.mask] = missing
                utils.verbose_out('writing %s' % fout, 4)
                imgout = gippy.GeoImage(fout, nx, ny, 1, gippy.GDT_Float32)
                imgout[0].Write(np.array(np.flipud(frland)).astype('float32'))
                imgout.SetBandName(prod, 1)
                imgout.SetUnits('fraction')
                imgout.SetNoData(missing)
                imgout.SetProjection(self._projection)
                imgout.SetAffine(np.array(self._geotransform))

            """
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

            elif val[0] == 'profile':
                pass
            """

            # add product to inventory
            archive_fp = self.archive_temp_path(fout)
            self.AddFile(sensor, key, archive_fp)
