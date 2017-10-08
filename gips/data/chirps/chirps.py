#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    Copyright (C) 2017 Applied Geosolutions
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
import sys
import datetime

from gips.data.core import Repository, Asset, Data
from gips import utils


class chirpsRepository(Repository):
    name = 'CHIRPS'
    description = 'Climate Hazards Group InfraRed Precipitation with Station data'
    _tile_attribute = 'tileid'

# sort of a singleton driver:  one asset, one sensor, one product, set them here
_tile_id = 'global'
_sensor = 'chirps'
_asset_type = 'global-daily'
_product_type = 'precip'

class chirpsAsset(Asset):
    """Asset class for CHIRPS, currently only supports africa daily precipitation."""
    Repository = chirpsRepository

    _sensors = {
        'chirps': {'description': 'Climate Hazards Group InfraRed Precipitation with Station data'}
    }

    _assets = {
        _asset_type: {

            # ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/global_daily/tifs/p05/2010/chirps-v2.0.2010.12.11.tif.gz

            'pattern': (r'^global-daily-chirps-v2\.0\.'
                        r'(?P<year>\d{4})\.(?P<month>\d{2})\.(?P<day>\d{2})\.tif\.gz$'),
            'ftp-basedir': '/pub/org/chg/products/CHIRPS-2.0/global_daily/tifs/p05/',
            'fn-prefix': _asset_type + '-',
            'startdate': datetime.date(1981, 1, 1), # used to prevent impossible searches
            'latency': 21, # latency appears to be approx. 3 weeks
        },
    }

    _host = 'ftp.chg.ucsb.edu' # for ftp asset fetch

    def __init__(self, filename):
        """Inspect a single filename and set some metadata."""
        super(chirpsAsset, self).__init__(filename)
        base_filename = os.path.basename(filename)

        match = re.match(self._assets[_asset_type]['pattern'], base_filename)
        if match is None:
            raise RuntimeError(
                    "File did not match asset naming pattern: '{}'".format(base_filename), filename)

        self.date = datetime.date(*[int(i) for i in match.group('year', 'month', 'day')])
        self.tile = _tile_id
        self.sensor = _sensor
        self.asset = _asset_type

    @classmethod
    def ftp_connect(cls, asset, date):
        """As super, but make the working dir out of (asset, date)."""
        # TODO this is exactly like prismAsset.ftp_connect, DRY out?
        wd = cls._assets[asset]['ftp-basedir'] + str(date.year)
        return super(chirpsAsset, cls).ftp_connect(wd)

    @classmethod
    def query_provider(cls, asset, tile, date):
        """Search for a matching asset in the CHIRPS ftp store.

        Returns (basename, None) on success; (None, None) otherwise."""
        conn = cls.ftp_connect(asset, date)
        filenames = [fn for fn in conn.nlst() if date.strftime('%Y.%m.%d') in fn]
        conn.quit()
        f_cnt = len(filenames)
        if f_cnt == 0:
            return None, None
        if f_cnt > 1:
            utils.verbose_out("Too many assets found:  " + repr(filenames), 1, sys.stderr)
            raise ValueError("Can't decide between {} assets".format(f_cnt), filenames)
        return filenames[0], None # URL is pointless for ftp, have to chdir manually anyway

    @classmethod
    def fetch(cls, asset, tile, date):
        """Fetch to the stage.

        Returns a list with one item, the full path to the staged asset.
        """
        utils.verbose_out('{}: fetch tile {} for {}'.format(asset, tile, date), 3)
        fn_prefix = cls._assets[asset]['fn-prefix']
        stage_dir = cls.Repository.path('stage')
        with utils.error_handler("Error downloading from " + cls._host, continuable=True), \
                utils.make_temp_dir(prefix='fetchtmp', dir=stage_dir) as td_name:
            remote_fn, _ = cls.query_provider(asset, tile, date)
            local_fn = fn_prefix + remote_fn
            temp_fp = os.path.join(td_name, local_fn)
            stage_fp = os.path.join(stage_dir, local_fn)
            utils.verbose_out("Downloading {}, local name {}".format(remote_fn, local_fn), 2)
            conn = cls.ftp_connect(asset, date)
            with open(temp_fp, "wb") as temp_fo:
                conn.retrbinary('RETR ' + remote_fn, temp_fo.write)
            conn.quit()
            os.rename(temp_fp, stage_fp)
            return [stage_fp]
        return []


class chirpsData(Data):
    name = chirpsRepository.name
    version = '0.1.0'
    Asset = chirpsAsset

    _products = {
        _product_type: {
            'description': 'Total rainfall for a period given by the asset',
            'assets': [_asset_type],
            'bands': [{'name': _product_type, 'units': 'mm'}],
        },
    }

    @Data.proc_temp_dir_manager
    def process(self, products=None, overwrite=False, **kwargs):
        """Produce data products and save them to files.

        Only one product; it's processed in the usual way, but for this
        driver, it's extracted from the gzip file and saved (not
        symlink/vsi).  Method signature is largely for campatibilty with
        the rest of gips, eg kwargs is unused.
        """
        needed_products = self.needed_products(products, overwrite)
        if len(needed_products) == 0:
            utils.verbose_out('No new processing required.')
            return

        # sanity check that requested product & asset look ok
        assert (needed_products.requested == {_product_type: [_product_type]}
                and _asset_type in self.assets)

        asset = self.assets[_asset_type]
        err_msg = 'Error creating product {} from {}'.format(
                        _product_type, os.path.basename(asset.filename))
        with utils.error_handler(err_msg, continuable=True):
            temp_fp = self.temp_product_filename(_sensor, _product_type)
            # make gdal/gippy-readable path to the inner file
            vsi_inner_path = '/vsigzip/' + asset.filename
            os.symlink(vsi_inner_path, temp_fp)
            archive_fp = self.archive_temp_path(temp_fp)
            self.AddFile(_sensor, _product_type, archive_fp)
