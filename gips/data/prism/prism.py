#!/usr/bin/env python
################################################################################
#    GIPS PRISM data module
#
#    AUTHOR: Ian Cooke
#    EMAIL:  ircwaves@gmail.com
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
#
#   PRISM Climate data is copyrighted by the PRISM Climate Group,
#   Oregon State University.  For information acceptable use, see
#   http://prism.oregonstate.edu/documents/PRISM_terms_of_use.pdf.
################################################################################

import os
from datetime import datetime, date, timedelta
from csv import DictReader
import re
import ftplib
import tempfile
from gips.core import SpatialExtent, TemporalExtent

from gips.data.core import Repository, Asset, Data

from gips.utils import settings, List2File
from gips import utils

from gippy import GeoImage

from numpy import mean

__author__ = "Ian Cooke <icooke@ags.io>"
__version__ = '0.1.1'


class prismRepository(Repository):
    name = 'PRISM'
    description = 'PRISM Gridded Climate Data'
    _datedir = '%Y%m%d'
    _tile_attribute = 'id'


class prismAsset(Asset):
    Repository = prismRepository
    _sensors = {
        'prism': {'description': 'Daily Gridded Climate Data'}
    }
    _defaultresolution = [4000.0, 4000.0]
    _startdate = date(1981, 1, 1)
    _latency = -7
    # LATENCY (approximate)
    # 6 months for stable
    # 1 month for early
    # 1 week for provisional
    _host = 'prism.nacse.org'
    _assets = {
        '_ppt': {
            'pattern': r'^PRISM_ppt_.+?\.zip$',
            'path': '/daily/ppt',
            'startdate': _startdate,
            'latency': _latency,
        },
        '_tmin': {
            'pattern': r'^PRISM_tmin_.+?\.zip$',
            'path': '/daily/tmin',
            'startdate': _startdate,
            'latency': _latency,
        },
        '_tmax': {
            'pattern': r'^PRISM_tmax_.+?\.zip$',
            'path': '/daily/tmax',
            'startdate': _startdate,
            'latency': _latency,
        },
    }
    _stab_score = {
        'stable': 3,
        'early': 2,
        'provisional': 1,
    }

    @classmethod
    def ftp_connect(cls, asset, date):
        """As super, but make the working dir out of (asset, date)."""
        wd = os.path.join(cls._assets[asset]['path'], date.strftime('%Y'))
        return super(prismAsset, cls).ftp_connect(wd)

    @classmethod
    def query_provider(cls, asset, tile, date, conn=None):
        """Determine availability of data for the given (asset, tile, date).

        Re-use the given ftp connection if possible; Returns (basename,
        None) on success; (None, None) otherwise."""
        if asset not in cls._assets:
            raise ValueError('{} has no defined asset for {}'.format(cls.Repository.name, asset))
        private_conn = conn is None # do we have a shared connection?
        if private_conn:
            conn = cls.ftp_connect(asset, date)
        # get the list of filenames for the year, filter down to the specific date
        filenames = [fn for fn in conn.nlst() if date.strftime('%Y%m%d') in fn]
        if private_conn:
            conn.quit()
        if 0 == len(filenames):
            return None, None
        # choose the one that has the most favorable stability & version values (usually only one)
        return max(filenames, key=(lambda x: prismAsset(x)._version)), None


    @classmethod
    def fetch(cls, asset, tile, date):
        """Fetch to the stage.

        Returns a list with one item, the full path to the staged asset.
        """
        utils.verbose_out('%s: fetch tile %s for %s' % (asset, tile, date), 3)
        with utils.error_handler("Error downloading from " + cls._host, continuable=True):
            ftp = cls.ftp_connect(asset, date) # starts chdir'd to the right directory
            asset_fn, _ = cls.query_provider(asset, tile, date, ftp)
            stage_dir_fp = cls.Repository.path('stage')
            stage_fp = os.path.join(stage_dir_fp, asset_fn)
            with utils.make_temp_dir(prefix='fetchtmp', dir=stage_dir_fp) as td_name:
                temp_fp = os.path.join(td_name, asset_fn)
                utils.verbose_out("Downloading " + asset_fn, 2)
                with open(temp_fp, "wb") as temp_fo:
                    ftp.retrbinary('RETR ' + asset_fn, temp_fo.write)
                ftp.quit()
                os.rename(temp_fp, stage_fp)
            return [stage_fp]
        return []

    def __init__(self, filename):
        """ Inspect a PRISM file """
        super(prismAsset, self).__init__(filename)
        bname = os.path.basename(filename)
        # 'PRISM_<var>_<stability>_<scale&version>_<time period>_bil.zip'
        _, variable, stability, scalever, date, bilzip = bname.split('_')
        assert bilzip == 'bil.zip', "didn't tokenize properly."
        scale = re.sub(r'(.+)D[0-9]+', r'\1', scalever)
        version = re.sub(r'.+(D[0-9]+)', r'\1', scalever)
        self.date = datetime.strptime(date, '%Y%m%d')
        self.asset = '_' + variable
        self.sensor = 'prism'
        self.scale = scale
        self.version = version
        self.stability = stability
        self._version = self._stab_score[self.stability] * .01 + int(self.version[1:])
        # only one tile
        self.tile = 'CONUS'

    def datafiles(self):
        datafiles = super(prismAsset, self).datafiles()
        datafiles = filter(lambda x: x.lower().endswith('.bil'), datafiles)
        if len(datafiles) > 0:
            indexfile = self.filename + '.index'
            utils.verbose_out('indexfile: {}'.format(indexfile), 3)
            List2File(datafiles, indexfile)
            return datafiles


class prismData(Data):
    """ A tile (CONUS State) of PRISM """
    name = 'PRISM'
    version = __version__
    Asset = prismAsset
    # Prism official docs at http://www.prism.oregonstate.edu/FAQ/ say:
    # "Dataset values are stored . . . precipitation as millimeters and
    # temperature as degrees Celsius."
    _products = {
        'ppt': {
            'description': 'Precipitate',
            'assets': ['_ppt'],
            'bands': [{'name': 'ppt', 'units': 'mm'}],
            'startdate': Asset._startdate,
            'latency': Asset._latency,
        },
        'pptsum': {
            'description': 'Cumulative Precipitate',
            'assets': ['_ppt'],
            'bands': [{'name': 'pptsum', 'units': 'mm'}],
            'arguments': ['days: temporal window width (default: 3 days)'],
            'startdate': Asset._startdate,
            'latency': Asset._latency,
        },
        'tmin': {
            'description': 'Daily Minimum Temperature',
            'assets': ['_tmin'],
            'bands': [{'name': 'tmin', 'units': 'degree Celcius'}],
            'startdate': Asset._startdate,
            'latency': Asset._latency,
        },
        'tmax': {
            'description': 'Daily Maximum Temperature',
            'assets': ['_tmax'],
            'bands': [{'name': 'tmin', 'units': 'degree Celcius'}],
            'startdate': Asset._startdate,
            'latency': Asset._latency,
        },
    }

    def process(self, *args, **kwargs):
        """Deduce which products need producing, then produce them."""
        products = super(prismData, self).process(*args, **kwargs)
        if len(products) == 0:
            return
        # overwrite = kwargs.get('overwrite', False)
        # utils.verbose_out('\n\noverwrite = {}\n'.format(overwrite), 2)
        # TODO: overwrite doesn't play well with pptsum -- wonder if it would
        #       if it was made into a composite product (which it is)
        assert len(prismAsset._sensors) == 1  # sanity check to force this code to stay current
        sensor = prismAsset._sensors.keys()[0]

        def get_bil_vsifile(d, a):
            with utils.error_handler('Error accessing asset {}'
                                     .format(d), continuable=True):
                return os.path.join(
                    '/vsizip/' + d.assets[a].filename,
                    d.assets[a].datafiles()[0])

        for key, val in products.requested.items():
            start = datetime.now()
            # check that we have required assets
            requiredassets = self.products2assets([val[0]])
            # val[0] s.b. key w/o product args
            description = self._products['pptsum']['description']
            missingassets = []
            availassets = []
            vsinames = {}

            for asset in requiredassets:
                bil = get_bil_vsifile(self, asset)
                if bil is None:
                    missingassets.append(asset)
                else:
                    availassets.append(asset)
                    vsinames[asset] = os.path.join(
                        '/vsizip/' + self.assets[asset].filename,
                        bil
                    )

            if not availassets:
                utils.verbose_out(
                    'There are no available assets ({}) on {} for tile {}'
                    .format(str(missingassets), str(self.date), str(self.id)),
                    5,
                )
                continue
            prod_fn = '{}_{}_{}.tif'.format(self.basename, 'prism', key)
            archived_fp = os.path.join(self.path, prod_fn) # final destination
            if val[0] in ['ppt', 'tmin', 'tmax']:
                with self.make_temp_proc_dir() as tmp_dir:
                    tmp_fp = os.path.join(tmp_dir, prod_fn)
                    os.symlink(vsinames[self._products[key]['assets'][0]], tmp_fp)
                    os.rename(tmp_fp, archived_fp)
            elif val[0] == 'pptsum':
                if len(val) < 2:
                    lag = 3 # no argument provided, use default lag of 3 days SB configurable.
                    prod_fn = re.sub(r'\.tif$', '-{}.tif'.format(lag), prod_fn)
                    archived_fp = os.path.join(self.path, prod_fn) # have to regenerate, sigh
                    utils.verbose_out('Using default lag of {} days.'.format(lag), 2)
                else:
                    with utils.error_handler("Error for pptsum lag value '{}').".format(val[1])):
                        lag = int(val[1])

                date_spec = '{},{}'.format(
                    datetime.strftime(
                        self.date - timedelta(days=lag), '%Y-%m-%d',
                    ),
                    datetime.strftime(self.date, '%Y-%m-%d'),
                )
                inv = self.inventory(dates=date_spec, products=['ppt'],)
                inv.process()
                # because DataInventory object doesn't update
                inv = self.inventory(dates=date_spec, products=['ppt'],)
                if len(inv.data) < lag:
                    utils.verbose_out(
                        '{}: requires {} preceding days ppt ({} found).'
                        .format(key, lag, len(inv.data)),
                        3,
                    )
                    continue  # go to next product to process
                imgs = []
                for tileobj in inv.data.values():
                    datobj = tileobj.tiles.values()[0]
                    imgs.append(GeoImage(get_bil_vsifile(datobj, '_ppt')))

                with self.make_temp_proc_dir() as tmp_dir:
                    tmp_fp = os.path.join(tmp_dir, prod_fn)
                    oimg = GeoImage(tmp_fp, imgs[0])
                    oimg.SetNoData(-9999)
                    oimg.SetBandName(
                        description + '({} day window)'.format(lag), 1
                    )
                    for chunk in oimg.Chunks():
                        oarr = oimg[0].Read(chunk) * 0.0 # wat
                        for img in imgs:
                            oarr += img[0].Read(chunk)
                        oimg[0].Write(oarr, chunk)
                    oimg.Process()
                    os.rename(tmp_fp, archived_fp)
                oimg = None  # help swig+gdal with GC
                products.requested.pop(key)
            self.AddFile(sensor, key, archived_fp)  # add product to inventory
        return products
