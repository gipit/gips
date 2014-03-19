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
import sys
import shutil
import errno
import argparse
import ogr
from datetime import datetime
import glob
from shapely.wkb import loads
from shapely.geometry import shape
import tarfile
import traceback
from pdb import set_trace

import gippy
from gippy.utils import VerboseOut, RemoveFiles, File2List, List2File
from gippy.data.datainventory import DataInventory
from gippy.settings import DATABASES


class Repository(object):
    """ Singleton (all classmethods) of file locations and sensor tiling system  """
    # root directory to data
    _rootpath = ''
    # Format code of date directories in repository
    _datedir = '%Y%j'

    _tilesdir = 'tiles'
    _qdir = 'quarantine'
    _vdir = 'vectors'
    _sdir = 'stage'

    _tiles_vector = 'tiles.shp'
    _tile_attribute = 'tile'

    @classmethod
    def feature2tile(cls, feature):
        """ Get tile designation from a geospatial feature (i.e. a row) """
        fldindex = feature.GetFieldIndex(cls._tile_attribute)
        return str(feature.GetField(fldindex))

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    @classmethod
    def path(cls, tile='', date=''):
        path = os.path.join(cls._rootpath, cls._tilesdir)
        if tile != '':
            path = os.path.join(path, tile)
        if date != '':
            path = os.path.join(path, str(date.strftime(cls._datedir)))
        return path

    @classmethod
    def find_tiles(cls):
        """ Get list of all available tiles """
        return os.listdir(os.path.join(cls._rootpath, cls._tilesdir))

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates available in repository for a tile """
        tdir = cls.path(tile=tile)
        if os.path.exists(tdir):
            return [datetime.strptime(os.path.basename(d), cls._datedir).date() for d in os.listdir(tdir)]
        else:
            return []

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    @classmethod
    def qpath(cls):
        """ quarantine path """
        return os.path.join(cls._rootpath, cls._qdir)

    @classmethod
    def vpath(cls):
        """ vectors path """
        return os.path.join(cls._rootpath, cls._vdir)

    @classmethod
    def spath(cls):
        """ staging path """
        return os.path.join(cls._rootpath, cls._sdir)

    @classmethod
    def tiles_vector(cls):
        """ Get GeoVector of sensor grid """
        fname = os.path.join(cls.vpath(), cls._tiles_vector)
        if os.path.isfile(fname):
            VerboseOut('%s: tiles vector %s' % (cls.__name__, fname), 4)
            tiles = gippy.GeoVector(fname)
        else:
            try:
                VerboseOut('%s: tiles vector %s' % (cls.__name__, cls._tiles_vector), 4)
                db = DATABASES['tiles']
                dbstr = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                        (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
                print dbstr
                tiles = gippy.GeoVector(dbstr, layer=cls._tiles_vector)
            except:
                raise Exception('unable to access %s tiles (file or database)' % cls.__name__)
        return tiles

    @classmethod
    def vector2tiles(cls, vector, pcov=0.0, ptile=0.0, **kwargs):
        """ Return matching tiles and coverage % for provided vector """
        start = datetime.now()
        import osr
        geom = vector.union()
        ogrgeom = ogr.CreateGeometryFromWkb(geom.wkb)
        tvector = cls.tiles_vector()
        tlayer = tvector.layer
        trans = osr.CoordinateTransformation(vector.layer.GetSpatialRef(), tlayer.GetSpatialRef())
        ogrgeom.Transform(trans)
        geom = loads(ogrgeom.ExportToWkb())
        tlayer.SetSpatialFilter(ogrgeom)
        tiles = {}
        tlayer.ResetReading()
        feat = tlayer.GetNextFeature()
        #fldindex = feat.GetFieldIndex(cls._tiles_attribute)
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkb())
            area = geom.intersection(tgeom).area
            if area != 0:
                tile = cls.feature2tile(feat)
                tiles[tile] = (area/geom.area, area/tgeom.area)
            feat = tlayer.GetNextFeature()
        remove_tiles = []
        for t in tiles:
            if (tiles[t][0] < pcov/100.0) or (tiles[t][1] < ptile/100.0):
                remove_tiles.append(t)
        for t in remove_tiles:
            tiles.pop(t, None)
        VerboseOut('%s: vector2tiles completed in %s' % (cls.__name__, datetime.now() - start), 4)
        return tiles


class Asset(object):
    """ Class for a single file asset (usually an original raw file or archive) """
    Repository = Repository

    # Sensors
    _sensors = {
        '': {'description': ''},
    }
    # dictionary of assets
    _assets = {
        '': {
            'pattern': '*',
        }
    }

    # TODO - move to be per asset ?
    _defaultresolution = [30.0, 30.0]

    def __init__(self, filename):
        """ Inspect a single file and populate variables. Needs to be extended """
        # full filename to asset
        self.filename = filename
        # the asset code
        self.asset = ''
        # if filename is archive, index of datafiles in archive...needed?
        #self.datafiles = []
        # tile designation
        self.tile = ''
        # full date
        self.date = datetime.now()
         # base/root name of this tile/date (used for product naming)
        self.basename = os.path.basename(filename)
        # sensor code (key used in cls.sensors dictionary)
        self.sensor = ''
        # dictionary of existing products in asset {'product name': [filename(s)]}
        self.products = {}

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    def sensor_meta(self):
        return self._sensors[self.sensor]

    def datafiles(self):
        """ Get list of datafiles from asset (if archive file) """
        if tarfile.is_tarfile(self.filename):
            tfile = tarfile.open(self.filename)
        else:
            raise Exception('%s is not a valid tar file' % self.filename)
        path = os.path.dirname(self.filename)
        indexfile = os.path.join(path, self.filename+'.index')
        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:
            tfile = tarfile.open(self.filename)
            datafiles = tfile.getnames()
            List2File(datafiles, indexfile)
        return datafiles

    def extract(self, filenames=[]):
        """ Extract filenames from asset (if archive file) """
        if tarfile.is_tarfile(self.filename):
            tfile = tarfile.open(self.filename)
        else:
            raise Exception('%s is not a valid tar file' % self.filename)
        path = os.path.dirname(self.filename)
        if len(filenames) == 0:
            filenames = self.datafiles()
        extracted_files = []
        for f in filenames:
            fname = os.path.join(path, f)
            if not os.path.exists(fname):
                tfile.extract(f, path)
            try:
                # this ensures we have permissions on extracted files
                if not os.path.isdir(fname):
                    os.chmod(fname, 0664)
            except:
                pass
            extracted_files.append(fname)
        return extracted_files

    ##########################################################################
    # Class methods
    ##########################################################################
    @classmethod
    def fetch(cls, asset, tile, date):
        """ Get this asset for this tile and date """
        raise Exception("Fetch not implemented for %s" % cls.__name__)

    @classmethod
    def dates(cls, asset, tile, dates, days):
        """ For a given asset get all dates possible (in repo or not) - used for fetch """
        from dateutil.rrule import rrule, DAILY
        # default assumes daily regardless of asset or tile
        datearr = rrule(DAILY, dtstart=dates[0], until=dates[1])
        dates = [dt for dt in datearr if days[0] <= int(dt.strftime('%j')) <= days[1]]
        return dates

    @classmethod
    def discover(cls, tile, date, asset=None):
        """ Factory function returns list of Assets """
        tpath = cls.Repository.path(tile, date)
        if asset is not None:
            assets = [asset]
        else:
            assets = cls._assets.keys()
        found = []
        for a in assets:
            files = glob.glob(os.path.join(tpath, cls._assets[a]['pattern']))
            # more than 1 asset??
            if len(files) > 1:
                VerboseOut(files, 3)
                raise Exception("Duplicate(?) assets found")
            if len(files) == 1:
                found.append(cls(files[0]))
        return found

    # TODO - not sure if this is needed
    @classmethod
    def sensor_names(cls):
        """ All possible sensor names """
        return sorted([s['description'] for s in cls._sensors.values()])

    @classmethod
    def archive(cls, path='', recursive=False, keep=False):
        """ Move assets from directory to archive location """
        start = datetime.now()

        try:
            os.makedirs(cls._qpath)
        except:
            pass

        fnames = []
        if recursive:
            for root, subdirs, files in os.walk(path):
                for a in cls._assets.values():
                    fnames.extend(glob.glob(os.path.join(root, a['pattern'])))
        else:
            for a in cls._assets.values():
                fnames.extend(glob.glob(os.path.join(path, a['pattern'])))

        numlinks = 0
        numfiles = 0
        for f in fnames:
            links = cls._archivefile(f)
            if links >= 0:
                if not keep:
                    RemoveFiles([f], ['.index', '.aux.xml'])
            if links > 0:
                numfiles = numfiles + 1
                numlinks = numlinks + links

        # Summarize
        VerboseOut('%s files (%s links) from %s added to archive in %s' %
                  (numfiles, numlinks, os.path.abspath(path), datetime.now()-start))
        if numfiles != len(fnames):
            VerboseOut('%s files not added to archive' % (len(fnames)-numfiles))

    @classmethod
    def _archivefile(cls, filename):
        """ archive specific file """
        bname = os.path.basename(filename)
        try:
            asset = cls(filename)
        except Exception, e:
            # if problem with inspection, move to quarantine
            qname = os.path.join(cls.Repository.qpath(), bname)
            if not os.path.exists(qname):
                os.link(os.path.abspath(filename), qname)
            VerboseOut('%s -> quarantine (file error)' % filename, 2)
            VerboseOut(traceback.format_exc(), 4)
            return 0

        if not hasattr(asset.date, '__len__'):
            asset.date = [asset.date]
        numlinks = 0
        otherversions = False
        for d in asset.date:
            tpath = cls.Repository.path(asset.tile, d)
            newfilename = os.path.join(tpath, bname)
            if not os.path.exists(newfilename):
                # check if another asset exists
                existing = cls.discover(asset.tile, d, asset.asset)
                if len(existing) > 0:
                    VerboseOut('%s: other version(s) already exists:' % bname, 2)
                    for ef in existing:
                        VerboseOut('\t%s' % os.path.basename(ef), 2)
                    otherversions = True
                else:
                    try:
                        os.makedirs(tpath)
                    except OSError as exc:
                        if exc.errno == errno.EEXIST and os.path.isdir(tpath):
                            pass
                        else:
                            raise Exception('Unable to make data directory %s' % tpath)
                    os.link(os.path.abspath(filename), newfilename)
                    #shutil.move(os.path.abspath(f),newfilename)
                    VerboseOut(bname + ' -> ' + newfilename, 2)
                    numlinks = numlinks + 1
            else:
                VerboseOut('%s already in archive' % filename, 2)
        if otherversions and numlinks == 0:
            return -1
        else:
            return numlinks
        # should return asset instance

    #def __str__(self):
    #    return os.path.basename(self.filename)


class Data(object):
    """ Collection of assets/products for one tile and date """
    name = 'Data'
    Asset = Asset

    _pattern = '*.tif'
    _products = {}
    _groups = {'': ''}

    def meta(self):
        """ Retrieve metadata for this tile """
        print '%s metadata!' % self.__name__
        #meta = self.Asset(filename)
        # add metadata to dictionary
        return {}

    def process(self, products):
        """ Make sure all products exist and process if needed """
        pass

    #def filter(self, **kwargs):
    #    """ Check if tile passes filter """
    #    return True

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    #@property
    #def path(self):
    #    """ Return repository path to this tile dir """
    #    return os.path.join(self.Data._rootpath, self.Data._tilesdir,
    #                        self.id, str(self.date.strftime(self.Data._datedir)))

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    def __init__(self, tile, date):
        """ Find all data and assets for this tile and date """
        self.path = self.Repository.path(tile, date)
        self.id = tile
        self.date = date
        self.assets = {}
        self.products = {}
        self.basename = ''
        # find all assets
        for asset in self.Asset.discover(tile, date):
            self.assets[asset.asset] = asset
            # sensor and basename assumes same value every time ?
            self.sensor = asset.sensor
            # should this be property of tile date, sensor?  same for all assets ?
            self.basename = asset.basename
            # products that come automatically with assets
            self.products.update(asset.products)
        # find all products
        prods = self.discover(os.path.join(self.path, self.basename))
        self.products.update(prods)
        if len(self.assets) == 0:
            raise Exception('no assets')

        #VerboseOut('%s %s: assets and products found' % (tile, date), 5)
        VerboseOut(self.assets, 5)
        VerboseOut(self.products, 5)

    @property
    def Repository(self):
        return self.Asset.Repository

    def open(self, product=''):
        if len(self.products) == 0:
            raise Exception("No products available to open!")
        if product == '':
            product = self.products.keys()[0]
        fname = self.products[product]
        if os.path.exists(fname):
            return gippy.GeoImage(fname)
        else:
            raise Exception('%s product does not exist' % product)

    def link(self, products, path='', copy=False):
        """ Create links in path to tile products """
        for p in products:
            fname = self.products[p]
            bname = os.path.basename(fname)
            fullbname = os.path.join(path, bname)
            if copy:
                try:
                    # TODO - copying doesn't seem to work
                    os.copy(fname, fullbname)
                    VerboseOut('%s: copying' % bname, 2)
                    return
                except:
                    VerboseOut('%s: Problem copying file' % bname, 2)
            # try hard link first, if it fails, soft link
            try:
                os.link(fname, fullbname)
                VerboseOut('%s: hard linking' % bname, 2)
            except:
                try:
                    os.symlink(fname, fullbname)
                    VerboseOut('%s: soft linking' % bname, 2)
                except:
                    VerboseOut('%s: Problem creating link' % bname, 2)

    ##########################################################################
    # Class methods
    ##########################################################################
    @classmethod
    def inventory(cls, **kwargs):
        return DataInventory(cls, **kwargs)

    # TODO - factory function of Tiles ?
    @classmethod
    def discover(cls, basefilename):
        """ Find products in path """
        badexts = ['.hdr', '.xml', 'gz', '.index']
        products = {}
        for p in cls._products:
            files = glob.glob(basefilename+'_'+p+cls._pattern)
            #if len(files) > 0:
            #    products[p] = files
            for f in files:
                rootf = os.path.splitext(f)[0]
                ext = os.path.splitext(f)[1]
                if ext not in badexts:
                    products[rootf[len(basefilename)+1:]] = f
        return products

    @classmethod
    def products2assets(cls, products):
        """ Get list of assets needed for these products """
        assets = []
        for p in products:
            if 'assets' in cls._products[p]:
                assets.extend(cls._products[p]['assets'])
            else:
                assets.append('')
        return set(assets)

    @classmethod
    def fetch(cls, products, tiles, dates, days):
        """ Download data for tiles and add to archive """
        assets = cls.products2assets(products)
        for a in assets:
            for t in tiles:
                asset_dates = cls.Asset.dates(a, t, dates, days)
                for d in asset_dates:
                    if not cls.Asset.discover(t, d, a):
                        status = cls.Asset.fetch(a, t, d)
                        VerboseOut("Fetch status: %s" % status, 2)

    @classmethod
    def products2groups(cls, products):
        """ Convert product list to groupings """
        groups = {}
        for group in cls._groups:
            groups[group] = {}
            for p in cls._groups[group]:
                if p in products:
                    groups[group][p] = products[p]
        return groups

    @classmethod
    def arg_parser(cls):
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        for gname in cls._groups:
            if gname == '':
                gname = 'Standard'
                products = cls._products.keys()
            else:
                products = cls._groups[gname]
            group = parser.add_argument_group('%s product arguments' % gname)
            for p in products:
                prod = cls._products[p]
                group.add_argument('--%s' % p, help=prod['description'], nargs='*')
        return parser
