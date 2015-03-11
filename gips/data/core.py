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
import sys
import errno
from osgeo import gdal, ogr
from datetime import datetime
import glob
from itertools import groupby
from shapely.wkt import loads
import tarfile
import traceback
import ftplib
import shutil
import commands

import gippy
from gips import __version__
from gips.GeoVector import GeoVector
from gips.utils import settings, VerboseOut, RemoveFiles, File2List, List2File, Colors, basename, mkdir, open_vector
from gippy.algorithms import CookieCutter

"""
The data.core classes are the base classes that are used by individual Data modules.
For a new dataset create children of Repository, Asset, and Data
"""


def repository_class(clsname):
    """ Get ClassRepository class object """
    exec('from gips.data.%s import %sRepository as cls' % (clsname.lower(), clsname))
    return cls


def asset_class(clsname):
    """ Get ClassAsset class object """
    exec('from gips.data.%s import %sAsset as cls' % (clsname.lower(), clsname))
    return cls


def data_class(clsname):
    """ Get ClassData class object """
    exec('from gips.data.%s import %sData as cls' % (clsname.lower(), clsname))
    return cls


class Repository(object):
    """ Singleton (all classmethods) of file locations and sensor tiling system  """
    # Description of the data source
    description = 'Data source description'
    # Format code of date directories in repository
    _datedir = '%Y%j'

    _tdir = 'tiles'
    _cdir = 'composites'
    _qdir = 'quarantine'
    _sdir = 'stage'
    _vdir = 'vectors'

    @classmethod
    def feature2tile(cls, feature):
        """ Get tile designation from a geospatial feature (i.e. a row) """
        att_name = cls.repo().get('tile_attribute', 'tile')
        fldindex = feature.GetFieldIndex(att_name)
        return str(feature.GetField(fldindex))

    ##########################################################################
    # Override these functions if not using a tile/date directory structure
    ##########################################################################
    @classmethod
    def path(cls, tile='', date=''):
        """ Get absolute data path for this tile and date """
        path = os.path.join(cls.rootpath(), cls._tdir)
        if tile != '':
            path = os.path.join(path, tile)
        if date != '':
            path = os.path.join(path, str(date.strftime(cls._datedir)))
        return path

    @classmethod
    def find_tiles(cls):
        """ Get list of all available tiles """
        return os.listdir(os.path.join(cls.rootpath(), cls._tdir))

    @classmethod
    def find_dates(cls, tile):
        """ Get list of dates available in repository for a tile """
        tdir = cls.path(tile=tile)
        if os.path.exists(tdir):
            return sorted([datetime.strptime(os.path.basename(d), cls._datedir).date() for d in os.listdir(tdir)])
        else:
            return []

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    @classmethod
    def repo(cls):
        """ Get dictionary of repository settings """
        return settings().REPOS[cls.name]

    @classmethod
    def rootpath(cls):
        """ Root path to repository """
        return cls.repo().get('rootpath', '')

    @classmethod
    def cpath(cls, dirs=''):
        """ Composites path """
        return cls._path(cls._cdir, dirs)

    @classmethod
    def qpath(cls):
        """ quarantine path """
        return cls._path(cls._qdir)

    @classmethod
    def spath(cls):
        """ staging path """
        return cls._path(cls._sdir)

    @classmethod
    def _path(cls, dirname, dirs=''):
        """ Get absolute path name to directory, creating if necessary """
        if dirs == '':
            path = os.path.join(cls.rootpath(), dirname)
        else:
            path = os.path.join(cls.rootpath(), dirname, dirs)
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except Exception:
                raise Exception("Repository: Error making directory %s" % path)
        return path

    @classmethod
    def vector(cls):
        """ Get GeoVector of sensor grid """
        # TODO = update to use gippy.GeoVector
        # check location from settings
        vpath = os.path.join('/etc/gips', cls.name)
        vector = open_vector(cls.repo().get('vector', 'tiles.shp'), path=vpath)
        return GeoVector(vector.Filename(), vector.LayerName())

    @classmethod
    def vector2tiles(cls, vector, pcov=0.0, ptile=0.0, tilelist=None):
        """ Return matching tiles and coverage % for provided vector """
        import osr
        # set spatial filter on tiles vector to speedup
        ogrgeom = ogr.CreateGeometryFromWkt(vector.WKT())
        tvector = cls.vector()
        tlayer = tvector.layer
        vsrs = osr.SpatialReference(vector.Projection())
        trans = osr.CoordinateTransformation(vsrs, tlayer.GetSpatialRef())
        ogrgeom.Transform(trans)
        tlayer.SetSpatialFilter(ogrgeom)

        # geometry of desired site (transformed)
        geom = loads(ogrgeom.ExportToWkt())

        tiles = {}
        # step through tiles vector
        tlayer.ResetReading()
        feat = tlayer.GetNextFeature()
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkt())
            area = geom.intersection(tgeom).area
            if area != 0:
                tile = cls.feature2tile(feat)
                tiles[tile] = (area / geom.area, area / tgeom.area)
            feat = tlayer.GetNextFeature()
        remove_tiles = []
        if tilelist is None:
            tilelist = tiles.keys()
        for t in tiles:
            if (tiles[t][0] < (pcov / 100.0)) or (tiles[t][1] < (ptile / 100.0)) or t not in tilelist:
                remove_tiles.append(t)
        for t in remove_tiles:
            tiles.pop(t, None)
        return tiles


class Asset(object):
    """ Class for a single file asset (usually an original raw file or archive) """
    Repository = Repository

    # Sensors
    _sensors = {
        # Does the data have multiple sensors possible for each asset? If not, a single sensor may be fine
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
        # tile designation
        self.tile = ''
        # full date
        self.date = datetime(1858, 4, 6)
        # sensor code (key used in cls.sensors dictionary)
        self.sensor = ''
        # dictionary of existing products in asset {'product name': [filename(s)]}
        self.products = {}

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    def datafiles(self):
        """ Get list of readable datafiles from asset (multiple filenames if tar or hdf file) """
        path = os.path.dirname(self.filename)
        indexfile = os.path.join(path, self.filename + '.index')
        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
            if len(datafiles) > 0:
                return datafiles
        try:
            if tarfile.is_tarfile(self.filename):
                tfile = tarfile.open(self.filename)
                tfile = tarfile.open(self.filename)
                datafiles = tfile.getnames()
            else:
                # Try subdatasets
                fh = gdal.Open(self.filename)
                sds = fh.GetSubDatasets()
                datafiles = [s[0] for s in sds]
            if len(datafiles) > 0:
                List2File(datafiles, indexfile)
                return datafiles
            else:
                return [self.filename]
        except:
            raise Exception('Problem accessing asset(s) in %s' % self.filename)

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
                VerboseOut("Extracting %s" % f, 3)
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
    def discover(cls, tile, date, asset=None):
        """ Factory function returns list of Assets for this tile and date """
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
                VerboseOut(files, 2)
                raise Exception("Duplicate(?) assets found")
            if len(files) == 1:
                found.append(cls(files[0]))
        return found

    @classmethod
    def start_date(cls, asset):
        """ Get starting date for this asset """
        return cls._assets[asset].get('startdate', None)

    @classmethod
    def end_date(cls, asset):
        """ Get ending date for this asset """
        edate = cls._assets[asset].get('enddate', None)
        if edate is None:
            latency = cls._assets[cls.asset].get('latency', None)
            edate = datetime.now() - datetime.timedelta(latency)
        return edate

    @classmethod
    def available(cls, asset, date):
        """ Check availability of an asset for given date """
        date1 = cls._assets[asset].get(['startdate'], None)
        date2 = cls._assets[asset].get(['enddate'], None)
        if date2 is None:
            date2 = datetime.now() - datetime.timedelta(cls._asssets[asset]['latency'])
        if date1 is None or date2 is None:
            return False
        if date < date1 or date > date2:
            return False
        return True

    # TODO - combine this with fetch to get all dates
    @classmethod
    def dates(cls, asset, tile, dates, days):
        """ For a given asset get all dates possible (in repo or not) - used for fetch """
        from dateutil.rrule import rrule, DAILY
        # default assumes daily regardless of asset or tile
        datearr = rrule(DAILY, dtstart=dates[0], until=dates[1])
        dates = [dt for dt in datearr if days[0] <= int(dt.strftime('%j')) <= days[1]]
        return dates

    @classmethod
    def fetch(cls, asset, tile, date):
        """ Fetch stub """
        raise Exception("Fetch not supported for this data source")

    @classmethod
    def fetch_ftp(cls, asset, tile, date):
        """ Fetch via FTP """
        url = cls._assets[asset].get('url', '')
        if url == '':
            raise Exception("%s: URL not defined for asset %s" % (cls.__name__, asset))
        VerboseOut('%s: fetch tile %s for %s' % (asset, tile, date), 3)
        ftpurl = url.split('/')[0]
        ftpdir = url[len(ftpurl):]
        try:
            ftp = ftplib.FTP(ftpurl)
            ftp.login('anonymous', settings().EMAIL)
            pth = os.path.join(ftpdir, date.strftime('%Y'), date.strftime('%j'))
            ftp.set_pasv(True)
            ftp.cwd(pth)

            filenames = []
            ftp.retrlines('LIST', filenames.append)

            for f in ftp.nlst('*'):
                VerboseOut("Downloading %s" % f, 2)
                ftp.retrbinary('RETR %s' % f, open(os.path.join(cls.Repository.spath(), f), "wb").write)
            ftp.close()
        except Exception, e:
            VerboseOut(traceback.format_exc(), 4)
            raise Exception("Error downloading: %s" % e)

    @classmethod
    def archive(cls, path='.', recursive=False, keep=False, **kwargs):
        """ Move assets from directory to archive location """
        start = datetime.now()

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
        assets = []
        for f in fnames:
            archived = cls._archivefile(f)
            if archived[1] >= 0:
                if not keep:
                    RemoveFiles([f], ['.index', '.aux.xml'])
            if archived[1] > 0:
                numfiles = numfiles + 1
                numlinks = numlinks + archived[1]
                assets.append(archived[0])

        # Summarize
        if numfiles > 0:
            VerboseOut('%s files (%s links) from %s added to archive in %s' %
                      (numfiles, numlinks, os.path.abspath(path), datetime.now() - start))
        if numfiles != len(fnames):
            VerboseOut('%s files not added to archive' % (len(fnames) - numfiles))
        return assets

    @classmethod
    def _archivefile(cls, filename):
        """ archive specific file """
        bname = os.path.basename(filename)
        try:
            asset = cls(filename)
        except Exception, e:
            # if problem with inspection, move to quarantine
            VerboseOut(traceback.format_exc(), 3)
            qname = os.path.join(cls.Repository.qpath(), bname)
            if not os.path.exists(qname):
                os.link(os.path.abspath(filename), qname)
            VerboseOut('%s -> quarantine (file error): %s' % (filename, e), 2)
            return (None, 0)

        dates = asset.date
        if not hasattr(dates, '__len__'):
            dates = [dates]
        numlinks = 0
        otherversions = False
        for d in dates:
            tpath = cls.Repository.path(asset.tile, d)
            newfilename = os.path.join(tpath, bname)
            if not os.path.exists(newfilename):
                # check if another asset exists
                existing = cls.discover(asset.tile, d, asset.asset)
                if len(existing) > 0:
                    VerboseOut('%s: other version(s) already exists:' % bname, 1)
                    for ef in existing:
                        VerboseOut('\t%s' % os.path.basename(ef.filename), 1)
                    otherversions = True
                else:
                    try:
                        os.makedirs(tpath)
                    except OSError as exc:
                        if exc.errno == errno.EEXIST and os.path.isdir(tpath):
                            pass
                        else:
                            raise Exception('Unable to make data directory %s' % tpath)
                    try:
                        os.link(os.path.abspath(filename), newfilename)
                        #shutil.move(os.path.abspath(f),newfilename)
                        VerboseOut(bname + ' -> ' + newfilename, 2)
                        numlinks = numlinks + 1
                    except Exception, e:
                        VerboseOut(traceback.format_exc(), 3)
                        raise Exception('Problem adding %s to archive: %s' % (filename, e))
            else:
                VerboseOut('%s already in archive' % filename, 2)
        if otherversions and numlinks == 0:
            return (asset, -1)
        else:
            return (asset, numlinks)
        # should return asset instance


class Data(object):
    """ Collection of assets/products for single date and spatial region """
    name = 'Data'
    version = '0.0.0'
    Asset = Asset

    _pattern = '*.tif'
    _products = {}
    _productgroups = {}

    def meta(self):
        """ Retrieve metadata for this tile """
        return {}

    def process(self, products, overwrite=False, **kwargs):
        """ Make sure all products exist and return those that need processing """
        # TODO - clean up this messy thing
        products = self.RequestedProducts(products)
        products = self.RequestedProducts([p for p in products.products if p not in self.products or overwrite])
        # TODO - this doesnt know that some products aren't available for all dates
        return products

    @classmethod
    def process_composites(cls, inventory, products, **kwargs):
        """ Process composite products using provided inventory """
        pass

    def copy(self, dout, products, site=None, res=None, interpolation=0, crop=False, overwrite=False, tree=False):
        """ Copy products to new directory, warp to projection if given site """
        # TODO - allow hard and soft linking options
        if res is None:
            res = self.Asset._defaultresolution
            #VerboseOut('Using default resolution of %s x %s' % (res[0], res[1]))
        dout = os.path.join(dout, self.id)
        if tree:
            dout = os.path.join(dout, self.date.strftime('%Y%j'))
        mkdir(dout)
        products = self.RequestedProducts(products)
        bname = '%s_%s' % (self.id, self.date.strftime('%Y%j'))
        for p in products.requested:
            if p not in self.sensors:
                # this product is not available for this day
                continue
            sensor = self.sensors[p]
            fin = self.filenames[(sensor, p)]
            fout = os.path.join(dout, "%s_%s_%s.tif" % (bname, sensor, p))
            if not os.path.exists(fout) or overwrite:
                try:
                    if site is not None:
                        # warp just this tile
                        resampler = ['near', 'bilinear', 'cubic']
                        cmd = 'gdalwarp %s %s -t_srs "%s" -tr %s %s -r %s' % \
                               (fin, fout, site.Projection(), res[0], res[1], resampler[interpolation])
                        print cmd
                        #result = commands.getstatusoutput(cmd)
                    else:
                        gippy.GeoImage(fin).Process(fout)
                        #shutil.copyfile(fin, fout)
                except Exception:
                    VerboseOut(traceback.format_exc(), 4)
                    VerboseOut("Problem creating %s" % fout)
        procstr = 'copied' if site is None else 'warped'
        VerboseOut('%s tile %s: %s files %s' % (self.date, self.id, len(products.requested), procstr))

    def filter(self, **kwargs):
        """ Check if tile passes filter - autofail if there are no assets or products """
        return True

    @classmethod
    def meta_dict(cls):
        return {
            'GIPS Version': __version__,
        }

    def find_files(self):
        """ Search path for non-asset files """
        filenames = glob.glob(os.path.join(self.path, self._pattern))
        assetnames = [a.filename for a in self.assets.values()]
        badexts = ['.index', '.xml']
        test = lambda x: x not in assetnames and os.path.splitext(f)[1] not in badexts
        filenames[:] = [f for f in filenames if test(f)]
        return filenames

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    def __init__(self, tile=None, date=None, path=None):
        """ Find all data and assets for this tile and date """
        self.id = tile
        self.date = date
        self.path = ''
        self.basename = ''
        self.assets = {}                # dict of asset name: Asset instance
        self.filenames = {}             # dict of (sensor, product): filename
        self.sensors = {}               # dict of asset/product: sensor
        if tile is not None and date is not None:
            self.path = self.Repository.path(tile, date)
            self.basename = self.id + '_' + self.date.strftime(self.Repository._datedir)
            # find all assets
            for asset in self.Asset.discover(tile, date):
                self.assets[asset.asset] = asset
                # products that come automatically with assets
                for p, val in asset.products.items():
                    self.filenames[(asset.sensor, p)] = val
                    self.sensors[p] = asset.sensor
                self.filenames.update({(asset.sensor, p): val for p, val in asset.products.items()})
                self.sensors[asset.asset] = asset.sensor
            # Find products
            self.ParseAndAddFiles()
        elif path is not None:
            self.path = path

    @property
    def Repository(self):
        """ The repository for this class """
        return self.Asset.Repository

    @classmethod
    def RequestedProducts(cls, *args, **kwargs):
        from gips.core import RequestedProducts
        return RequestedProducts(cls, *args, **kwargs)

    def __getitem__(self, key):
        """ Get filename for product key """
        if type(key) == tuple:
            return self.filenames[key]
        else:
            return self.filenames[(self.sensor_set[0], key)]

    def __str__(self):
        """ Text representation """
        return '%s: %s: %s' % (self.name, self.date, ' '.join(self.product_set))

    def __len__(self):
        """ Number of products """
        return len(self.filenames)

    @property
    def valid(self):
        return False if len(self.filenames) == 0 and len(self.assets) == 0 else True

    @property
    def day(self):
        return self.date.strftime('%j')

    @property
    def sensor_set(self):
        """ Return list of sensors used """
        return list(set(sorted(self.sensors.values())))

    @property
    def products(self):
        """ Get list of products """
        return sorted([k[1] for k in self.filenames.keys()])

    @property
    def product_set(self):
        """ Return list of products available """
        return list(set(self.products))

    def ParseAndAddFiles(self, filenames=None):
        """ Parse and Add filenames to existing filenames """
        if filenames is None:
            filenames = self.find_files()
        datedir = self.Repository._datedir
        for f in filenames:
            bname = basename(f)
            parts = bname.split('_')
            if len(parts) < 3 or len(parts) > 4:
                # Skip this file
                VerboseOut('Unrecognizable file: %s' % f, 3)
                continue
            offset = 1 if len(parts) == 4 else 0
            try:
                if self.date is None:
                    # First time through
                    self.date = datetime.strptime(parts[0 + offset], datedir).date()
                else:
                    date = datetime.strptime(parts[0 + offset], datedir).date()
                    if date != self.date:
                        raise Exception('Mismatched dates: %s' % ' '.join(filenames))
                sensor = parts[1 + offset]
                product = parts[2 + offset]
                self.AddFile(sensor, product, f)
            except Exception:
                # This was just a bad file
                VerboseOut('Unrecognizable file: %s' % f, 3)
                continue

    def AddFile(self, sensor, product, filename):
        """ Add products (dictionary  (sensor, product): filename) to instance """
        self.filenames[(sensor, product)] = filename
        # TODO - currently assumes single sensor for each product
        self.sensors[product] = sensor

    def asset_filenames(self, product):
        assets = self._products[product]['assets']
        filenames = []
        for asset in assets:
            filenames.extend(self.assets[asset].datafiles())
        if len(filenames) == 0:
            VerboseOut('There are no available assets on %s for tile %s' % (str(self.date), str(self.id), ), 3)
            return None
        return filenames

    def open(self, product, sensor=None, update=False):
        """ Open and return a GeoImage """
        if sensor is None:
            sensor = self.sensors[product]
        try:
            fname = self.filenames[(sensor, product)]
            return gippy.GeoImage(fname)
        except:
            raise Exception('error reading product (%s, %s)' % (sensor, product))

    def open_assets(self, product):
        """ Open and return a GeoImage of the assets """
        return gippy.GeoImage(self.asset_filenames(product))

    # TODO - make general product_filter function
    def masks(self, patterns=None):
        """ List all products that are masks """
        if not patterns:
            patterns = ['acca', 'fmask', 'mask']
        m = []
        for p in self.products:
            if any(pattern in p for pattern in patterns):
                m.append(p)
        return m

    @classmethod
    def pprint_header(cls):
        """ Print product inventory header showing product coverage"""
        header = Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE')
        for a in sorted(cls._products.keys()):
            header = header + ('{:^10}'.format(a if a != '' else 'Coverage'))
        return header + '{:^10}'.format('Product') + Colors.OFF

    @classmethod
    def pprint_asset_header(cls):
        """ Print header info for asset coverage """
        header = Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE')
        for a in sorted(cls.Asset._assets.keys()):
            header = header + ('{:^10}'.format(a if a != '' else 'Coverage'))
        header = header + '{:^10}'.format('Product') + Colors.OFF
        print header

    def pprint(self, dformat='%j', colors=None):
        """ Print product inventory for this date """
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        if colors is None:
            sys.stdout.write('  '.join(sorted(self.products)))
        else:
            for p in sorted(self.products):
                sys.stdout.write(colors[self.sensors[p]] + p + Colors.OFF + '  ')
        sys.stdout.write('\n')

    ##########################################################################
    # Class methods
    ##########################################################################
    @classmethod
    def discover(cls, path):
        """ Find products in path and return Data object for each date """
        files0 = glob.glob(os.path.join(path, cls._pattern))
        files = []
        datedir = cls.Asset.Repository._datedir
        for f in files0:
            parts = basename(f).split('_')
            if len(parts) == 3 or len(parts) == 4:
                try:
                    datetime.strptime(parts[len(parts) - 3], datedir)
                    files.append(f)
                except:
                    pass

        datas = []
        if len(files) == 0:
            return datas

        # Group by date
        sind = len(basename(files[0]).split('_')) - 3

        func = lambda x: datetime.strptime(basename(x).split('_')[sind], datedir).date()
        for date, fnames in groupby(sorted(files), func):
            dat = cls(path=path)
            dat.ParseAndAddFiles(list(fnames))
            datas.append(dat)
        return datas

    @classmethod
    def inventory(cls, site=None, key='', where=[], tiles=None, pcov=0.0, ptile=0.0, 
                  dates=None, days=None, **kwargs):
        """ Return list of inventories (size 1 if not looping through geometries) """
        from gips.inventory import DataInventory
        from gips.core import SpatialExtent, TemporalExtent
        spatial = SpatialExtent.factory(cls, site=site, key=key, where=where, tiles=tiles, pcov=pcov, ptile=ptile)
        temporal = TemporalExtent(dates, days)
        return DataInventory(cls, spatial[0], temporal, **kwargs)

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
    def fetch(cls, products, tiles, textent):
        """ Download data for tiles and add to archive """
        assets = cls.products2assets(products)
        fetched = []
        for a in assets:
            for t in tiles:
                asset_dates = cls.Asset.dates(a, t, textent.datebounds, textent.daybounds)
                for d in asset_dates:
                    # if we don't have it already
                    if not cls.Asset.discover(t, d, a):
                        try:
                            cls.Asset.fetch(a, t, d)
                            fetched.append((a, t, d))
                        except Exception, e:
                            VerboseOut(traceback.format_exc(), 4)
                            VerboseOut('Problem fetching asset: %s' % e, 3)
        return fetched

    @classmethod
    def product_groups(cls):
        """ Return dict of groups and products in each one """
        groups = cls._productgroups
        groups['Standard'] = []
        grouped_products = [x for sublist in cls._productgroups.values() for x in sublist]
        for p in cls._products:
            if p not in grouped_products:
                groups['Standard'].append(p)
        if len(groups['Standard']) == 0:
            del groups['Standard']
        return groups

    @classmethod
    def products2groups(cls, products):
        """ Convert product list to groupings """
        p2g = {}
        groups = {}
        allgroups = cls.product_groups()
        for g in allgroups:
            groups[g] = {}
            for p in allgroups[g]:
                p2g[p] = g
        for p, val in products.items():
            g = p2g[val[0]]
            groups[g][p] = val
        return groups

    @classmethod
    def print_products(cls):
        print Colors.BOLD + "\n%s Products v%s" % (cls.name, cls.version) + Colors.OFF
        groups = cls.product_groups()
        opts = False
        txt = ""
        for group in groups:
            txt = txt + Colors.BOLD + '\n%s Products\n' % group + Colors.OFF
            for p in sorted(groups[group]):
                h = cls._products[p]['description']
                txt = txt + '   {:<12}{:<40}\n'.format(p, h)
                if 'arguments' in cls._products[p]:
                    opts = True
                    #sys.stdout.write('{:>12}'.format('options'))
                    args = [['', a] for a in cls._products[p]['arguments']]
                    for a in args:
                        txt = txt + '{:>12}     {:<40}\n'.format(a[0], a[1])
        if opts:
            print "  Optional qualifiers listed below each product."
            print "  Specify by appending '-option' to product (e.g., ref-toa)"
        sys.stdout.write(txt)

    @classmethod
    def extra_arguments(cls):
        return {}
