#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  mhanson@ags.io
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
from shapely.wkb import loads
import tarfile
import traceback
import ftplib

import gippy
from gips import settings, GeoVector, __version__
from gips.utils import VerboseOut, RemoveFiles, File2List, List2File, Colors, basename

"""
The data.core classes are the base classes that are used by individual Data modules.
For a new dataset create children of Repository, Asset, and Data
"""


class Repository(object):
    """ Singleton (all classmethods) of file locations and sensor tiling system  """
    _rootpath = ''
    _tiles_vector = 'tiles.shp'
    _tile_attribute = 'tile'
    # Format code of date directories in repository
    _datedir = '%Y%j'

    _tilesdir = 'tiles'
    _cdir = 'composites'
    _qdir = 'quarantine'
    _sdir = 'stage'
    _vdir = 'vectors'

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
            return sorted([datetime.strptime(os.path.basename(d), cls._datedir).date() for d in os.listdir(tdir)])
        else:
            return []

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
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
    def vpath(cls):
        """ vectors path """
        return cls._path(cls._vdir)

    @classmethod
    def _path(cls, dirname, dirs=''):
        if dirs == '':
            path = os.path.join(cls._rootpath, dirname)
        else:
            path = os.path.join(cls._rootpath, dirname, dirs)
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except Exception:
                raise Exception("Repository: Error making directory %s" % path)
        return path

    @classmethod
    def tiles_vector(cls):
        """ Get GeoVector of sensor grid """
        fname = os.path.join(cls.vpath(), cls._tiles_vector)
        if os.path.isfile(fname):
            tiles = GeoVector(fname)
            VerboseOut('%s: tiles vector %s' % (cls.__name__, fname), 4)
        else:
            try:
                db = settings.DATABASES['tiles']
                dbstr = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                        (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
                tiles = GeoVector(dbstr, layer=cls._tiles_vector)
                VerboseOut('%s: tiles vector %s' % (cls.__name__, cls._tiles_vector), 4)
            except:
                raise Exception('unable to access %s tiles (file or database)' % cls.__name__)
        return tiles

    @classmethod
    def vector2tiles(cls, vector, pcov=0.0, ptile=0.0):
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
        while feat is not None:
            tgeom = loads(feat.GetGeometryRef().ExportToWkb())
            area = geom.intersection(tgeom).area
            if area != 0:
                tile = cls.feature2tile(feat)
                tiles[tile] = (area / geom.area, area / tgeom.area)
            feat = tlayer.GetNextFeature()
        remove_tiles = []
        for t in tiles:
            if (tiles[t][0] < (pcov / 100.0)) or (tiles[t][1] < (ptile / 100.0)):
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
        # if filename is archive, index of datafiles in archive...needed?
        #self.datafiles = []
        # tile designation
        self.tile = ''
        # full date
        self.date = datetime(1900, 1, 1)
        # sensor code (key used in cls.sensors dictionary)
        self.sensor = ''
        # dictionary of existing products in asset {'product name': [filename(s)]}
        self.products = {}

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
    def datafiles(self):
        """ Get list of datafiles from asset (if tar or hdf file) """
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
        except:
            raise Exception('Unable to get datafiles from %s' % self.filename)

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
    def fetch(cls, asset, tile, date):
        """ Get this asset for this tile and date """
        url = cls._assets[asset].get('url', '')
        if url == '':
            raise Exception("%s: URL not defined for asset %s" % (cls.__name__, asset))

        # TODO - extend to be base functionatlity called with super().fetch()
        # TODO - move ftp to separate fetch_ftp function

        ftpurl = url.split('/')[0]
        ftpdir = url[len(ftpurl):]

        try:
            ftp = ftplib.FTP(ftpurl)
            ftp.login('anonymous', settings.EMAIL)
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
            VerboseOut(traceback.format_exc(), 3)
            raise Exception("Error downloading: %s" % e)

    @classmethod
    def available(cls, asset, date):
        """ Check availability of an asset for given date """
        date1 = cls._assets[asset]['startdate']
        date2 = cls._assets[asset]['enddate']
        if date2 is None:
            date2 = datetime.now() - datetime.timedelta(cls._asssets[asset]['latency'])
        if date < date1 or date > date2:
            return False
        return True

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
                VerboseOut(files, 2)
                raise Exception("Duplicate(?) assets found")
            if len(files) == 1:
                found.append(cls(files[0]))
        return found

    @classmethod
    def archive(cls, path='.', recursive=False, keep=False):
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

    #def __str__(self):
    #    return os.path.basename(self.filename)


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

    def process(self, products=None, overwrite=False, **kwargs):
        """ Make sure all products exist and process if needed """
        # Check to see what needs to be processed
        if products is None:
            products = self.RequestedProducts()
        else:
            products = self.RequestedProducts([p for p in products if p not in self.products or overwrite])
        if len(products) > 0:
            VerboseOut("Processing products for tile %s: %s" % (self.id, products), 2)
        return products

    @classmethod
    def process_composites(cls, inventory, products, **kwargs):
        """ Process composite products using provided inventory """
        pass

    def filter(self, **kwargs):
        """ Check if tile passes filter - autofail if there are no assets or products """
        return True

    @classmethod
    def meta_dict(cls):
        return {
            'GIPS Version': __version__,
            'GIPPY Version': gippy.__version__,
        }

    def find_files(self, path):
        """ Search path for valid files """
        filenames = glob.glob(os.path.join(path, self._pattern))
        assetnames = [a.filename for a in self.assets.values()]
        badexts = ['.index', '.xml']
        test = lambda x: x not in assetnames and os.path.splitext(f)[1] not in badexts
        filenames[:] = [f for f in filenames if test(f)]
        return filenames

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

    def available(self, product):
        """ Check availability of product (and required assets) for this date """
        assets = self._products[product]['assets']
        missing = []
        avail = []
        filenames = []
        for asset in assets:
            try:
                fnames = self.assets[asset].datafiles()
            except Exception:
                missing.append(asset)
            else:
                avail.append(asset)
                filenames.extend(fnames)
        if not avail:
            VerboseOut('There are no available assets (%s) on %s for tile %s'
                       % (str(missing), str(self.date), str(self.id), ), 3)
            return None
        return filenames

    def ParseAndAddFiles(self, filenames=None):
        """ Parse and Add filenames to existing filenames """
        if filenames is None:
            filenames = self.find_files(self.path)
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

    def update(self, files):
        """ Merge files dictionary into this instance """
        if files.numfiles > 0:
            self.AddFiles(files.filenames)

    def open(self, product, sensor=None, update=False):
        """ Open and return a GeoImage """
        if sensor is None:
            sensor = self.sensors[product]
        fname = self.filenames[(sensor, product)]
        if os.path.exists(fname):
            return gippy.GeoImage(fname)
        else:
            raise Exception('%s_%s product does not exist' % (sensor, product))

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

    def pprint_header(self):
        """ Print product inventory header """
        return Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE') + '{:^10}'.format('Products') + Colors.OFF

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
    def inventory(cls, **kwargs):
        from gips.inventory import DataInventory
        return DataInventory(cls, **kwargs)

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
        fetched = []
        for a in assets:
            for t in tiles:
                asset_dates = cls.Asset.dates(a, t, dates, days)
                for d in asset_dates:
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

    @classmethod
    def test(cls):
        import random
        VerboseOut(Colors.BOLD + "%s Tests" % cls.name + Colors.OFF)
        #VerboseOut(Colors.BOLD + "Complete Inventory Test" % cls.name + Colors.OFF)
        #inv = cls.inventory()
        # Check for test shapefile
        test_vectors = glob.glob(os.path.join(cls.Asset.Repository.vpath(), 'test*.shp'))
        for testv in test_vectors:
            VerboseOut(Colors.BOLD + "Shapefile Inventory Test using %s" % basename(testv) + Colors.OFF)
            inv = cls.inventory(site=testv)
            # Pick random date
            date = random.choice(inv.dates)
            inv = cls.inventory(site=testv, dates=date.strftime('%Y-%j'))
            # Pick random product from those that have coverage
            pcov = inv[inv.dates[0]].product_coverage()
            for p, cov in pcov.items():
                if cov == 0:
                    del pcov[p]
            products = [random.choice(pcov.keys()) + "-TEST"]
            VerboseOut(Colors.BOLD + "Random product test: %s %s" % (date, ' '.join(products)) + Colors.OFF)
            try:
                inv = cls.inventory(site=testv, dates=date.strftime('%Y-%j'), products=products)
                #inv = cls.inventory(site=testv, dates='1987-282', products=[product])
                inv.pprint()
                print Colors.BOLD + "\nProcessing test %s" % ' '.join(products) + Colors.OFF
                inv.process(overwrite=False)

                print Colors.BOLD + "\nProject test" + Colors.OFF
                inv.project(suffix='TEST')

                print Colors.BOLD + "\nProject test (nomosaic)" + Colors.OFF
                inv.project(suffix='TEST', nomosaic=True)

                print Colors.BOLD + "\nProject test (nowarp)" + Colors.OFF
                inv.project(suffix='TEST-nowarp', nowarp=True)

                # Cleanup
                newfiles = [f for t in inv[inv.dates[0]].tiles for f in inv[inv.dates[0]][t].filenames.values()]
                RemoveFiles(newfiles, ['hdr', 'xml'])
                inv = None
            except:
                VerboseOut(traceback.format_exc(), 3)
                raise Exception("Errors running tests")
