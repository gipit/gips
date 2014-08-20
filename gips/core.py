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
import tarfile
import traceback
import ftplib
import inspect

import gippy
import gips
from gips.utils import VerboseOut, RemoveFiles, File2List, List2File
from gips.inventory import DataInventory, ProjectInventory

from gips.GeoVector import GeoVector
from gips.version import __version__


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
            return [datetime.strptime(os.path.basename(d), cls._datedir).date() for d in os.listdir(tdir)]
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
            os.makedirs(path)
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
                db = gips.settings.DATABASES['tiles']
                dbstr = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                        (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
                tiles = GeoVector(dbstr, layer=cls._tiles_vector)
                VerboseOut('%s: tiles vector %s' % (cls.__name__, cls._tiles_vector), 4)
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
        self.date = datetime(1900, 1, 1)
        # sensor code (key used in cls.sensors dictionary)
        self.sensor = ''
        # dictionary of existing products in asset {'product name': [filename(s)]}
        self.products = {}

    ##########################################################################
    # Child classes should not generally have to override anything below here
    ##########################################################################
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
        ftpurl = url.split('/')[0]
        ftpdir = url[len(ftpurl):]

        try:
            ftp = ftplib.FTP(ftpurl)
            ftp.login('anonymous', gips.settings.EMAIL)
            pth = os.path.join(ftpdir, date.strftime('%Y'), date.strftime('%j'))
            ftp.set_pasv(True)
            try:
                ftp.cwd(pth)
            except Exception, e:
                raise Exception("Error downloading")

            filenames = []
            ftp.retrlines('LIST', filenames.append)

            for f in ftp.nlst('*'):
                VerboseOut("Downloading %s" % f, 2)
                ftp.retrbinary('RETR %s' % f, open(os.path.join(cls.Repository.spath(), f), "wb").write)

            ftp.close()
        except Exception, e:
            VerboseOut(traceback.format_exc(), 3)
            raise Exception("Error downloading")

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

    # TODO - not sure if this is needed
    @classmethod
    def sensor_names(cls):
        """ All possible sensor names """
        return sorted([s['description'] for s in cls._sensors.values()])

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
                      (numfiles, numlinks, os.path.abspath(path), datetime.now()-start))
        if numfiles != len(fnames):
            VerboseOut('%s files not added to archive' % (len(fnames)-numfiles))
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
            VerboseOut('%s -> quarantine (file error)' % filename, 2)
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
                    except:
                        VerboseOut(traceback.format_exc(), 3)
                        VerboseOut('%s: probem adding to archive' % filename)
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
    """ Collection of assets/products for one tile and date """
    name = 'Data'
    Asset = Asset

    _pattern = '*.tif'
    _products = {}

    def meta(self):
        """ Retrieve metadata for this tile """
        print '%s metadata!' % self.__class__.__name__
        #meta = self.Asset(filename)
        # add metadata to dictionary
        return {}

    def process(self, products, **kwargs):
        """ Make sure all products exist and process if needed """
        pass

    @classmethod
    def process_composites(cls, inventory, products, **kwargs):
        """ Process composite products using provided inventory """
        pass

    def filter(self, **kwargs):
        """ Check if tile passes filter """
        return True

    @classmethod
    def meta_dict(cls):
        return {
            'GIPS Version': __version__,
            'GIPPY Version': gippy.__version__,
        }

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
        self.sensor = ''
        # find all assets
        for asset in self.Asset.discover(tile, date):
            self.assets[asset.asset] = asset
            # sensor and basename assumes same value every time ?
            self.sensor = asset.sensor
            # products that come automatically with assets
            self.products.update(asset.products)
        # find all products
        self.basename = self.id + '_' + self.date.strftime(self.Repository._datedir) + '_' + self.sensor
        prods = self.discover(os.path.join(self.path, self.basename))
        self.products.update(prods)
        if len(self.assets) == 0:
            raise Exception('no assets')
        #VerboseOut('%s %s: assets and products found' % (tile, date), 5)
        #VerboseOut(self.assets, 5)
        #VerboseOut(self.products, 5)

    @property
    def Repository(self):
        return self.Asset.Repository

    def open(self, product=''):
        if len(self.products) == 0:
            raise Exception("No products available to open!")
        if product == '':
            product = self.products.keys()[0]
        fname = self.products[product]
        try:
            return gippy.GeoImage(fname)
        except:
            raise Exception('%s problem reading' % product)

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
        fetched = []
        for a in assets:
            for t in tiles:
                asset_dates = cls.Asset.dates(a, t, dates, days)
                for d in asset_dates:
                    if not cls.Asset.discover(t, d, a):
                        try:
                            status = cls.Asset.fetch(a, t, d)
                            fetched.append((a, t, d))
                        except:
                            pass
        return fetched

    @classmethod
    def products2groups(cls, products):
        """ Convert product list to groupings """
        groups = {}
        for p, val in cls._products.items():
            group = val.get('group', 'Standard')
            groups[group] = {}
        for p, val in products.items():
            group = cls._products[val[0]].get('group', 'Standard')
            groups[group][p] = val
        return groups

    @classmethod
    def arg_parser(cls):
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        groups = {}
        for p in cls._products.values():
            group = p.get('group', 'Standard')
            if group == 'Standard' and p.get('composite'):
                group = 'Composites'
            if group not in groups:
                groups[group] = parser.add_argument_group('%s product arguments' % group)
        for p, product in cls._products.items():
            if p != '':
                group = product.get('group', 'Standard')
                if group == 'Standard' and product.get('composite'):
                    group = 'Composites'
                nargs = product.get('args', None)
                choices = product.get('choices', None)
                if choices is not None:
                    groups[group].add_argument('--%s' % p, help=product['description'], choices=choices, nargs='?', const=[])
                elif nargs == '?':
                    groups[group].add_argument('--%s' % p, help=product['description'], nargs=nargs, const=[])
                elif nargs == '*':
                    groups[group].add_argument('--%s' % p, help=product['description'], nargs=nargs)
                else:
                    groups[group].add_argument('--%s' % p, help=product['description'], action='store_true')
        return parser

    @classmethod
    def extra_arguments(cls):
        return {}

    @classmethod
    def test(cls):
        VerboseOut("%s: running tests" % cls.name)
        # archive
        # inventory
        # process


class Algorithm(object):
    name = 'Algorithm Name'
    __version__ = '0.0.0'

    def __init__(self, **kwargs):
        """ Calls "run" function, or "command" if algorithm uses subparser """
        start = datetime.now()
        if 'projdir' in kwargs:
            self.inv = ProjectInventory(kwargs['projdir'], kwargs.get('products'))
        if 'command' not in kwargs:
            command = 'run'
        else:
            command = kwargs['command']
            VerboseOut('Running %s' % command, 2)
        exec('self.%s(**kwargs)' % command)
        VerboseOut('Completed %s in %s' % (command, datetime.now()-start), 2)

    def run(self, **kwargs):
        pass

    @classmethod
    def info(cls):
        """ Name and versions of algorithm and GIPS library """
        return 'GIPS (v%s): %s (v%s)' % (gips.__version__, cls.name, cls.__version__)

    @classmethod
    def parser(cls):
        """ Parser for algorithm specific options """
        parser = argparse.ArgumentParser(add_help=False)
        return parser

    @classmethod
    def project_parser(cls):
        """ Parser for using GIPS project directory """
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('projdir', help='GIPS Project directory')
        parser.add_argument('-p', '--products', help='Products to operate on', nargs='*')
        return parser

    @classmethod
    def vparser(cls):
        """ Parser for adding verbose keyword """
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2+: debug', default=1, type=int)
        return parser

    @classmethod
    def main(cls):
        """ Main for algorithm classes """
        dhf = argparse.ArgumentDefaultsHelpFormatter

        # Top level parser
        parser = argparse.ArgumentParser(formatter_class=dhf, parents=[cls.parser(), cls.vparser()], description=cls.info())

        args = parser.parse_args()
        gippy.Options.SetVerbose(args.verbose)
        VerboseOut(cls.info())

        try:
            cls(**vars(args))
        except Exception, e:
            VerboseOut('Error in %s: %s' % (cls.name, e))
            VerboseOut(traceback.format_exc(), 3)

"""
    //! Rice detection algorithm
    GeoImage RiceDetect(const GeoImage& image, string filename, vector<int> days, float th0, float th1, int dth0, int dth1) {
        if (Options::Verbose() > 1) cout << "RiceDetect(" << image.Basename() << ") -> " << filename << endl;

        GeoImage imgout(filename, image, GDT_Byte, image.NumBands());
        imgout.SetNoData(0);
        imgout[0].SetDescription("rice");
        for (unsigned int b=1;b<image.NumBands();b++) {
            imgout[b].SetDescription("day"+to_string(days[b]));
        }

        CImg<float> cimg;
        CImg<unsigned char> cimg_datamask, cimg_dmask;
        CImg<int> cimg_th0, cimg_th0_prev, cimg_flood, cimg_flood_prev;
        int delta_day;

        for (unsigned int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
            if (Options::Verbose() > 3) cout << "Chunk " << iChunk << " of " << image[0].NumChunks() << endl;
            cimg = image[0].Read<float>(iChunk);
            cimg_datamask = image[0].DataMask(iChunk);
            CImg<int> DOY(cimg.width(), cimg.height(), 1, 1, 0);
            CImg<int> cimg_rice(cimg.width(), cimg.height(), 1, 1, 0);
            cimg_th0_prev = cimg.get_threshold(th0);
            cimg_flood = (cimg_th0_prev^1).mul(cimg_datamask);

            for (unsigned int b=1;b<image.NumBands();b++) {
                if (Options::Verbose() > 3) cout << "Day " << days[b] << endl;
                delta_day = days[b]-days[b-1];
                cimg = image[b].Read<float>(iChunk);
                cimg_datamask = image[b].DataMask(iChunk);
                cimg_th0 = cimg.get_threshold(th0);
                // Replace any nodata values with the last value
                cimg_forXY(cimg_datamask, x, y) {
                    if (cimg_datamask(x,y) == 0) { cimg_th0(x,y) = cimg_th0_prev(x,y); }
                }

                DOY += delta_day;                                       // running total of days
                DOY.mul(cimg_flood);                                    // reset if it hasn't been flooded yet
                DOY.mul(cimg_th0);                                      // reset if in hydroperiod

                cimg_dmask = DOY.get_threshold(dth1,false,true)^=1;      // mask of where past high date
                DOY.mul(cimg_dmask);

                // locate (and count) where rice criteria met
                CImg<unsigned char> newrice = cimg.threshold(th1,false,true) & DOY.get_threshold(dth0,false,true);
                cimg_rice = cimg_rice + newrice;

                // update flood map
                cimg_flood |= (cimg_th0^1);
                // remove new found rice pixels, and past high date
                cimg_flood.mul(newrice^=1).mul(cimg_dmask);

                //imgout[b].Write(DOY, iChunk);
                imgout[b].Write(DOY, iChunk);
                cimg_th0_prev = cimg_th0;
            }
            imgout[0].Write(cimg_rice,iChunk);              // rice map count
        }
        return imgout;
    }
"""
