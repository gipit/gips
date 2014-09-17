#!/usr/bin/env python

import sys
import os
from datetime import datetime
import traceback
import shutil

import gippy
from gips.utils import VerboseOut, Colors
from gips.GeoVector import GeoVector
import commands
import tempfile


class Tiles(object):
    """ Collection of tiles for a single date and sensor """

    def __init__(self, dataclass, site=None, tiles=None, date=None,
                 products=None, sensors=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage - dictionary of tile id and % coverage with site
        self.tiles - dictionary of tile id and a Tile instance
        """
        self.dataclass = dataclass
        self.site = site
        # Calculate spatial extent
        if isinstance(tiles, list):
            tiles = dict((t, 1) for t in tiles)
        self.tile_coverage = tiles
        #if tiles is not None:
        #    self.tile_coverage = dict((t, 1) for t in tiles)
        #elif site is not None:
        #    self.tile_coverage = self.Repository.vector2tiles(GeoVector(site), **kwargs)
        #else:
        #    self.tile_coverage = dict((t, (1, 1)) for t in self.Repository.find_tiles())
        self.date = date
        self.requested_products = products
        if sensors is None:
            sensors = dataclass.Asset._sensors.keys()

        # For each tile locate files/products
        VerboseOut('%s: searching %s tiles for products and assets' % (self.date, len(self.tile_coverage)), 4)
        self.tiles = {}
        for t in self.tile_coverage.keys():
            try:
                tile = dataclass(t, self.date)
                # Custom filter based on dataclass
                good = tile.filter(**kwargs)
                if good and tile.sensor in sensors:
                    self.tiles[t] = tile
                # TODO - check all tiles - should be same sensor?
                self.sensor = tile.sensor
            except:
                #VerboseOut(traceback.format_exc(), 5)
                continue
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    @property
    def numfiles(self):
        return len(self.tiles)

    def coverage(self):
        """ Calculates % coverage of site for each asset """
        asset_coverage = {}
        for a in self.dataclass.Asset._assets:
            cov = 0.0
            norm = float(len(self.tile_coverage)) if self.site is None else 1.0
            for t in self.tiles:
                if a in self.tiles[t].assets:
                    cov = cov + (self.tile_coverage[t][0] / norm)
            asset_coverage[a] = cov * 100
        return asset_coverage

    def process(self, overwrite=False, **kwargs):
        """ Determines what products need to be processed for each tile and calls Data.process """
        for tileid, tile in self.tiles.items():
            toprocess = {}
            for pname, args in self.requested_products.items():
                if pname not in tile.products or overwrite:
                    toprocess[pname] = args
            if len(toprocess) != 0:
                VerboseOut('Processing products for tile %s: %s' % (tileid, ' '.join(toprocess.keys())), 2)
                self.tiles[tileid].process(toprocess, **kwargs)

    def project(self, datadir, res=None, crop=False, nowarp=False, **kwargs):
        """ Create image of final product (reprojected/mosaiced) """
        res = res if res else self.dataclass.Asset._defaultresolution
        if not hasattr(res, "__len__"):
            res = [res, res]
        start = datetime.now()
        bname = self.date.strftime('%Y%j')
        sensor = self.sensor if self.sensor != '' else ''
        if self.site is None:
            for t in self.tiles:
                for p in self.requested_products:
                    fout = os.path.join(datadir, t + '_' + bname + ('_%s_%s.tif' % (sensor, p)))
                    if not os.path.exists(fout):
                        try:
                            VerboseOut("Creating %s" % os.path.basename(fout))
                            shutil.copy(self.tiles[t].products[p], fout)
                        except Exception:
                            VerboseOut("Problem copying %s" % fout, 2)
                            VerboseOut(traceback.format_exc(), 3)
        else:
            for product in self.requested_products:
                filename = os.path.join(datadir, bname + ('_%s_%s.tif' % (sensor, product)))
                if not os.path.exists(filename):
                    try:
                        filenames = [self.tiles[t].products[product] for t in self.tiles]
                        # TODO - cookiecutter should validate pixels in image.  Throw exception if not
                        if nowarp:
                            self._mosaic(filenames, filename, self.site)
                        else:
                            gippy.CookieCutter(filenames, filename, self.site, res[0], res[1], crop)
                    except:
                        VerboseOut("Problem projecting %s" % filename, 2)
        t = datetime.now() - start
        VerboseOut('%s: created project files for %s tiles in %s' % (self.date, len(self.tiles), t), 2)

    def _mosaic(self, infiles, outfile, vectorfile):
        """ Mosaic multiple files together, but do not warp """
        img = gippy.GeoImage(infiles[0])
        nd = img[0].NoDataValue()
        srs = img.Projection()
        img = None
        for f in range(1, len(infiles)):
            img = gippy.GeoImage(infiles[f])
            if img.Projection() != srs:
                raise Exception("Input files have non-matching projections and must be warped")
        # transform vector to image projection
        vector = GeoVector(vectorfile)
        vsrs = vector.proj()
        from gips.GeoVector import transform_shape
        geom = transform_shape(vector.union(), vsrs, srs)
        extent = geom.bounds
        ullr = "%f %f %f %f" % (extent[0], extent[3], extent[2], extent[1])
        # run command
        nodatastr = '-n %s -a_nodata %s -init %s' % (nd, nd, nd)
        cmd = 'gdal_merge.py -o %s -ul_lr %s %s %s' % (outfile, ullr, nodatastr, " ".join(infiles))
        result = commands.getstatusoutput(cmd)
        imgout = gippy.GeoImage(outfile)
        for b in range(0, img.NumBands()):
            imgout[b].CopyMeta(img[b])
        # warp and rasterize vector
        vec1 = vector.transform(srs)
        vec1name = vec1.filename
        td = tempfile.mkdtemp()
        mask = gippy.GeoImage(os.path.join(td, vector.layer.GetName()), imgout, gippy.GDT_Byte, 1)
        maskname = mask.Filename()
        mask = None
        cmd = 'gdal_rasterize -at -burn 1 -l %s %s %s' % (vec1.layer.GetName(), vec1.filename, maskname)
        result = commands.getstatusoutput(cmd)
        VerboseOut(result, 4)
        mask = gippy.GeoImage(maskname)
        imgout.AddMask(mask[0]).Process().ClearMasks()
        vec1 = None
        mask = None
        shutil.rmtree(os.path.dirname(maskname))
        shutil.rmtree(os.path.dirname(vec1name))
        return imgout

    def pprint_header(self):
        header = Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE')
        for a in sorted(self.dataclass.Asset._assets.keys()):
            header = header + ('{:^10}'.format(a if a != '' else 'Coverage'))
        return header + '{:^10}'.format('Products') + Colors.OFF

    def pprint(self, dformat='%j', color=''):
        """ Print coverage for each and every asset """
        #assets = [a for a in self.dataclass.Asset._assets]
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        asset_coverage = self.coverage()
        for a in sorted(asset_coverage):
            sys.stdout.write(color + '  {:>4.1f}%   '.format(asset_coverage[a]) + Colors.OFF)
        products = [p for t in self.tiles for p in self.tiles[t].products]
        prods = []
        for p in set(products):
            if products.count(p) == len(self.tiles):
                prods.append(p)
        #prods = []
        #for t in self.tiles:
        #    for p in self.tiles[t].products:
                #prods.append(p)
        for p in sorted(set(prods)):
            sys.stdout.write('  ' + p)
        sys.stdout.write('\n')

    #def print_products(self, dformat='%j'):
    #    """ Print products that have been processed """
    #    sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
