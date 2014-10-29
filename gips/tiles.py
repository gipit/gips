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

    def __init__(self, dataclass, site=None, tiles=None, date=None, products=None, sensors=None, **kwargs):
        """ Locate data matching vector location (or tiles) and date
        self.tile_coverage      dict of tile id: %coverage with site
        self.tiles              dict of tile id: tile instance
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
                good = tile.filter(**kwargs)
                # TODO - fix custom filter based on dataclass...sensor should pass to 'dataclass'
                if good:  # and tile.sensor in sensors:
                    self.tiles[t] = tile
                # assume same
            except:
                #VerboseOut(traceback.format_exc(), 5)
                continue
        if len(self.tiles) == 0:
            raise Exception('No valid data found')

    @property
    def sensor_set(self):
        """ Return list of sensors used in all the tiles """
        s = set()
        for t in self.tiles:
            s.update(self.tiles[t].sensor_set)
        return list(s)

    def which_sensor(self, key):
        """ Get sensor code used for provided asset or product key """
        for t in self.tiles:
            if key in self.tiles[t].sensors:
                return self.tiles[t].sensors[key]

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

    def project(self, datadir, res, crop=False, nowarp=False, nomosaic=False, **kwargs):
        """ Create image of final product (reprojected/mosaiced) """
        start = datetime.now()
        bname = self.date.strftime('%Y%j')

        if self.site is None:
            nomosaic = True
            nowarp = True
        if not hasattr(res, "__len__"):
            res = [res, res]

        if nomosaic:
            # Tile project
            # TODO - allow hard and soft link options
            for t in self.tiles:
                tiledir = datadir.replace('TILEID', t)
                if not os.path.exists(tiledir):
                    os.makedirs(tiledir)
                for p in self.requested_products:
                    sensor = self.which_sensor(p)
                    filename = self.tiles[t].products[p]
                    fout = os.path.join(tiledir, t + '_' + bname + ('_%s_%s.tif' % (sensor, p)))
                    if not os.path.exists(fout):
                        try:
                            VerboseOut("Creating %s" % os.path.basename(fout))
                            if nowarp:
                                gippy.GeoImage(filename).Process(fout)
                            else:
                                # Warp each tile
                                gippy.CookieCutter([filename], fout, self.site, res[0], res[1], crop)
                        except Exception:
                            VerboseOut("Problem creating %s" % fout, 2)
                            VerboseOut(traceback.format_exc(), 3)
        else:
            # Shapefile project
            if not os.path.exists(datadir):
                os.makedirs(datadir)
            for product in self.requested_products:
                sensor = self.which_sensor(product)
                fout = os.path.join(datadir, bname + ('_%s_%s.tif' % (sensor, product)))
                if not os.path.exists(fout):
                    try:
                        filenames = [self.tiles[t].products[product] for t in self.tiles]
                        # TODO - cookiecutter should validate pixels in image.  Throw exception if not
                        if nowarp:
                            self._mosaic(filenames, fout, self.site)
                        else:
                            gippy.CookieCutter(filenames, fout, self.site, res[0], res[1], crop)
                    except:
                        VerboseOut("Problem projecting %s" % fout, 2)
                        VerboseOut(traceback.format_exc(), 3)
        t = datetime.now() - start
        VerboseOut('%s: created project files for %s tiles in %s' % (self.date, len(self.tiles), t), 2)

    def _mosaic(self, infiles, outfile, vectorfile):
        """ Mosaic multiple files together, but do not warp """
        img = gippy.GeoImage(infiles[0])
        nd = img[0].NoDataValue()
        srs = img.Projection()
        for f in range(1, len(infiles)):
            _img = gippy.GeoImage(infiles[f])
            if _img.Projection() != srs:
                raise Exception("Input files have non-matching projections and must be warped")
            _img = None
        # transform vector to image projection
        vector = GeoVector(vectorfile)
        vsrs = vector.proj()
        from gips.GeoVector import transform_shape
        geom = transform_shape(vector.union(), vsrs, srs)
        extent = geom.bounds
        ullr = "%f %f %f %f" % (extent[0], extent[3], extent[2], extent[1])

        # run merge command
        nodatastr = '-n %s -a_nodata %s -init %s' % (nd, nd, nd)
        cmd = 'gdal_merge.py -o %s -ul_lr %s %s %s' % (outfile, ullr, nodatastr, " ".join(infiles))
        result = commands.getstatusoutput(cmd)
        VerboseOut('%s: %s' % (cmd, result), 4)
        imgout = gippy.GeoImage(outfile, True)
        for b in range(0, img.NumBands()):
            imgout[b].CopyMeta(img[b])
        img = None
        return self.crop2vector(imgout, vector)

    def crop2vector(self, img, vector):
        """ Crop a GeoImage down to a vector """
        start = datetime.now()
        # transform vector to srs of image
        srs = img.Projection()
        vec_t = vector.transform(srs)
        vecname = vec_t.filename
        # rasterize the vector
        td = tempfile.mkdtemp()
        mask = gippy.GeoImage(os.path.join(td, vector.layer.GetName()), img, gippy.GDT_Byte, 1)
        maskname = mask.Filename()
        mask = None
        cmd = 'gdal_rasterize -at -burn 1 -l %s %s %s' % (vec_t.layer.GetName(), vecname, maskname)
        result = commands.getstatusoutput(cmd)
        VerboseOut('%s: %s' % (cmd, result), 4)
        mask = gippy.GeoImage(maskname)
        img.AddMask(mask[0]).Process().ClearMasks()
        vec_t = None
        mask = None
        shutil.rmtree(os.path.dirname(maskname))
        shutil.rmtree(os.path.dirname(vecname))
        VerboseOut('Cropped to vector in %s' % (datetime.now() - start))
        return img

    def pprint_header(self):
        header = Colors.BOLD + Colors.UNDER + '{:^12}'.format('DATE')
        for a in sorted(self.dataclass.Asset._assets.keys()):
            header = header + ('{:^10}'.format(a if a != '' else 'Coverage'))
        return header + '{:^10}'.format('Products') + Colors.OFF

    def pprint(self, dformat='%j', colors=None):
        """ Print coverage for each and every asset """
        #assets = [a for a in self.dataclass.Asset._assets]
        sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
        asset_coverage = self.coverage()
        for a in sorted(asset_coverage):
            color = ['', '']
            if colors is not None:
                s = self.which_sensor(a)
                if s is not None:
                    color = [colors[s], Colors.OFF]
            cov = asset_coverage[a]
            if cov > 0:
                sys.stdout.write(color[0] + '  {:>4.1f}%   '.format(cov) + color[1])
            else:
                sys.stdout.write('          ')
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
            color = colors[self.which_sensor(p)]
            sys.stdout.write('  ' + color + p + Colors.OFF)
        sys.stdout.write('\n')

    #def print_products(self, dformat='%j'):
    #    """ Print products that have been processed """
    #    sys.stdout.write('{:^12}'.format(self.date.strftime(dformat)))
