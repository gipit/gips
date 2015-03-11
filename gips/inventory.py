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

import sys
import os
from datetime import datetime as dt
import traceback
import numpy
from copy import deepcopy

import gippy
from gips.core import SpatialExtent, TemporalExtent
from gips.tiles import Tiles
from gips.utils import VerboseOut, Colors, basename, mkdir
from gips.data.core import Data
from gips.mapreduce import MapReduce

class Inventory(object):
    """ Base class for inventories """
    _colors = [Colors.PURPLE, Colors.RED, Colors.GREEN, Colors.BLUE]

    def __init__(self):
        pass

    def __getitem__(self, date):
        """ Indexing operator for class """
        return self.data[date]

    def __len__(self):
        """ Length of inventory (# of dates) """
        return len(self.dates)

    def get_subset(self, dates):
        """ Return subset of inventory """
        inv = deepcopy(self)
        for d in inv.dates:
            if d not in dates:
                del inv.data[d]
        return inv

    @property
    def sensor_set(self):
        sset = set()
        for date in self.dates:
            sset.update(self.data[date].sensor_set)
        return sorted(sset)

    @property
    def dates(self):
        """ Get sorted list of dates """
        return sorted(self.data.keys())

    @property
    def numfiles(self):
        """ Total number of files in inventory """
        return sum([len(dat) for dat in self.data.values()])

    @property
    def datestr(self):
        return '%s dates (%s - %s)' % (len(self.dates), self.dates[0], self.dates[-1])

    def color(self, sensor):
        """ Return color for sensor """
        return self._colors[list(self.sensor_set).index(sensor)]

    def pprint(self, md=False):
        """ Print the inventory """
        if len(self.data) == 0:
            print 'No matching files in inventory'
            return
        self.data[self.data.keys()[0]].pprint_asset_header()
        dformat = '%m-%d' if md else '%j'
        oldyear = 0
        formatstr = '{:<12}\n'
        colors = {k: self.color(k) for k in self.sensor_set}
        for date in self.dates:
            # if new year then write out the year
            if date.year != oldyear:
                sys.stdout.write(Colors.BOLD + formatstr.format(date.year) + Colors.OFF)
            self.data[date].pprint(dformat, colors)
            oldyear = date.year
        if self.numfiles != 0:
            VerboseOut("\n\n%s files on %s dates" % (self.numfiles, len(self.dates)), 1)
            self.print_legend()

    def print_legend(self):
        print Colors.BOLD + '\nSENSORS' + Colors.OFF
        for key in sorted(self.sensor_set):
            try:
                desc = self.dataclass.Asset._sensors[key]['description']
                scode = key + ': ' if key != '' else ''
            except:
                desc = ''
                scode = key
            print self.color(key) + '%s%s' % (scode, desc) + Colors.OFF


class ProjectInventory(Inventory):
    """ Inventory of project directory (collection of Data class) """

    def __init__(self, projdir='', products=[]):
        """ Create inventory of a GIPS project directory """
        self.projdir = os.path.abspath(projdir)
        if not os.path.exists(self.projdir):
            raise Exception('Directory %s does not exist!' % self.projdir)

        self.data = {}
        product_set = set()
        sensor_set = set()
        try:
            for dat in Data.discover(self.projdir):
                self.data[dat.date] = dat
                # All products and sensors used across all dates
                product_set = product_set.union(dat.product_set)
                sensor_set = sensor_set.union(dat.sensor_set)

            if not products:
                products = list(product_set)
            self.requested_products = products
            self.sensors = sensor_set
        except:
            VerboseOut(traceback.format_exc(), 4)
            raise Exception("%s does not appear to be a GIPS project directory" % self.projdir)

    def products(self, date):
        """ Intersection of available products and requested products for this date """
        return set(self.data[date].products).intersection(set(self.requested_products))

    def new_image(self, filename, dtype=gippy.GDT_Byte, numbands=1, nodata=None):
        """ Create new image with the same template as the files in project """
        img = gippy.GeoImage(self.data[self.dates[0]].open(self.requested_products[0]))
        imgout = gippy.GeoImage(filename, img, dtype, numbands)
        img = None
        if nodata is not None:
            imgout.SetNoData(nodata)
        return imgout

    def data_size(self):
        """ Get 'shape' of inventory: #products x rows x columns """
        img = gippy.GeoImage(self.data[self.dates[0]].open(self.requested_products[0]))
        sz = (len(self.requested_products), img.YSize(), img.XSize())
        return sz

    def get_data(self, dates=None, products=None, chunk=None):
        """ Read all files as time series, stacking all products """
        # TODO - change to absolute dates
        days = numpy.array([int(d.strftime('%j')) for d in dates])
        imgarr = []
        if dates is None:
            dates = self.dates
        if products is None:
            products = self.requested_products
        for p in products:
            gimg = self.get_timeseries(p, dates=dates)
            # TODO - move numpy.squeeze into swig interface file?
            ch = gippy.Recti(chunk[0], chunk[1], chunk[2], chunk[3])
            arr = numpy.squeeze(gimg.TimeSeries(days.astype('float64'), ch))
            arr[arr == gimg[0].NoDataValue()] = numpy.nan
            if len(days) == 1:
                dims = arr.shape
                arr = arr.reshape(1, dims[0], dims[1])
            imgarr.append(arr)
        data = numpy.vstack(tuple(imgarr))
        return data

    def get_timeseries(self, product='', dates=None):
        """ Read all files as time series """
        if dates is None:
            dates = self.dates
        # TODO - multiple sensors
        filenames = [self.data[date][product] for date in dates]
        img = gippy.GeoImage(filenames)
        return img

    def map_reduce(self, func, numbands=1, products=None, readfunc=None, nchunks=100, **kwargs):
        """ Apply func to inventory to generate an image with numdim output bands """
        if products is None:
            products = self.requested_products
        if readfunc is None:
            readfunc = lambda x: self.get_data(products=products, chunk=x)
        inshape = self.data_size()
        outshape = [numbands, inshape[1], inshape[2]]
        mr = MapReduce(inshape, outshape, readfunc, func, **kwargs)
        mr.run(nchunks=nchunks)
        return mr.assemble()


class DataInventory(Inventory):
    """ Manager class for data inventories (collection of Tiles class) """

    def __init__(self, dataclass, spatial, temporal, products=None, fetch=False, **kwargs):
        """ Create a new inventory
        :dataclass: The Data class to use (e.g., LandsatData, ModisData)
        :spatial: The spatial extent requested
        :temporal: The temporal extent requested
        :products: List of requested products of interest
        :fetch: bool indicated if missing data should be downloaded
        """
        VerboseOut('Retrieving inventory for site %s' % spatial.sitename, 2)

        self.dataclass = dataclass
        Repository = dataclass.Asset.Repository
        self.spatial = spatial
        self.temporal = temporal
        self.products = dataclass.RequestedProducts(products)

        if fetch:
            try:
                dataclass.fetch(self.products.base, self.spatial.tiles, self.temporal)
            except Exception, e:
                raise Exception('Error downloading %s: %s' % (dataclass.name, e))
            dataclass.Asset.archive(Repository.spath())

        # find data
        self.data = {}
        for date in self.temporal.prune_dates(self.spatial.available_dates):
            try:
                dat = Tiles(dataclass, self.spatial, date, self.products, **kwargs)
                if len(dat) > 0:
                    self.data[date] = dat
            except Exception, e:
                raise Exception('DataInventory: Error accessing tiles in %s repository' % dataclass.name)

    @property
    def sensor_set(self):
        """ The set of all sensors used in this inventory """
        return sorted(self.dataclass.Asset._sensors.keys())

    def process(self, *args, **kwargs):
        """ Process assets into requested products """
        # TODO - some check on if any processing was done
        start = dt.now()
        VerboseOut('Processing [%s] on %s dates (%s files)' % (self.products, len(self.dates), self.numfiles), 3)
        if len(self.products.standard) > 0:
            for date in self.dates:
                try:
                    self.data[date].process(*args, **kwargs)
                except:
                    VerboseOut(traceback.format_exc(), 4)
        if len(self.products.composite) > 0:
            self.dataclass.process_composites(self, self.products.composite, **kwargs)
        VerboseOut('Processing completed in %s' % (dt.now() - start), 2)

    def mosaic(self, datadir='./', tree=False, **kwargs):
        """ Create project files for data in inventory """
        # make sure products have been processed first
        self.process(overwrite=False)
        start = dt.now()
        VerboseOut('Creating mosaic project %s' % datadir, 2)
        VerboseOut('  Dates: %s' % self.datestr)
        VerboseOut('  Products: %s' % self.products)

        dout = datadir
        for d in self.dates:
            if tree:
                dout = os.path.join(datadir, d.strftime('%Y%j'))
            self.data[d].mosaic(dout, **kwargs)

        VerboseOut('Completed mosaic project in %s' % (dt.now() - start), 2)

        # Print project inventory when completed
        if not tree:
            inv = ProjectInventory(datadir)
            inv.pprint()
            return inv

    #def warptiles(self):
        """ Just copy or warp all tiles in the inventory """

    def pprint(self, **kwargs):
        """ Print inventory """
        print
        if self.spatial.site is not None:
            print Colors.BOLD + 'Asset Coverage for site %s' % (self.spatial.sitename) + Colors.OFF
            self.spatial.print_tile_coverage()
            print
        else:
            print Colors.BOLD + 'Asset Holdings' + Colors.OFF
        super(DataInventory, self).pprint(**kwargs)
