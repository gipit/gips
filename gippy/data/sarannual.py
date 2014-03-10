#!/usr/bin/env python

import os
import datetime
import glob
import tarfile
import copy
import numpy
from collections import OrderedDict
from pdb import set_trace

import gippy
from gippy.data.core import Data
from gippy.utils import VerboseOut, File2List, List2File, RemoveFiles


class SARAnnualData(Data):
    """ Represents a single date and temporal extent along with (existing) product variations """
    name = 'SARannual'
    sensors = {
        #'AFBS': 'PALSAR FineBeam Single Polarization',
        'PALSAR': 'PALSAR Mosaic (FineBeam Dual Polarization)',
        #'AWB1': 'PALSAR WideBeam (ScanSAR Short Mode)',
        #'JFBS': 'JERS-1 FineBeam Single Polarization'
    }
    _defaultresolution = [0.00044444444, 0.00044444444]
    _rootdir = '/titan/data/SARannual'
    _tiledir = os.path.join(_rootdir, 'tiles')
    _stagedir = os.path.join(_rootdir, 'stage')
    _datedir = '%Y'

    _tiles_vector = os.path.join(_rootdir, 'vectors', 'tiles.shp')

    _prodpattern = '*'
    _metapattern = '.hdr'

    _assets = {
        'MOS': {
            'pattern': '???????_??_MOS.tar.gz'
        },
        'FNF': {
            'pattern': '???????_??_FNF.tar.gz'
        },        
    }

    _products = OrderedDict([
        ('sign', {
            'description': 'Sigma nought (radar backscatter coefficient)',
            'assets': 'MOS',
        }),
        ('fnf', {
            'description': 'Forest/NonForest Mask',
            'assets': 'FNF',
        })
    ])

    # SAR specific constants
    # launch dates for PALSAR (A) and JERS-1 (J)
    #_launchdate = {'PALSAR': datetime.date(2006, 1, 24)}
    _databands = ["sl_HH", "sl_HV"]

    @classmethod
    def inspect(cls, filename):
        """ Inspect a single file and get some basic info """
        path, fname = os.path.split(filename)
        tile = fname[0:7]
        date = datetime.datetime.strptime(fname[8:10],'%y')

        indexfile = os.path.join(path, fname+'.index')
        if os.path.exists(indexfile):
            datafiles = File2List(indexfile)
        else:
            tfile = tarfile.open(filename)
            try:
                datafiles = tfile.getnames()
            except:
                raise Exception('%s is not a valid tar file' % filename)
            List2File(datafiles, indexfile)

        return {
            'asset': fname[11:14],
            'filename': filename,
            'datafiles': datafiles,
            'tile': tile,
            'date': date,
            'basename': fname[0:10],
            'sensor': 'PALSAR',
            'products': {},
        }

    def meta(self, tile):
        """ Get metadata for this tile """
        return {'CF': -83.0}

    def _extract(self,tile):
        tdata = self.tiles[tile]
        meta = self.meta(tile)
        datafiles = {}
        for a in tdata['assets']:
            if tarfile.is_tarfile(a):
                tfile = tarfile.open(a)
            else: raise Exception('%s is not a valid tar file' % os.path.basename(a))
            tfile.extractall(tdata['path'])
            for f in tdata['datafiles'][a]:
                filename = os.path.join(tdata['path'],f)
                try:
                    os.chmod(filename,0664)
                except: pass
                #if not os.path.exists(fname+'.hdr'):
                if f[-3:] != 'hdr':
                    bandname = f[len(tdata['basename'])+1:]
                    datafiles[bandname] = filename
        return datafiles

    def processtile(self, tile, products):
        """ Make sure all products have been pre-processed """
        if len(products) == 0: raise Exception('Tile %s: No products specified' % tile)
        # extract all data from archive
        tdata = self.tiles[tile]
        datafiles = self._extract(tile)
        meta = self.meta(tile)
        if 'sign' in products.keys():
            bands = [b for b in self._databands if b in datafiles]
            if len(bands) > 0:
                img = gippy.GeoImage(datafiles[bands[0]]) 
                del bands[0]
                for b in bands: img.AddBand(gippy.GeoImage(datafiles[b])[0])
                img.SetNoData(0)
                mask = gippy.GeoImage(datafiles['mask'], False)
                img.AddMask(mask[0] == 255)
                imgout = gippy.SigmaNought(img, products['sign'], meta['CF'])
                tdata['products']['sign'] = imgout.Filename()
                img = None
                imgout = None

        if 'fnf' in products.keys():
            if 'C' in datafiles:
                tdata['products']['fnf'] = datafiles['C']
        #if 'angle' in products.keys():
        #    tdata['products']['linci'] = datafiles['linci']
        # Remove unused stuff
        #for key, f in datafiles.items():
        #    if key not in tdata['products'] and key != 'hdr': RemoveFiles([f],['.hdr','.aux.xml'])

    @classmethod
    def feature2tile(cls,feature):
        """ Get tile designation from a geospatial feature (i.e. a row) """
        fldindex_lat = feature.GetFieldIndex("lat")
        fldindex_lon = feature.GetFieldIndex("lon")
        lat = int(feature.GetField(fldindex_lat)+0.5)
        lon = int(feature.GetField(fldindex_lon)-0.5)
        if lat < 0:
            lat_h = 'S'
        else: lat_h = 'N'
        if lon < 0:
            lon_h = 'W'
        else: lon_h = 'E'
        tile = lat_h + str(abs(lat)).zfill(2) + lon_h + str(abs(lon)).zfill(3)
        return tile

def main(): SARAnnualData.main()
