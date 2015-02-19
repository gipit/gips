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
import errno
import gippy
from datetime import datetime
import tempfile
import commands
import shutil
import traceback


class Colors():
    _c = '\033['
    OFF     = _c + '0m'
    # Font styles
    BOLD    = _c + '1m'
    UNDER   = _c + '4m'
    REV     = _c + '7m'
    # Text colors
    BLACK   = _c + '30m'
    RED     = _c + '31m'
    GREEN   = _c + '32m'
    YELLOW  = _c + '33m'
    BLUE    = _c + '34m'
    PURPLE  = _c + '35m'
    CYAN    = _c + '36m'
    WHITE   = _c + '37m'
    # Background colors
    _BLACK  = _c + '40m'
    _RED    = _c + '41m'
    _GREEN  = _c + '42m'
    _YELLOW = _c + '43m'
    _BLUE   = _c + '44m'
    _PURPLE = _c + '45m'
    _CYAN   = _c + '46m'
    _WHITE  = _c + '47m'


def VerboseOut(obj, level=1):
    if gippy.Options.Verbose() >= level:
        #pprint.PrettyPrinter().pprint(obj)
        if not isinstance(obj, (list, tuple)):
            obj = [obj]
        for o in obj:
            print o


def File2List(filename):
    f = open(filename)
    txt = f.readlines()
    txt2 = []
    for t in txt:
        txt2.append(t.rstrip('\n'))
    return txt2


def List2File(lst, filename):
    f = open(filename, 'w')
    f.write('\n'.join(lst) + '\n')
    f.close()


def RemoveFiles(filenames, extensions=['']):
    for f in filenames:
        for ext in ([''] + extensions):
            try:
                os.remove(f + ext)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise
                continue


def basename(str):
    return os.path.splitext(os.path.basename(str))[0]


def mkdir(dname):
    """ Create a directory if it doesn't exist """
    if not os.path.exists(dname):
        os.makedirs(dname)
    return dname


def link(src, dst):
    """ Create link in this directory """
    if os.path.lexists(dst):
        os.remove(dst)
    # link path path relative to dst
    os.symlink(os.path.relpath(src, os.path.dirname(dst)), os.path.abspath(dst))
    return dst


def settings():
    """ Retrieve GIPS settings """
    import imp
    try:
        import gips.settings
        return gips.settings
    except:
        import imp
        return imp.load_source('gips.settings', '/etc/gips/settings.py')


import time
from functools import wraps


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print "%s time: %s seconds" % (function.func_name, str(t1 - t0))
        return result
    return function_timer


def parse_vectorname(fname, path=''):
    """ Parse name to determine if database or filename and return (shortname, filename, layer, feature) """
    parts = fname.split(':')
    shortname = basename(parts[0]).replace('_', '').replace(':', '-')
    if len(parts) == 1:
        return (shortname, os.path.join(path, fname), '', 0)
    if parts[0] in settings().DATABASES.keys():
        try:
            db = settings().DATABASES[parts[0]]
            filename = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                        (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
            layer = parts[1]
            feature = 0 if len(parts) < 3 else parts[2]
            shortname = '%s-%s-%s' % (shortname, layer.replace('_', ''), feature)
            return (shortname, filename, layer, int(feature))
        except Exception, e:
            VerboseOut(traceback.format_exc(), 4)
            VerboseOut('Error accessing database vector %s: %s' % (fname, e))
    else:
        feature = 0 if len(parts) < 2 else parts[1]
        shortname = '%s-%s' % (shortname, feature)
        return (shortname, os.path.join(path, parts[0]), '', int(feature))


from shapely.wkt import loads as wktloads
from osr import SpatialReference, CoordinateTransformation
from ogr import CreateGeometryFromWkt

def transform_shape(shape, ssrs, tsrs):
    """ Transform shape from ssrs to tsrs (all wkt) and return as wkt """
    ogrgeom = CreateGeometryFromWkt(shape)
    trans = CoordinateTransformation(SpatialReference(ssrs), SpatialReference(tsrs))
    ogrgeom.Transform(trans)
    wkt = ogrgeom.ExportToWkt()
    ogrgeom = None
    return wkt


def crop2vector(img, vector):
    """ Crop a GeoImage down to a vector - only used by mosaic """
    # TODO - update to use gippy.GeoVector
    # TODO - incorporate into GIPPY?
    #start = datetime.now()
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
    #VerboseOut('Cropped to vector in %s' % (datetime.now() - start), 3)
    return img


def mosaic(images, outfile, vector):
    """ Mosaic multiple files together, but do not warp """
    nd = images[0][0].NoDataValue()
    srs = images[0].Projection()
    # check they all have same projection
    filenames = [images[0].Filename()]
    for f in range(1, images.NumImages()):
        if images[f].Projection() != srs:
            raise Exception("Input files have non-matching projections and must be warped")
        filenames.append(images[f].Filename())
    # transform vector to image projection
    geom = wktloads(transform_shape(vector.WKT(), vector.Projection(), srs))

    extent = geom.bounds
    ullr = "%f %f %f %f" % (extent[0], extent[3], extent[2], extent[1])

    # run merge command
    nodatastr = '-n %s -a_nodata %s -init %s' % (nd, nd, nd)
    cmd = 'gdal_merge.py -o %s -ul_lr %s %s %s' % (outfile, ullr, nodatastr, " ".join(filenames))
    result = commands.getstatusoutput(cmd)
    VerboseOut('%s: %s' % (cmd, result), 4)
    imgout = gippy.GeoImage(outfile, True)
    for b in range(0, images[0].NumBands()):
        imgout[b].CopyMeta(images[0][b])
    img = None
    return imgout
    # TODO - update to use gippy.GeoVector
    #return crop2vector(imgout, vector)


# old code utilizing shared memory array
# Chunk it up
#chunksz = int(data.shape[0] / nproc)
#extra = data.shape[0] - chunksz * nproc
#chunks = [chunksz] * (nproc - extra) + [chunksz + 1] * extra

#queue = multiprocessing.Queue()
#from agspy.contrib import shmarray
#classmap = shmarray.create_copy(classmap)
#tmp = numpy.ctypeslib.as_ctypes(classmap)
#cmap = sharedctypes.Array(tmp._type_, tmp, lock=False)
