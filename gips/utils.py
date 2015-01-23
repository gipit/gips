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
import multiprocessing
import numpy
import tempfile
import commands
import shutil
from gips import GeoVector


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


def crop2vector(img, vector):
    """ Crop a GeoImage down to a vector - only used by mosaic """
    # TODO - incorporate into GIPPY?
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
    VerboseOut('Cropped to vector in %s' % (datetime.now() - start), 3)
    return img


def mosaic(infiles, outfile, vectorfile):
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
    return imgout
    # TODO - fix crop2vector, throwing out of bounds error in GDAL read
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
