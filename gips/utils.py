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


def atmospheric_model(doy, lat):
    """ Determine atmospheric model
    1 - Tropical
    2 - Mid-Latitude Summer
    3 - Mid-Latitude Winter
    4 - Sub-Arctic Summer
    5 - Sub-Arctic Winter
    6 - US Standard Atmosphere
    """
    # Determine season
    if doy < 121 or doy > 274:
        if lat < 0:
            summer = True
        else:
            summer = False
    else:
        if lat < 0:
            summer = False
        else:
            summer = True
    # Determine model
    if abs(lat) <= 15:
        model = 1
    elif abs(lat) >= 60:
        if summer:
            model = 4
        else:
            model = 5
    else:
        if summer:
            model = 2
        else:
            model = 3
    return model


def chunk_data(datasz, nchunks=100):
    """ Create chunks given input data size """
    if len(datasz) == 3:
        chaxis = 1
    else:
        chaxis = 0
    chunksz = int(datasz[chaxis] / nchunks)
    remainder = datasz[chaxis] - chunksz * nchunks
    chszs = [chunksz] * (nchunks - remainder) + [chunksz + 1] * remainder
    chunks = []
    for ichunk in range(nchunks):
        # This is being inverted because gippy is X x Y, whereas numpy is Y x X
        #chunks.append(gippy.Recti(0, sum(chszs[:ichunk]), datasz[2], chszs[ichunk]))
        chunks.append([0, sum(chszs[:ichunk]), datasz[2], chszs[ichunk]])
    return chunks


def _mr_init(_readfunc, _func):
    """ Put functions witin global namespace for each process """
    global readfunc, func
    readfunc = _readfunc
    func = _func


def _mr_worker_3d_to_2d(chunk):
    """ Reduces multiple band image (bands x rows x cols) to single band image (rows x cols) """
    data = readfunc(gippy.Recti(chunk[0], chunk[1], chunk[2], chunk[3]))
    valid = numpy.all(~numpy.isnan(data), axis=0)
    output = numpy.zeros((data.shape[1], data.shape[2]))
    output[valid] = func(data[:, valid])
    data = None
    return output


def map_reduce(datasz, readfunc, func, nchunks=100, nproc=2):
    """ Chunk up data read from readfunc, apply func, then reassemble into output array """
    chunks = chunk_data(datasz, nchunks=nchunks)
    pool = multiprocessing.Pool(nproc, initializer=_mr_init, initargs=(readfunc, func))
    dataparts = pool.map(_mr_worker_3d_to_2d, chunks)
    # reassemble data
    dataout = numpy.zeros((datasz[1], datasz[2]))
    for i, ch in enumerate(chunks):
        dataout[ch[1]:ch[1] + ch[3], ch[0]:ch[0] + ch[2]] = dataparts[i]
    return dataout


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
    return crop2vector(imgout, vector)


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
