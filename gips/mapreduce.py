#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
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

import multiprocessing
import numpy
from pdb import set_trace 

# TODO - check dimensions and call appropriate function

def init(_data, _func):
    global data, func
    data = _data
    func = _func


def worker_2d_to_1d(locs):
    """ Apply a function to a piece of the 2-D datain array and return 1-D output """
    dat = data[locs[0]:locs[1], :]
    valid = numpy.all(~numpy.isnan(dat), axis=1)
    chunk = numpy.zeros(dat.shape[0])
    chunk[valid] = func(dat[valid, :])
    return chunk


def map_reduce(data, func, nchunks=100, nproc=2):
    # chunk up input data
    chunksz = int(data.shape[0] / nchunks)
    extra = data.shape[0] - chunksz * nchunks
    chunks = [chunksz] * (nchunks - extra) + [chunksz + 1] * extra
    inputs = []
    for ichunk in range(nchunks):
        x0 = sum(chunks[:ichunk])
        x1 = x0 + chunks[ichunk]
        inputs.append((x0, x1))
    pool = multiprocessing.Pool(nproc, initializer=init, initargs=(data, func))
    data_parts = pool.map(worker_2d_to_1d, inputs)
    # reassemble data
    dataout = numpy.zeros(data.shape[0])
    for i, chunk in enumerate(inputs):
        dataout[chunk[0]:chunk[1]] = data_parts[i]
    return dataout


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
