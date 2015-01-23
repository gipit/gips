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

import numpy
import multiprocessing


def _worker(chunk):
    """ Worker function (has access to global variables set in _mr_init """
    # read chunk of data and make sure it is 3-D: BxYxX
    data = rfunc(chunk)
    shape = data.shape
    if len(shape) == 2:
        data = data.reshape((1, shape[0], shape[1]))
        shape = data.shape

    # make output array for this chunk
    output = numpy.empty((outshape[0], shape[1], shape[2]))
    output[:] = numpy.nan

    # only run on valid pixel signatures unless keepnodata set
    if keepnodata:
        valid = numpy.ones((shape[1], shape[2])).astype('bool')
    else:
        valid = numpy.all(~numpy.isnan(data), axis=0)

    # run processing function
    output[:, valid] = pfunc(data[:, valid])

    # write using write function if provided
    if wfunc is not None:
        wfunc((output, chunk))
        return None
    else:
        return output


class MapReduce(object):
    """ General purpose class for performing map reduction functions """

    def __init__(self, inshape, outshape, rfunc, pfunc, wfunc=None, nproc=2, keepnodata=False):
        """ Create multiprocessing pool """
        self.inshape = inshape
        self.outshape = outshape
        self.pool = multiprocessing.Pool(nproc, initializer=self._mr_init,
                                         initargs=(inshape, outshape, rfunc, pfunc, wfunc, keepnodata))

    def run(self, nchunks=100, chunks=None):
        """ Run the multiprocessing pool """
        if chunks is None:
            self.chunks = self.chunk(self.inshape, nchunks=nchunks)
        else:
            self.chunks = chunks
        self.dataparts = self.pool.map(_worker, self.chunks)

    def assemble(self):
        """ Reassemble output parts into single array """
        dataout = numpy.empty(self.outshape)
        for i, ch in enumerate(self.chunks):
            dataout[:, ch[1]:ch[1] + ch[3], ch[0]:ch[0] + ch[2]] = self.dataparts[i]
        return dataout.squeeze()

    @staticmethod
    def _mr_init(_inshape, _outshape, _rfunc, _pfunc, _wfunc, _keepnodata):
        """ Initializer sets globals for processes """
        global inshape, outshape, rfunc, pfunc, wfunc, keepnodata
        inshape = _inshape
        outshape = _outshape
        rfunc = _rfunc
        pfunc = _pfunc
        wfunc = _wfunc
        keepnodata = _keepnodata

    @staticmethod
    def chunk(shape, nchunks=100):
        """ Create chunks given input data size (B x Y x X) """
        chunksz = int(shape[1] / nchunks)
        remainder = shape[1] - chunksz * nchunks
        chszs = [chunksz] * (nchunks - remainder) + [chunksz + 1] * remainder
        chunks = []
        for ichunk in range(nchunks):
            # This is being inverted because gippy is X x Y, whereas numpy is Y x X
            #chunks.append(gippy.Recti(0, sum(chszs[:ichunk]), datasz[2], chszs[ichunk]))
            chunks.append([0, sum(chszs[:ichunk]), shape[2], chszs[ichunk]])
        return chunks

    @staticmethod
    def get_shapes(arrin, numbands):
        """ Create in and out shapes based on input array and output numbands) """
        inshape = arrin.shape
        if len(inshape) == 2:
            arrin = arrin.reshape((1, inshape))
            inshape = arrin.shape
        outshape = (numbands, inshape[1], inshape[2])
        return (inshape, outshape)


def map_reduce_array(arrin, pfunc, numbands=1, nchunks=100, nproc=2, keepnodata=False):
    """ Apply user defined pfunc to a numpy array using multiple processors """
    (inshape, outshape) = MapReduce.get_shapes(arrin, numbands)

    # read data from global input array
    rfunc = lambda chunk: arrin[:, chunk[1]:chunk[1] + chunk[3], chunk[0]:chunk[0] + chunk[2]]

    mr = MapReduce(inshape, outshape, rfunc=rfunc, pfunc=pfunc, nproc=nproc, keepnodata=keepnodata)
    mr.run(nchunks=nchunks)
    return mr.assemble()


def _test_map_reduce_array(arrin, pfunc, numbands=1, nchunks=100, nproc=2, keepnodata=False):
    """ Test map_reduce_array functions without using multiprocessing """

    (inshape, outshape) = MapReduce.get_shapes(arrin, numbands)

    # read data from global input array
    rfunc = lambda chunk: arrin[:, chunk[1]:chunk[1] + chunk[3], chunk[0]:chunk[0] + chunk[2]]

    chunks = MapReduce.chunk(inshape, nchunks=nchunks)

    MapReduce._mr_init(inshape, outshape, rfunc, pfunc, None, keepnodata)

    dataout = numpy.empty(outshape)
    for ch in chunks:
        dataout[:, ch[1]:ch[1] + ch[3], ch[0]:ch[0] + ch[2]] = _worker(ch)
    return dataout
