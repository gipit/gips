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

import numpy
import multiprocessing
import gippy


class MapReduce(object):

    def __init__(self, inshape, outshape, rfunc, pfunc, wfunc=None, nchunks=100, nproc=2, keepnodata=False):
        self.inshape = inshape
        self.outshape = outshape
        self.chunks = self.chunk(inshape, nchunks=nchunks)
        pool = multiprocessing.Pool(nproc, initializer=self._mr_init,
                                    initargs=(inshape, outshape, rfunc, pfunc, wfunc, keepnodata))
        self.dataparts = pool.map(self._worker, self.chunks)

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
    def _worker(chunk):
        """ Worker function (has access to global variables set in _mr_init """
        ch = gippy.Recti(chunk[0], chunk[1], chunk[2], chunk[3])
        data = rfunc(ch)
        shape = data.shape
        output = numpy.empty((outshape, shape[1], shape[2]))
        output[:] = numpy.nan
        if keepnodata:
            valid = numpy.ones((data.shape[1], data.shape[2]))
        else:
            valid = numpy.all(~numpy.isnan(data), axis=0)
        output[:, valid] = pfunc(data[:, valid])
        data = None
        if wfunc is not None:
            wfunc((output, ch))
            return None
        else:
            return output

    @classmethod
    def chunk(cls, shape, nchunks=100):
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

    def assemble(self):
        """ Reassmble output parts into single array """
        dataout = numpy.empty(self.outshape)
        for i, ch in enumerate(self.chunks):
            dataout[:, ch[1]:ch[1] + ch[3], ch[0]:ch[0] + ch[2]] = self.dataparts[i]
        return dataout.squeeze()


def map_reduce_array(arrin, pfunc, numbands=1, nchunks=100, nproc=2):
    inshape = arrin.shape
    if len(inshape) == 2:
        inshape = (1, inshape[0], inshape[1])
    outshape = (numbands, inshape[1], inshape[2])

    def reader(chunk):
        """ Default reader - read data from global input array """
        ch = gippy.Recti(chunk[0], chunk[1], chunk[2], chunk[3])
        data = arrin[:, ch[1]:ch[1] + ch[3], ch[0]:ch[0] + ch[2]]
        return data

    mr = MapReduce(inshape, outshape, rfunc=reader, pfunc=pfunc, nchunks=nchunks, nproc=nproc)
    return mr.assemble()


#class MapReduceImage(MapReduce):
