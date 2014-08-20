#!/usr/bin/env python

import os
import argparse
import random
import numpy

import gippy
from gips.core import Algorithm
from gips.utils import VerboseOut


class Truth(Algorithm):
    name = 'Truth Utilities'
    __version__ = '1.0.0'

    def merge(self, files, output, **kwargs):
        """ Merge multiple truth files into one showing agreement over all """
        VerboseOut('Merging truth files', 2)
        for i, filename in enumerate(files):
            img = gippy.GeoImage(filename)
            if i == 0:
                imgarr = img[0].Read()
                imgout = gippy.GeoImage(output, img)
                imgout.CopyColorTable(img)
                imgout.SetNoData(0)
            else:
                imgarr2 = img[0].Read()
                it = numpy.nditer(imgarr, flags=['multi_index'])
                while not it.finished:
                    (x, y) = it.multi_index
                    if imgarr[x, y] != imgarr2[x, y]:
                        imgarr[x, y] = 0
                    it.iternext()
                imgarr2 = None
            img = None
        imgout[0].Write(imgarr)
        VerboseOut('Created new truth map %s' % os.path.basename(output), 2)

    def prune(self, truthfile='', samples=10, **kwargs):
        """ Prune down a truth image to a max number of samples """
        VerboseOut('Pruning truth file %s' % os.path.basename(truthfile), 2)
        random.seed()
        suffix = '_pruned%sk' % str(samples)
        samples = samples * 1000
        gimg_truth = gippy.GeoImage(truthfile)
        timg = gimg_truth[0].Read()
        for i in range(1, numpy.max(timg)+1):
            locs = numpy.where(timg == i)
            num = len(locs[0])
            if num > samples:
                sample = random.sample(range(0, num-1), num-samples)
                locs = (locs[0][sample], locs[1][sample])
                timg[locs] = 0
            elif num < 10:
                timg[locs] = 0
            if num != 0:
                locs = numpy.where(timg == i)
                VerboseOut('Class %s: %s -> %s samples' % (i, num, len(locs[0])), 2)
        gimg_pruned = gippy.GeoImage(os.path.splitext(truthfile)[0] + suffix, gimg_truth)
        gimg_pruned.SetNoData(0)
        gimg_pruned.CopyColorTable(gimg_truth)
        gimg_pruned[0].Write(timg)

    @classmethod
    def parser(cls):
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(add_help=False)
        subparser = parser0.add_subparsers(dest='command')

        h = 'Merge multiple truth maps into single map of agreement'
        parser = subparser.add_parser('merge', parents=[cls.vparser()], help=h, formatter_class=dhf)
        parser.add_argument('files', help='List of truth files', nargs='+')
        parser.add_argument('-o', '--output', help='Output file', required=True)

        h = 'Prune Truth to S samples per class'
        parser = subparser.add_parser('prune', parents=[cls.vparser()], help=h, formatter_class=dhf)
        parser.add_argument('truthfile', help='Truth map to prune')
        parser.add_argument('-s', '--samples', help='Number of samples (in 1K units) to extract', default=10, type=int)

        return parser0


def main():
    Truth.main()
