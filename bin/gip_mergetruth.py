#!/usr/bin/env python

import argparse
import random

import gippy
from gips.utils import VerboseOut


def merge_truth(files, output):
    """ Merge multiple truth files into one showing agreement over all """
    random.seed()
    VerboseOut('Merging truth files', 2)
    for i, filename in files:
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
    VerboseOut('Created new truth map %s', os.path.basename(output), 2)


if __name__ == "__main__":
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description='GIP: Merge Truth files', 
        help='Merge multiple truth maps into single map where agreement', formatter_class=dhf)

    parser.add_argument('files', help='List of truth files', nargs='+')
    parser.add_argument('-o', '--output', help='Output file', default='truth.tif')

    args = parser.parse_args()

    merge_truth(files=args.files, output=args.output)