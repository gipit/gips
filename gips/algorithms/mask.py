#!/usr/bin/python

import os
import argparse
from datetime import datetime
import numpy as np

import gippy
from gips.algorithms import core
from gips.utils import VerboseOut


class Mask(object):
    name = 'Mask'
    __version__ = '0.1.0'
    suffix = '_masked'

    def __init__(self, inv, file_mask='', product_mask='', overwrite=False, *args, **kwargs):
        if file_mask == '' and product_mask == '':
            raise Exception('No masks supplied!')
        if file_mask != '':
            mask_file = gippy.GeoImage(file_mask)
        for date in sorted(inv):
            VerboseOut('%s' % date)
            if product_mask != '':
                mask_product = gippy.GeoImage(inv[date][product_mask])
            for p in inv[date]:
                if product_mask != p:
                    fname = inv[date][p]
                    img = gippy.GeoImage(fname)
                    maskit = False
                    if file_mask != '':
                        img.AddMask(mask_file[0])
                        maskit = True
                    if product_mask != '':
                        img.AddMask(mask_product[0])
                        maskit = True
                    if maskit:
                        VerboseOut('Masking %s' % fname, 2)
                        if overwrite:
                            img.Process()
                        else:
                            fout = os.path.splitext(fname)[0]
                            img.Process(fout+self.suffix)
                    img.ClearMasks()
                    img = None
            mask_product = None
        mask_file = None

    @classmethod
    def arg_parser(cls):
        parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-f', '--file_mask', help='Mask files with this file (of matching dimensions)', default='')
        parser.add_argument('-p', '--product_mask', help='Mask files with this product (in project directory)', default='')
        parser.add_argument('--overwrite', help='Overwrite existing files', default=False, action='store_true')
        #parser.add_argument('-i', '--invert', help='Invert mask (0->1, 1->0)', default=False, action='store_true')
        #parser.add_argument('--value', help='Mask == val', default=1)
        return parser


def main():
    core.main(Mask)


if __name__ == "__main__":
    main()
