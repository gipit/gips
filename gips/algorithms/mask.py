#!/usr/bin/python

import os
import argparse
from datetime import datetime
import numpy as np

import gippy
from gips.core import Algorithm
from gips.utils import VerboseOut
from gips.inventory import project_inventory


class Mask(Algorithm):
    name = 'Mask'
    __version__ = '0.1.0'
    suffix = '_masked'

    def __init__(self, project, fmask='', pmask='', overwrite=False, **kwargs):
        inv = project_inventory(project)
        if fmask == '' and pmask == '':
            raise Exception('No masks supplied!')
        if fmask != '':
            mask_file = gippy.GeoImage(fmask)
        for date in sorted(inv):
            VerboseOut('%s' % date)
            if pmask != '':
                mask_product = gippy.GeoImage(inv[date][pmask])
            for p in inv[date]:
                if pmask != p:
                    fname = inv[date][p]
                    img = gippy.GeoImage(fname)
                    maskit = False
                    if fmask != '':
                        img.AddMask(mask_file[0])
                        maskit = True
                    if pmask != '':
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
    def parser(cls):
        parser = argparse.ArgumentParser(add_help=False, parents=[cls.project_parser()])
        parser.add_argument('--fmask', help='Mask files with this file (of matching dimensions)', default='')
        parser.add_argument('--pmask', help='Mask files with this product (in project directory)', default='')
        parser.add_argument('--overwrite', help='Overwrite existing files', default=False, action='store_true')
        #parser.add_argument('-i', '--invert', help='Invert mask (0->1, 1->0)', default=False, action='store_true')
        #parser.add_argument('--value', help='Mask == val', default=1)
        return parser


def main():
    Mask.main()


#if __name__ == "__main__":
#    main()
