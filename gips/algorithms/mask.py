#!/usr/bin/python

import os
import argparse
from datetime import datetime
import numpy as np

import gippy
from gips.algorithms.core import Algorithm
from gips.utils import VerboseOut, basename
from gips.inventory import ProjectInventory
from pdb import set_trace


class Mask(Algorithm):
    name = 'Mask'
    __version__ = '1.0.0'
    suffix = '-masked'

    def run(self, fmask='', pmask='', overwrite=False, **kwargs):
        if fmask == '' and pmask == '':
            raise Exception('No masks supplied!')
        if fmask != '':
            mask_file = gippy.GeoImage(fmask)
        for date in self.inv.dates:
            VerboseOut('Masking files from %s' % date)
            available_masks = self.inv[date].masks(pmask)
            for p in self.inv[date].products:
                # don't mask any masks
                if p in available_masks:
                    continue
                meta = ''
                img = self.inv[date].open(p, True)
                if fmask != '':
                    img.AddMask(mask_file[0])
                    meta = basename(fmask) + ' '
                for mask in available_masks:
                    img.AddMask(self.inv[date].open(mask)[0])
                    meta = meta + basename(self.inv[date][mask]) + ' '
                if meta != '':
                    VerboseOut('Masking %s' % img.Basename(), 2)
                    if overwrite:
                        img.Process()
                        imgout.SetMeta('MASKS', meta)
                    else:
                        fout = os.path.splitext(img.Filename())[0]
                        imgout = img.Process(fout+self.suffix)
                        imgout.SetMeta('MASKS', meta)
                        imgout = None
                img = None
        mask_file = None

    @classmethod
    def parser(cls, parser0):
        cls.add_project_parser(parser0)
        parser0.add_argument('--fmask', help='Mask files with this file (of matching dimensions)', default='')
        parser0.add_argument('--pmask', help='Mask files with this product', nargs='*', default=[])
        parser0.add_argument('--overwrite', help='Overwrite existing files', default=False, action='store_true')
        #parser.add_argument('-i', '--invert', help='Invert mask (0->1, 1->0)', default=False, action='store_true')
        #parser.add_argument('--value', help='Mask == val', default=1)
        return parser0


def main():
    Mask.main()
