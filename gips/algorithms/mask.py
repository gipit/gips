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

import gippy
from gips.algorithms.core import Algorithm
from gips.utils import VerboseOut, basename


class Mask(Algorithm):
    name = 'Mask'
    __version__ = '1.0.0'
    suffix = '-masked'

    def run(self, fmask='', pmask='', original=False, overwrite=False, **kwargs):
        if fmask == '' and pmask == '':
            raise Exception('No masks supplied!')
        if fmask != '':
            mask_file = gippy.GeoImage(fmask)
        for date in self.inv.dates:
            VerboseOut('Masking files from %s' % date)
            available_masks = self.inv[date].masks(pmask)
            for p in self.inv.product_list(date):
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
                    if original:
                        VerboseOut('  %s' % (img.Basename()), 2)
                        img.Process()
                        img.SetMeta('MASKS', meta)
                    else:
                        fout = os.path.splitext(img.Filename())[0] + self.suffix + '.tif'
                        if not os.path.exists(fout) or overwrite:
                            VerboseOut('  %s -> %s' % (img.Basename(), basename(fout)), 2)
                            imgout = img.Process(fout)
                            imgout.SetMeta('MASKS', meta)
                            imgout = None
                img = None
        mask_file = None

    @classmethod
    def parser(cls, parser0):
        cls.add_project_parser(parser0)
        parser0.add_argument('--fmask', help='Mask files with this file (of matching dimensions)', default='')
        parser0.add_argument('--pmask', help='Mask files with this product', nargs='*', default=[])
        h = 'Write mask to original image instead of creating new image'
        parser0.add_argument('--original', help=h, default=False, action='store_true')
        h = 'Overwrite existing files when creating new'
        parser0.add_argument('--overwrite', help=h, default=False, action='store_true')
        #parser.add_argument('-i', '--invert', help='Invert mask (0->1, 1->0)', default=False, action='store_true')
        #parser.add_argument('--value', help='Mask == val', default=1)
        return parser0


def main():
    Mask.main()
