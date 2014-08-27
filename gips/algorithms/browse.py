#!/usr/bin/env python

import argparse

import gippy
from gips.algorithms.core import Algorithm


class Browse(Algorithm):
    """ Create browse imagery of a GIPS project """
    name = 'GIPS Browse'
    __version__ = '1.0.0'

    def run(self, quality=75, **kwargs):
        for date in self.inv.dates:
            for p in self.inv.product_list(date):
                try:
                    filename = self.inv[date]
                    gippy.BrowseImage(self.inv[date].open(p), quality)
                except:
                    pass

    @classmethod
    def parser(cls, parser0):
        cls.add_project_parser(parser0)
        parser0.add_argument('-q', '--quality', help='JPG Quality', default=75)
        return parser0


def main():
    Browse.main()
