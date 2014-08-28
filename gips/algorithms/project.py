#!/usr/bin/env python

import argparse

import gippy
from gips.algorithms.core import Algorithm
from pdb import set_trace


class Project(Algorithm):
    name = 'GIPS Project Utilities'
    __version__ = '1.0.0'

    def inventory(self, **kwargs):
        self.inv.pprint()

    def browse(self, quality=75, **kwargs):
        for date in self.inv.dates:
            for p in self.inv.product_list(date):
                try:
                    filename = self.inv[date]
                    gippy.BrowseImage(self.inv[date].open(p), quality)
                except:
                    pass

    @classmethod
    def parser(cls, parser0):
        subparser, kwargs = cls.subparser(parser0, project=True)

        parser = subparser.add_parser('browse', help='Create browse imagery', **kwargs)
        parser.add_argument('-q', '--quality', help='JPG Quality', default=75)

        parser = subparser.add_parser('inventory', help='Print project inventory', **kwargs)

        return parser0


def main():
    Project.main()
