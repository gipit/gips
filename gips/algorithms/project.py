#!/usr/bin/env python

import argparse

import gippy
from gips.algorithms.core import Algorithm
from pdb import set_trace


class Project(Algorithm):
    name = 'GIPS Project Utilities'
    __version__ = '1.0.0'

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
        subparser = cls.subparser(parser0)
        kwargs = cls.parser_kwargs(project=True)

        parser = subparser.add_parser('browse', help='Create browse imagery', **kwargs)
        parser0.add_argument('-q', '--quality', help='JPG Quality', default=75)

        return parser0


def main():
    Project.main()
