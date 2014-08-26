#!/usr/bin/env python

import argparse

from gips.algorithms.core import Algorithm


class Project(Algorithm):
    name = 'GIPS Project Utilities'
    __version__ = '1.0.0'

    def browse(self, **kwargs):
        """ Create browse imagery for project """
        for date in self.inv.dates:
            for p in self.inv.requested_products:
                fname = self.inv[date][p]
                print fname

    @classmethod
    def parser(cls, parser0):
        subparser = cls.subparser(parser0)
        kwargs = cls.parser_kwargs(project=True)

        parser = subparser.add_parser('browse', help='Create browse imagery', **kwargs)

        return parser0


def main():
    Project.main()
