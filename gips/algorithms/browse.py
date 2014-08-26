#!/usr/bin/env python

import argparse

from gips.algorithms.core import Algorithm


class Browse(Algorithm):
    """ Create browse imagery of a GIPS project """
    name = 'GIPS Browse'
    __version__ = '1.0.0'

    def run(self, **kwargs):
        for date in self.inv.dates:
            for p in self.inv.requested_products:
                

    @classmethod
    def parser(cls, parser0):
        cls.add_project_parser(parser0)
        return parser0


def main():
    Browse.main()
