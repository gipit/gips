#!/usr/bin/env python

import argparse

from gips.algorithms.core import Algorithm


class Project(Algorithm):
    name = 'GIPS Project Utilities'
    __version__ = '1.0.0'


    @classmethod
    def parser(cls):
        parser0 = argparse.ArgumentParser(add_help=False, parents=[cls.project_parser()])

        return parser0


def main():
    Project.main()
