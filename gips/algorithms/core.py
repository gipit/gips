#!/usr/bin/env python

import argparse
from datetime import datetime
import traceback

import gippy
import gips
from gips.utils import VerboseOut
from gips.inventory import ProjectInventory


class Algorithm(object):
    name = 'Algorithm Name'
    __version__ = '0.0.0'

    def __init__(self, **kwargs):
        if 'projdir' in kwargs:
            self.inv = ProjectInventory(kwargs['projdir'], kwargs.get('products'))

    def run_command(self, **kwargs):
        """ Calls "run" function, or "command" if algorithm uses subparser """
        start = datetime.now()
        if 'command' not in kwargs:
            command = 'run'
        else:
            command = kwargs['command']
            #VerboseOut('Running %s' % command, 2)
        exec('self.%s(**kwargs)' % command)
        #VerboseOut('Completed %s in %s' % (command, datetime.now()-start), 2)
        pass

    def run(self, **kwargs):
        pass

    @classmethod
    def info(cls):
        """ Name and versions of algorithm and GIPS library """
        return 'GIPS (v%s): %s (v%s)' % (gips.__version__, cls.name, cls.__version__)

    @classmethod
    def parser(cls, parser0):
        """ Parser for algorithm specific options (defined by children) """
        return parser0

    @classmethod
    def subparser(cls, parser0, project=False):
        """ Add subparser to parser and return keywords user to include """
        subparser = parser0.add_subparsers(dest='command')
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2+: debug', default=1, type=int)
        if project:
            parser = cls.add_project_parser(parser)
        kwargs = {
            'formatter_class': argparse.ArgumentDefaultsHelpFormatter,
            'parents': [parser],
        }
        return (subparser, kwargs)

    @classmethod
    def add_project_parser(cls, parser):
        group = parser.add_argument_group('project options')
        group.add_argument('projdir', help='GIPS Project directory')
        group.add_argument('-p', '--products', help='Products to operate on', nargs='*')
        return parser

    @classmethod
    def main(cls):
        """ Main for algorithm classes """
        dhf = argparse.ArgumentDefaultsHelpFormatter

        # Top level parser
        parser = argparse.ArgumentParser(formatter_class=dhf, description=cls.info())
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2+: debug', default=1, type=int)
        parser = cls.parser(parser)

        args = parser.parse_args()
        gippy.Options.SetVerbose(args.verbose)
        VerboseOut(cls.info())

        try:
            alg = cls(**vars(args))
            alg.run_command(**vars(args))
        except Exception, e:
            VerboseOut('Error in %s: %s' % (cls.name, e))
            VerboseOut(traceback.format_exc(), 3)
