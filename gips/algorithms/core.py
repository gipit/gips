#!/usr/bin/env python

import argparse
from datetime import datetime

import gippy
import gips
from gips.utils import VerboseOut


class Algorithm(object):
    name = 'Algorithm Name'
    __version__ = '0.0.0'

    def __init__(self, **kwargs):
        if 'projdir' in kwargs:
            self.inv = ProjectInventory(kwargs['projdir'], kwargs.get('products'))

    def _run(self, **kwargs):
        """ Calls "run_all" function, or "command" if algorithm uses subparser """
        start = datetime.now()
        if 'command' not in kwargs:
            command = 'run_all'
        else:
            command = kwargs['command']
            VerboseOut('Running %s' % command, 2)
        exec('self.%s(**kwargs)' % command)
        VerboseOut('Completed %s in %s' % (command, datetime.now()-start), 2)
        pass

    def run(self, **kwargs):
        pass

    @classmethod
    def info(cls):
        """ Name and versions of algorithm and GIPS library """
        return 'GIPS (v%s): %s (v%s)' % (gips.__version__, cls.name, cls.__version__)

    @classmethod
    def parser(cls):
        """ Parser for algorithm specific options """
        parser = argparse.ArgumentParser(add_help=False)
        return parser

    @classmethod
    def project_parser(cls):
        """ Parser for using GIPS project directory """
        parser = argparse.ArgumentParser(add_help=False)
        group = parser.add_argument_group('project options')
        group.add_argument('projdir', help='GIPS Project directory')
        group.add_argument('-p', '--products', help='Products to operate on', nargs='*')
        #parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2+: debug', default=1, type=int)
        return parser

    @classmethod
    def vparser(cls):
        """ Parser for adding verbose keyword """
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2+: debug', default=1, type=int)
        return parser

    @classmethod
    def main(cls):
        """ Main for algorithm classes """
        dhf = argparse.ArgumentDefaultsHelpFormatter

        # Top level parser
        p0 = [cls.parser(), cls.vparser()]
        parser = argparse.ArgumentParser(formatter_class=dhf, parents=p0, description=cls.info())

        args = parser.parse_args()
        gippy.Options.SetVerbose(args.verbose)
        VerboseOut(cls.info())

        try:
            alg = cls(**vars(args))
            alg.run(**vars(args))
        except Exception, e:
            VerboseOut('Error in %s: %s' % (cls.name, e))
            VerboseOut(traceback.format_exc(), 3)
