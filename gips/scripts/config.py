#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  matt.a.hanson@gmail.com
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
import gips
from gips import __version__ as version
from gips.parsers import GIPSParser
from gips.utils import VerboseOut, create_environment_settings, create_user_settings
import pprint
import traceback


def main():
    import gips
    title = 'GIPS Configuration Utility (v%s)' % (version)

    try:
        parser = GIPSParser(description=title, datasources=False)
        subparser = parser.add_subparsers(dest='command')
        subparser.add_parser('print', help='Print current settings')
        p = subparser.add_parser('env', help='Configure GIPS repositories in this environment')
        p.add_argument('-r', '--repos', help='Top level directory for repositories', default='/data/repos')
        p.add_argument('-e', '--email', help='Set email address (used for anonymous FTP sources)', required=True)
        p = subparser.add_parser('user', help='Configure GIPS repositories for this user (for per user customizations)')
        #p.add_argument('-e', '--email', help='Set email address (used for anonymous FTP sources)')
        h = 'Install full configuration file without inheriting from environment settings'
        p.add_argument('-f', '--full', help=h, default=False, action='store_true')
        args = parser.parse_args()
        print title

        if args.command == 'print':
            try:
                from gips.utils import settings
                s = settings()
                for v in dir(s):
                    if not v.startswith('__') and v != 'gips':
                        print
                        print v
                        exec('pprint.pprint(s.%s)' % v)
            except Exception, e:
                #print traceback.format_exc()
                print 'Unable to access settings. Run `gips_config env`'

        elif args.command == 'env':
            try:
                create_environment_settings(args.repos, email=args.email)
            except Exception, e:
                print 'Could not create environment settings: %s' % e

        elif args.command == 'user':
            try:
                # first try importing environment settings
                import gips.settings
                cfgfile = create_user_settings()
                print 'User settings file: %s' % cfgfile
            except Exception, e:
                # could not import gips.settings, TODO - see if user wants to proceed with user-only config
                #raise Exception('No environment settings found...run `gips_config env` to install in environment first')
                print 'Could not create user settings: %s' % e
                
    except Exception, e:

        VerboseOut(traceback.format_exc(), 1)
        print 'GIPS error: %s' % e


if __name__ == "__main__":
    main()

