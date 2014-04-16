#!/usr/bin/env python
################################################################################
#    GIPPY: Geospatial Image Processing library for Python
#
#    Copyright (C) 2014 Matthew A Hanson
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

import os, argparse
import datetime
import gippy

def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='Geospatial Image Process It', formatter_class=dhf)
    subparser = parser0.add_subparsers(dest='command', help='Operations')

    # Global options
    gparser = argparse.ArgumentParser(add_help=False, formatter_class=dhf)
    gparser.add_argument('files', nargs='*', help='Image files to process')
    gparser.add_argument('-v','--verbose', help='Verbosity level (0-4)', default=1, type=int)

    # Output file options
    oparser = argparse.ArgumentParser(add_help=False)
    group = oparser.add_argument_group('Output file options')
    group.add_argument('--format', help='Output format', default='GTiff')
    #group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_algname')

    # ALGORITHMS
    #parser = subparser.add_parser('CreateMask',help='Create valid pixel mask based on all bands', parents=[gparser,oparser], formatter_class=dhf)
    #group = parser.add_argument_group('algorithm arguments')
    
    #parser = subparser.add_parser('FixBadPixels',help='Replace Inf/NaN results with NoData value', parents=[gparser], formatter_class=hformat)
   
    #parser = subparser.add_parser('ApplyMask', help='Apply Mask to existing image', parents=gparser, formatter_class=hformat)
    #parser.add_argument('-m','--mask', help='Mask to apply', required=True)

    # TODO - automate reading in all algorithms in gippy.algorithms
    algs = ['kmeans','thresh','acca','acloud','fmask']
    for a in algs:
        exec('import gippy.algorithms.%s as %s' % (a, a))
        eval('%s.add_options(subparser,[gparser,oparser])' % a)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetDefaultFormat(args.format)

    for f in args.files:
        start = datetime.datetime.now()
        fbase,ext = os.path.splitext(os.path.basename(f))
        #if args.suffix == '_algname': args.suffix = '_' + args.command
        #fout = fbase + '_' + args.suffix
        imgout = eval('%s.process(f, **args.__dict__)' % args.command)
        #if args.command == 'CreateMask':
        #    fout = fbase + args.suffix # + ext
        #    exec("gippy.CreateMask(gippy.GeoImage('%s'),'%s')" % (f,fout))
        #    
        #elif args.command == 'FixBadPixels':
        #    exec("gippy.FixBadPixels(gippy.GeoImage('%s'))" % (f))

        print "%s -> %s: %s" % (fbase, imgout.Basename(), datetime.datetime.now()-start)
            
