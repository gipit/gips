#!/usr/bin/env python

import os, argparse
import datetime
import gippy

if __name__ == "__main__":
    prog = os.path.split(__file__)[1]
    hformat = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(prog=prog, formatter_class=hformat, description='Geospatial Image Process IT!')
    subparser = parser0.add_subparsers(help='Operations', dest='command')

    # Global options
    gparser = argparse.ArgumentParser(add_help=False, formatter_class=hformat)
    gparser.add_argument('files', nargs='*', help='Image files to process')
    gparser.add_argument('-v','--verbose', help='Verbosity level (0-3)', default=1, type=int)

    # Output file options
    oparser = argparse.ArgumentParser(add_help=False)
    group = oparser.add_argument_group('Output file options')
    group.add_argument('--format', help='Output format', default='GTiff')

    # ALGORITHMS
    parser = subparser.add_parser('CreateMask',help='Create valid pixel mask based on all bands', parents=[gparser,oparser], formatter_class=hformat)
    group = parser.add_argument_group('algorithm arguments')
    group.add_argument('-s','--suffix',help='Append suffix to filename for output', default='_mask')

    parser = subparser.add_parser('FixBadPixels',help='Replace Inf/NaN results with NoData value', parents=[gparser], formatter_class=hformat)
   
    #parser = subparser.add_parser('ApplyMask', help='Apply Mask to existing image', parents=gparser, formatter_class=hformat)
    #parser.add_argument('-m','--mask', help='Mask to apply', required=True)

    # TODO - automate reading in all algorithms in gippy.algorithms
    algs = ['kmeans','thresh','index2prob','permutations']
    for a in algs:
        exec('import gippy.algorithms.%s as %s' % (a, a))
        eval('%s.add_options(subparser,[gparser,oparser])' % a)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetDefaultFormat(args.format)

    for f in args.files:
        start = datetime.datetime.now()
        fbase,ext = os.path.splitext(os.path.basename(f))
        if args.command == 'CreateMask':
            fout = fbase + args.suffix # + ext
            exec("gippy.CreateMask(gippy.GeoImage('%s'),'%s')" % (f,fout))
            print "%s -> %s: %s" % (fbase, os.path.basename(fout), datetime.datetime.now()-start)
        elif args.command == 'FixBadPixels':
            exec("gippy.FixBadPixels(gippy.GeoImage('%s'))" % (f))
            print "%s: %s" % (fbase, datetime.datetime.now()-start)
        else:
            eval('%s.process(f, **args.__dict__)' % args.command)