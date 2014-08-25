#!/usr/bin/env python

import os
import glob
import gips.settings as settings


def rename_products(repo, pattern_in, pattern_out):
    path = os.path.join(settings.REPOS[repo]['rootpath'], 'tiles')
    for tile in os.listdir(path):
        d = os.path.join(path, tile)
        if os.path.isdir(d):
            for date in os.listdir(os.path.join(path, tile)):
                filenames = glob.glob(os.path.join(path, tile, date, '*' + pattern_in + '*'))
                for f in filenames:
                    newf = f.replace(pattern_in, pattern_out)
                    try:
                        print '%s -> %s' % (os.path.basename(f), os.path.basename(newf))
                        os.rename(f, newf)
                    except:
                        print 'problem renaming %s' % f


# GIPS v 0.7.0
# Renaming files to replace _ with - in product designations
for p in ['ref', 'rad', 'acca', 'fmask', 'bi', 'ndvi', 'evi', 'lswi', 'ndsi', 'satvi', 'satvi2', 'msavi2']:
    rename_products('landsat', p+'_', p+'-')

print 'Remaining files that were not renamed:'
print os.system('find . -name *_*_*_*_*')
