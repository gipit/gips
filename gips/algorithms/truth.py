#!/usr/bin/env python

import os
import argparse
import random
import numpy

import gippy
from gips.algorithms.core import Algorithm
from gips.utils import VerboseOut, basename
from sklearn.metrics import confusion_matrix as cmatrix
from pdb import set_trace

def confusion_matrix(truth, pred, nodata=None):
    """ Caculated confusion matrix and common metrics
        producer accuracy - % of class correctly identified (inverse is false negative rate)
        consumer accuracy - % liklihood pixel of that class is correct (reliability) inverse is false positive rate
    """
    if nodata is not None:
        valid = numpy.where(truth != nodata)
        truth = truth[valid]
        pred = pred[valid]
        valid = numpy.where(pred != nodata)
        truth = truth[valid]
        pred = pred[valid]
    cm = cmatrix(truth, pred).astype('float')
    cm_diag = numpy.diag(cm).astype('float')
    rowsums = cm.sum(axis=0)
    colsums = cm.sum(axis=1)
    total = rowsums.sum().astype('float')
    accuracy = cm_diag.sum() / total
    prod_accuracy = cm_diag / colsums
    user_accuracy = cm_diag / rowsums

    #err = (rowsums/total * colsums/total).sum()
    #kappa = (accuracy - err)/(1. - err)
    err = (rowsums/total * colsums/total).sum()
    kappa = (cm_diag.sum()/total - err)/(1. - err)
    classes = numpy.unique(numpy.concatenate((truth, pred)))
    df = pandas.DataFrame(cm, index=classes, columns=classes)
    return {'cm': df, 'accuracy': accuracy,  'kappa': kappa,
            'prod_accuracy': prod_accuracy, 'user_accuracy': user_accuracy}

    """
    # Simplifying confusion matrix
    # extracting just rice?
    loc = numpy.where(cmatrix['cm'].columns == 3)
    cm = cmatrix['cm'].values
    ricecm = numpy.zeros((2, 2))
    ricecm[1, 1] = cm[loc, loc]
    ricecm[0, 1] = (cm[loc, :].sum() - ricecm[1, 1])
    ricecm[1, 0] = (cm[:, loc].sum() - ricecm[1, 1])
    ricecm[0, 0] = cm.sum().sum() - ricecm.sum()
    rowsums = cm.sum(axis=0)
    colsums = cm.sum(axis=1)
    total = rowsums.sum().astype('float')
    err = (rowsums/total * colsums/total).sum()
    kappa = (numpy.diag(ricecm).sum()/total - err)/(1. - err)
    """

"""
def calculate_stats(cmi):  # ytrue, ypred):
    #cmi = confusion_matrix(ytrue, ypred)
    cmi = cmi.values
    npts = cmi.sum(axis=1)
    cm = cmi.copy().astype('float32')/cmi.sum()
    overall = numpy.sum(numpy.diag(cm))
    producers = numpy.diag(cm) / numpy.sum(cm, axis=1)
    consumers = numpy.diag(cm) / numpy.sum(cm, axis=0)
    ave_producers = producers.mean()
    ave_consumers = producers.mean()
    err = numpy.sum(cm.sum(axis=0) * cm.sum(axis=1))
    kappa = (overall - err)/(1. - err)
    nrm = cmi.astype('float32')/numpy.diag(cmi)
    maxmiss = numpy.triu(nrm, 1).max()
    return {'kappa': kappa, 'overall': overall, 'maxmiss': maxmiss,
            'producers': producers, 'consumers': consumers, 'npts': npts}  # , 'cm': npts }
"""


class Truth(Algorithm):
    name = 'Truth Utilities'
    __version__ = '1.0.0'

    @classmethod
    def stats(cls, files, counts=False, **kwargs):
        """ Print out counts of each class for entire array """
        fheader = '{:^8}'.format('')
        header = '{0:^8}'.format('Class')
        for f in files:
            fieldwidth = 10
            if counts:
                fieldwidth = 20
                header = header + '{0:^10}'.format('#Pixels')
            fheader = fheader + ('{0:^%s} ' % fieldwidth).format(basename(f))
            header = header + '{0:^10} '.format('Scene%')

        data = []
        numclasses = 255
        for filename in files:
            VerboseOut('Calculating class stats for %s' % basename(filename), 2)
            gimg = gippy.GeoImage(filename)
            img = gimg[0].Read()
            totalsz = gimg.Size()
            nodatasz = len(numpy.where(img == gimg[0].NoDataValue())[0])
            datasz = totalsz - nodatasz
            dat = [totalsz, nodatasz, datasz]
            #numclass = int(numpy.max(img))
            """
            if args.vector:
                #dbftable = zonal.get_dbftable(args.vector)
                #attrmap = dbftable.create_map('FID',args.field)
                rvector = zonal.rasterize_vector(args.vector, f, field=args.field)
                fh_vec = gdal.Open(rvector)
                imgv = fh_vec.GetRasterBand(1).ReadAsArray()
                numclass_vec = numpy.max(imgv)
                totalpix_vec = []
                for iclass in range(1,numclass_vec+1):
                    head = head + ' {0:^6}'.format('TS%s'%iclass)
                    totalpix_vec.append( len(numpy.where(imgv == iclass)[0]) )
            """
            for i in range(1, numclasses+1):
                inds = numpy.where(img == i)
                num = float(len(inds[0]))
                dat.append(num)
                """
                if args.vector:
                    img_vec_sub = imgv[inds]
                    for iTS in range(1, numclass_vec+1):
                        vinds = numpy.where(img_vec_sub == iTS)
                        numpix = float(len(vinds[0]))
                        line = line + ' {ts:^6.1f}'.format(ts=numpix/totalpix_vec[iTS-1]*100)
                """
            data.append(dat)
            gimg = None

        print fheader
        print header
        #nodatastr = ['%4.2f' % float(numnodata)/sz*100]
        #print 'NoData: %4.2f%% of scene' % (float(numnodata)/sz*100)
        line = '{0:^8}'.format('NoData')
        for f in range(0, len(files)):
            if counts:
                line = line + '{0:<10}'.format(data[f][1])
            line = line + '{0:^10.2f} '.format(data[f][1]/data[f][0])
        print line
        for i in range(1, numclasses+1):
            totalsum = 0
            line = '{0:^8}'.format(i)
            for f in range(0, len(files)):
                num = data[f][2+i]
                totalsum = totalsum + num
                if counts:
                    line = line + '{0:<10}'.format(int(num))
                line = line + '{0:^10.2f} '.format(100*num/float(data[f][2]))
            if totalsum > 0:
                print line

    @classmethod
    def merge(cls, files, output, group=None, **kwargs):
        """ Merge multiple truth files into one showing agreement over all """
        def _group(cls, arr, groups):
            for row in groups:
                for n in row[1:]:
                    arr[arr == n] = row[0]
            return arr
        VerboseOut('Merging truth files', 2)
        for i, filename in enumerate(files):
            img = gippy.GeoImage(filename)
            if i == 0:
                imgarr = img[0].Read()
                if group is not None:
                    imgarr = cls._group(imgarr, group)
                imgout = gippy.GeoImage(output, img)
                imgout.CopyColorTable(img)
                imgout.SetNoData(0)
            else:
                imgarr2 = img[0].Read()
                if group is not None:
                    imgarr2 = cls._group(imgarr2, group)
                it = numpy.nditer(imgarr, flags=['multi_index'])
                while not it.finished:
                    (x, y) = it.multi_index
                    if imgarr[x, y] != imgarr2[x, y]:
                        imgarr[x, y] = 0
                    it.iternext()
                imgarr2 = None
            img = None
        imgout[0].Write(imgarr)
        VerboseOut('Created new truth map %s' % os.path.basename(output), 2)

    @classmethod
    def prune(cls, truthfile='', samples=10, **kwargs):
        """ Prune down a truth image to a max number of samples """
        VerboseOut('Pruning truth file %s' % os.path.basename(truthfile), 2)
        random.seed()
        suffix = '_pruned%sk' % str(samples)
        samples = samples * 1000
        gimg_truth = gippy.GeoImage(truthfile)
        timg = gimg_truth[0].Read()
        for i in range(1, numpy.max(timg)+1):
            locs = numpy.where(timg == i)
            num = len(locs[0])
            if num > samples:
                sample = random.sample(range(0, num-1), num-samples)
                locs = (locs[0][sample], locs[1][sample])
                timg[locs] = 0
            elif num < 10:
                timg[locs] = 0
            if num != 0:
                locs = numpy.where(timg == i)
                VerboseOut('Class %s: %s -> %s samples' % (i, num, len(locs[0])), 2)
        gimg_pruned = gippy.GeoImage(os.path.splitext(truthfile)[0] + suffix, gimg_truth)
        gimg_pruned.SetNoData(0)
        gimg_pruned.CopyColorTable(gimg_truth)
        gimg_pruned[0].Write(timg)

    def analyze(cls, datadir, truth, **kwargs):
        """ Compare output map(s) with a truth map """

        gimg_truth = gippy.GeoImage(truth)
        truthimg = gimg_truth[0].Read()

        header = 'Day, Overall Accuracy, Rice Kappa, Rice Accuracy, Rice Reliability, Rice (ha)'
        fout = open(os.path.join(datadir, 'summary.csv'), 'w')
        fout.write(header+'\n')

        for filename in glob.glob(os.path.join(datadir, '*master.tif')):
            gimg = gippy.GeoImage(filename)
            predimg = gimg[0].Read()

            # TODO - better extraction of day
            day = gimg.Basename()[0:3]

            cmatrix = confusion_matrix(truthimg, predimg, nodata=gimg[0].NoDataValue())
            cmatrix.to_csv(os.path.join(datadir, '%s.cmatrix.csv' % day))

            fout.write('%s, %4.2f, %4.2f, %4.2f\n' % (day, cmatrix['accuracy'], cmatrix['kappa']))

        fout.close()

        """
        if args.truth is None:
            header = 'Day, Rice (ha)'
                        # get total rice, in hectares
            ricetotal = len(numpy.where(pred == 3)[0]) * 0.09

            if args.truth is None:
                line = '%s, %4.2f' % (day, ricetotal)
        """

    @classmethod
    def parser(cls, parser0):
        # Add a subparser and get keywords to add to commands
        subparser, kwargs = cls.subparser(parser0)

        h = 'Merge multiple truth maps and/or multiple values into single map of agreement'
        parser = subparser.add_parser('merge', help=h, **kwargs)
        parser.add_argument('files', help='List of truth files', nargs='+')
        parser.add_argument('-o', '--output', help='Output file', required=True)
        parser.add_argument('-g', '--group', nargs='*', action='append', default=None, type=int)

        h = 'Prune truth to S samples per class'
        parser = subparser.add_parser('prune', help=h, **kwargs)
        parser.add_argument('truthfile', help='Truth map to prune')
        parser.add_argument('-s', '--samples', help='Number of samples (in 1K units) to extract', default=10, type=int)

        h = 'Compare truth maps output map(s)'
        parser = subparser.add_parser('analyze', help=h, **kwargs)
        parser.add_argument('-d', '--datadir', help='Input data directory of classmaps', required=True)
        parser.add_argument('-t', '--truth', help='Truth image (not pruned)', required=True)

        h = 'Print stats for each class in classification map'
        parser = subparser.add_parser('stats', help=h, **kwargs)
        parser.add_argument('files', help='List of truth files', nargs='+')
        parser.add_argument('--counts', help='Include total pixel count', default=False, action='store_true')
        #parser0.add_argument('-v','--vector',help='Count pixels in each class for provided shapefile training data')
        #parser0.add_argument('-f','--field', help='The field in vector file to use as training set', default='FID')

        return parser0


def main():
    Truth.main()
