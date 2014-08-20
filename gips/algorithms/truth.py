#!/usr/bin/env python

import os
import argparse
import random
import numpy

import gippy
from gips.core import Algorithm
from gips.utils import VerboseOut


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
    cm = confusion_matrix(truth, pred).astype('float')
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

    def merge(self, files, output, **kwargs):
        """ Merge multiple truth files into one showing agreement over all """
        VerboseOut('Merging truth files', 2)
        for i, filename in enumerate(files):
            img = gippy.GeoImage(filename)
            if i == 0:
                imgarr = img[0].Read()
                imgout = gippy.GeoImage(output, img)
                imgout.CopyColorTable(img)
                imgout.SetNoData(0)
            else:
                imgarr2 = img[0].Read()
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

    def prune(self, truthfile='', samples=10, **kwargs):
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

    def analyze(self, datadir, truth, **kwargs):
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
    def parser(cls):
        dhf = argparse.ArgumentDefaultsHelpFormatter
        parser0 = argparse.ArgumentParser(add_help=False)
        subparser = parser0.add_subparsers(dest='command')

        h = 'Merge multiple truth maps into single map of agreement'
        parser = subparser.add_parser('merge', parents=[cls.vparser()], help=h, formatter_class=dhf)
        parser.add_argument('files', help='List of truth files', nargs='+')
        parser.add_argument('-o', '--output', help='Output file', required=True)

        h = 'Prune truth to S samples per class'
        parser = subparser.add_parser('prune', parents=[cls.vparser()], help=h, formatter_class=dhf)
        parser.add_argument('truthfile', help='Truth map to prune')
        parser.add_argument('-s', '--samples', help='Number of samples (in 1K units) to extract', default=10, type=int)

        # Results
        h = 'Compare truth maps output map(s)'
        parser = subparser.add_parser('analyze', parents=[cls.vparser()], help=h, formatter_class=dhf)
        parser.add_argument('-d', '--datadir', help='Input data directory of classmaps', required=True)
        parser.add_argument('-t', '--truth', help='Truth image (not pruned)', required=True)

        return parser0


def main():
    Truth.main()
