#!/usr/bin/env python

import os
import argparse
from datetime import datetime
import numpy
import pandas
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

import gippy
from gips.core import Algorithm
from gips.utils import VerboseOut
from gips.inventory import ProjectInventory

# temp imports
import glob
from pdb import set_trace
from scipy.misc import toimage
import gc
import warnings
import traceback

#from skll.metrics import kappa

warnings.filterwarnings('ignore', category=DeprecationWarning, module='pandas', lineno=570)
#warnings.filterwarnings('ignore', category=UserWarning, module='sklearn', lineno=41)


class Tclass(Algorithm):
    name = 'Temporal Image Classifier'
    __version__ = '1.0.0'

    def __init__(self, truth, train, **kwargs):
        super(Tclass, self).__init__(**kwargs)
        if os.path.exists(truth):
            self.timg = gippy.GeoImage(truth)
        self.traindir = train.rstrip('/')

    def run(self, **kwargs):
        """ Train using truth, then classify """
        self.train(**kwargs)
        self.classify(**kwargs)

    def train(self, **kwargs):
        """ Extract training data from data given truth map """
        start = datetime.now()
        VerboseOut('Extracting training data to form spectral vectors', 2)
        tinv = ProjectInventory(self.traindir, self.inv.requested_products)
        days = [int(d.strftime('%j')) for d in tinv.dates]

        prefix = os.path.join(tinv.projdir, self.timg.Basename() + '.')

        for p in tinv.requested_products:
            fout = prefix + p + '.pkl'
            if os.path.exists(fout):
                continue
            # Extract data from training images and truth image
            img = tinv.get_timeseries(p)
            data = img.Extract(self.timg[0])
            truth = data[0, :]
            data = data[1:, :]
            data[data == img[0].NoDataValue()] = numpy.nan

            # Complete signature with interpolation and filling
            df = pandas.DataFrame(data, index=days)
            df.dropna(axis=1, how='all', inplace=True)
            truth = truth[df.columns.values]
            # aggregate identical days
            df = df.groupby(lambda x: x).mean()
            df = df.reindex(index=[d for d in range(min(days), max(days)+1)])
            VerboseOut('Interpolating and filling', 2)
            for col in df.columns:
                df[col] = df[col].interpolate()
            df.fillna(method='bfill', inplace=True)
            df.fillna(method='ffill', inplace=True)

            # Write output files
            VerboseOut('Writing out %s' % fout, 2)
            df.T.to_pickle(fout)
            img = None
            data = None
            df = None
        fout = prefix + 'truth.pkl'
        if not os.path.exists(fout):
            VerboseOut('Writing out %s' % fout, 2)
            pandas.DataFrame(truth).to_pickle(fout)
        VerboseOut('Done extracting training data in %s' % (datetime.now()-start))

    def _readtrain(self):
        """ Read training data """
        # Read training data
        prefix = os.path.join(self.traindir, self.timg.Basename())
        datafiles = ['%s.%s.pkl' % (prefix, p) for p in self.inv.requested_products]
        self.dfs = {}
        for f in datafiles:
            basename = os.path.splitext(os.path.basename(f))[0][len(self.timg.Basename())+1:]
            self.dfs[basename] = pandas.read_pickle(f)
            days = self.dfs[basename].columns
        self.truth = numpy.squeeze(pandas.read_pickle(prefix + '.truth.pkl'))
        set_trace()
        #self.classes = numpy.unique(truth.values).astype(int)

    def classify(self, series=False, **kwargs):
        """ Classify data given training data """

        # Create output directory and link to training data
        self.dout = '%s.tclass' % self.inv.projdir
        if not os.path.exists(self.dout):
            os.makedirs(self.dout)
        truth_link = os.path.abspath(os.path.join(self.dout, 'truth'))
        training_link = os.path.abspath(os.path.join(self.dout, 'train'))
        if os.path.exists(truth_link):
            os.remove(truth_link)
        if os.path.exists(training_link):
            os.remove(training_link)
        os.symlink(os.path.abspath(self.timg.Filename()), truth_link)
        os.symlink(os.path.abspath(self.traindir), training_link)

        if series:
            date_series = [self.inv.dates[0].strftime('%Y-%j') + ',' + d.strftime('%Y-%j') for d in self.inv.dates]
        else:
            date_series = [self.inv.dates]

        self._readtrain()
        self.series = series

        if series:
            days = [self.inv.dates[0]]
            self._train_and_classify(days)
            for d in self.inv.dates[1:]:
                days.append(d)
                self._train_and_classify(days)
        else:
            self._train_and_classify(self.inv.dates)

    def _train_and_classify(self, dates):
        start = datetime.now()
        VerboseOut('Classifying data from %s dates between %s and %s' % (len(dates), dates[0], dates[-1]), 2)
        days = numpy.array([int(d.strftime('%j')) for d in dates])

        # extract days from training data
        i = 0
        for p in self.inv.requested_products:
            tdata = self.dfs[p][days]
            if i == 0:
                tdata = self.dfs[p][days]
            else:
                tdata = tdata.join(self.dfs[p][days], rsuffix=p)
            i = i+1

        VerboseOut('Extracting data from images', 2)
        # extract data from images - TODO - improve assembling time series
        imgarr = []
        for p in self.inv.requested_products:
            gimg = self.inv.get_timeseries(p, dates=dates)
            # TODO - move numpy.squeeze into swig interface file?
            arr = numpy.squeeze(gimg.TimeSeries(days.astype('float64')))

            if len(days) == 1:
                dims = arr.shape
                arr = arr.reshape(1, dims[0], dims[1])
            imgarr.append(arr)
        data = numpy.vstack(tuple(imgarr))

        # train classifier for this many days
        VerboseOut('Training and running classifier', 2)
        classifier = RandomForestClassifier()
        classifier.fit(tdata, self.truth)
        classmap = self._classify(classifier, data, gimg[0].NoDataValue())

        # Output most recently added
        fout = os.path.join(self.dout, '%s' % str(days[-1]).zfill(3))
        VerboseOut('Writing classmap to %s' % fout)
        imgout = gippy.GeoImage(fout, gimg, gippy.GDT_Byte, 1)
        imgout.SetNoData(0)
        imgout.CopyColorTable(self.timg)
        try:
            imgout[0].Write(classmap)
        except Exception, e:
            VerboseOut("problem writing")
            VerboseOut(traceback.format_exc(), 3)
        imgout = None

        if self.series:
            if len(days) == 1:
                self.classmap_todate = classmap
            else:
                locs = numpy.where(classmap != 0)
                self.classmap_todate[locs] = classmap[locs]
                fout = os.path.join(self.dout, '%s_master' % str(days[-1]).zfill(3))
                VerboseOut('Writing classmap to %s' % fout)
                imgout = gippy.GeoImage(fout, gimg, gippy.GDT_Byte, 1)
                imgout.SetNoData(0)
                if self.colortable != '':
                    imgout.CopyColorTable(self.colortable)
                imgout[0].Write(self.classmap_todate)
        imgout = None
        VerboseOut('Trained and classified in %s' % (datetime.now() - start), 2)

    def _classify(self, model, data, nodata=-32768):
        """ Uses model to classify data """
        VerboseOut('Classifying', 2)
        #data = numpy.squeeze(gippy.GeoImage_float(img).Read())
        dims = data.shape
        data = data.reshape(dims[0], dims[1]*dims[2]).T
        valid = numpy.all(data != nodata, axis=1)
        classmap = numpy.zeros((dims[1]*dims[2]), numpy.uint8)
        classmap[valid] = model.predict(data[valid, :]).astype(numpy.uint8)
        classmap = classmap.reshape(dims[1], dims[2])
        return classmap

    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(add_help=False, parents=[cls.project_parser()])
        parser.add_argument('-t', '--truth', help='Truth image used with training data', required=True)
        parser.add_argument('--train', help='Directory of training data', required=True)
        parser.add_argument('--series', help='Run days as series', default=False, action='store_true')

        """
        parser = subparser.add_parser('show', help='Display animation of classmaps vs time',
            parents=[gparser], formatter_class=dhf)
        parser.add_argument('-i', '--input', help='Input data directory of classmaps', required=True)

        parser = subparser.add_parser('test', help='Use training data to test and output metrics',
            formatter_class=dhf)
        parser.add_argument('-i', '--input', help='Input training data directory', default='train')
        #parser.add_argument('-o','--output', help='Output directory', default='test')
        parser.add_argument('--start', help='Starting day of year', default=None, type=int)
        parser.add_argument('--step', help='Number of days per step', default=8, type=int)
        """
        return parser


def main():
    Tclass.main()

    """
    def test(self):
        from sklearn.cross_validation import train_test_split
        start = datetime.now()

        # Day vector
        if args.start is None: args.start = min(days)
        testdays = range(args.start,max(days),args.step)
        cur_testdays = []
        accuracy = numpy.zeros((len(classes),len(testdays)),float)

        outputdir = 'test_d%s_s%s' % (args.start, args.step)
        if not os.path.exists(outputdir): os.makedirs(outputdir)

        for iday,day in enumerate(testdays):
            cur_testdays.append(day)
            i = 0
            for key in dfs:
                if i == 0:
                    data = dfs[key][cur_testdays]
                else:
                    data = data.join(dfs[key][cur_testdays],rsuffix=key)
                i = i+1

            # Split into train and test
            X_train, X_test, y_train, y_test = train_test_split(data,truth,train_size=0.25)

            # Classify
            classifier = RandomForestClassifier()
            classifier.fit(X_train, numpy.squeeze(y_train))
            y_pred = classifier.predict(X_test)

            # Generate confusion matrix
            cm = pandas.DataFrame(confusion_matrix(y_test, y_pred), index=classes, columns=classes)
            totals = numpy.sum(cm.values,1)
            for i,c in enumerate(classes): accuracy[i,iday] = cm.values[i,i]/float(totals[i])
            cm.to_csv(os.path.join(outputdir,'confusion_matrix_%s.csv' % day))
            stats = calculate_stats(cm)
            f = open(os.path.join(outputdir,'confusion_matrix_%s_stats.csv' % day), 'w')
            f.write('kappa = %s\n' % stats['kappa'])
            f.write('overall = %s\n' % stats['overall'])
            f.write('maxmiss = %s\n' % stats['maxmiss'])
            f.write('producers = ')
            for item in stats['producers']: f.write(str(item))
            f.write('\n')
            f.write('consumers = ')
            for item in stats['consumers']: f.write(str(item))
            f.write('\n')
            f.write('npts = ')
            for item in stats['npts']: f.write(str(item))
            f.write('\n')
            f.close()
            VerboseOut("Day %s: %4.2f%%" % (day, cm.values[2,2]/float(totals[2]) * 100), 1)

        print 'Done: %s' % (datetime.now() - start)
        filename = os.path.join(outputdir,'accuracy.csv')
        pandas.DataFrame(accuracy, index=classes, columns=cur_testdays).to_csv(filename)
    """
    """
    def show(self):
        filenames = glob.glob(os.path.join(args.input,'classmap*tif'))
        for filename in filenames:
            gimg = gippy.GeoImage(filename)
            img = gippy.GeoRaster_byte(gimg[0]).Read()
            dim = img.shape
            ar = float(dim[0])/dim[1]
            #img.resize()
            set_trace()
            img = numpy.resize(img,(int(800*ar),800))
            toimage(img).show()
            set_trace()
    """