#!/usr/bin/env python

import os
import sys
import argparse
from datetime import datetime
import agspy.utils.dateparse as dateparse
import pandas, numpy
import gippy
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix


# temp imports
import glob
from pdb import set_trace
from scipy.misc import toimage
import gc
import warnings

from gips.inventory import project_inventory

__version__ = '0.1.0'

#from skll.metrics import kappa

warnings.filterwarnings('ignore', category=DeprecationWarning, module='pandas', lineno=570)
#warnings.filterwarnings('ignore', category=UserWarning, module='sklearn', lineno=41)

def cmatrix(truth, pred, output=None, nodata=None):
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
    classes = numpy.unique(numpy.concatenate((truth,pred)))
    df = pandas.DataFrame(cm, index=classes, columns=classes)
    return {'cm': df, 'accuracy': accuracy, #'kappa': kappa,
        'prod_accuracy': prod_accuracy, 'user_accuracy': user_accuracy}

def calculate_stats(cmi): #ytrue, ypred):
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
        'producers': producers, 'consumers': consumers, 'npts': npts} #, 'cm': npts }

def VerboseOut(txt, level=1):
    if gippy.Options.Verbose() >= level: print txt

def classify(model, data, nodata=-32768):
    VerboseOut('Classifying', 2)
    #data = numpy.squeeze(gippy.GeoImage_float(img).Read())
    dims = data.shape
    data = data.reshape(dims[0], dims[1]*dims[2]).T
    valid = numpy.all(data != nodata, axis=1)
    classmap = numpy.zeros((dims[1]*dims[2]), numpy.byte)
    classmap[valid] = model.predict(data[valid,:]).astype(numpy.byte)
    classmap = classmap.reshape(dims[1],dims[2])
    return classmap

def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='Rice DSS Mapper', formatter_class=argparse.RawTextHelpFormatter)
    subparser = parser0.add_subparsers(dest='command')

    gparser = argparse.ArgumentParser(add_help=False, formatter_class=dhf)
    group = gparser.add_argument_group('General Options')
    group.add_argument('-v','--verbose', help='Verbosity level', default=1, type=int)

    parser = subparser.add_parser('truth', help='Make single truth map from multiple maps')
    parser.add_argument('files', help='List of truth files', nargs='+')

    parser = subparser.add_parser('train', help='Generate temporal vectors from data', 
        parents=[gparser], formatter_class=dhf)
    parser.add_argument('indir', help='Input data directory of training data')
    #parser.add_argument('product', help='Product to train on (will find matching files)')
    parser.add_argument('-d','--days', help='Beginning and end of interval (doy1,doy2)',default=(1,365))
    parser.add_argument('--max', help='Maximum number of samples', default=10000, type=int)
    parser.add_argument('-o','--output', help='Output directory', default='train')
    
    parser = subparser.add_parser('test', help='Use training data to test and output metrics', 
        parents=[gparser], formatter_class=dhf)
    parser.add_argument('-i','--input',help='Input training data directory',default='train')
    #parser.add_argument('-o','--output', help='Output directory', default='test')
    parser.add_argument('--start',help='Starting day of year', default=None, type=int)
    parser.add_argument('--step',help='Number of days per step', default=8, type=int)

    parser = subparser.add_parser('class', help='Classify', 
        parents=[gparser], formatter_class=dhf)
    #parser.add_argument('files', nargs='*', help='Image file to classify')
    parser.add_argument('-s','--site',help='Vector file for region of interest', default='site.shp')
    parser.add_argument('-d','--dates',help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    parser.add_argument('--days',help='Include data within these days of year (doy1,doy2)',default=None)
    parser.add_argument('-i','--input',help='Input training data directory',default='train')
    parser.add_argument('-o','--output', help='Output directory', default='classify')
    parser.add_argument('--series', help='Run for every day in data inventory', default=False, action='store_true')

    parser = subparser.add_parser('analyze', help='Compare directory of outputs vs truth', 
        parents=[gparser], formatter_class=dhf)
    parser.add_argument('-i','--input',help='Input data directory of classmaps',required=True)
    parser.add_argument('-t','--truth',help='Truth map to compare against',default=None)

    parser = subparser.add_parser('show', help='Display animation of classmaps vs time', 
        parents=[gparser], formatter_class=dhf)
    parser.add_argument('-i','--input', help='Input data directory of classmaps',required=True)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)
    gippy.Options.SetChunkSize(512)

    products = ['ndvi', 'lswi']

    VerboseOut('Temporal Classifier v%s' % __version__)

    if args.command == 'test' or args.command == 'class':
        # Read training data
        datafiles = [os.path.join(args.input, p+'.pkl') for p in products]
        dfs = {}
        for f in datafiles:
            basename = os.path.splitext(os.path.basename(f))[0]
            dfs[basename] = pandas.read_pickle(f)
            days = dfs[basename].columns
        truth = numpy.squeeze(pandas.read_pickle(os.path.join(args.input, 'truth.pkl')))
        classes = numpy.unique(truth.values).astype(int)
        cdl = gippy.GeoImage(os.path.join(args.input, 'cdl.tif'))

    if args.command == 'truth':
        # Merge multiple truth files into one showing agreement over all
        # Read in truth files from directory
        print 'Input files: ', args.files
        random.seed()
        VerboseOut('Merging truth files', 2)
        for i, filename in enumerate(glob.glob('CDL/CDL*mask.tif')):
            cdl = gippy.GeoImage(filename)
            if i == 0:
                cdlimg = gippy.GeoRaster_byte(cdl[0]).Read()
                cdlout = gippy.GeoImage(os.path.join(args.output, 'cdl'), cdl)
                cdlout.CopyColorTable(cdl)
                cdlout.SetNoData(0)
            else:
                img = gippy.GeoRaster_byte(cdl[0]).Read()
                it = numpy.nditer(cdlimg, flags=['multi_index'])
                while not it.finished:
                    (x, y) = it.multi_index
                    if cdlimg[x, y] != img[x, y]:
                        cdlimg[x, y] = 0
                    it.iternext()
                img = None
            cdl = None
        # Pruning down to max samples
        VerboseOut('Pruning truth map', 2)
        for i in range(1,numpy.max(cdlimg)+1):
            locs = numpy.where(cdlimg == i)
            num = len(locs[0])
            if num > args.max:
                sample = random.sample(range(0,num-1),num-args.max)
                locs = (locs[0][sample], locs[1][sample])
                cdlimg[locs] = 0
            elif num < 10:
                cdlimg[locs] = 0
            if num != 0:
                locs = numpy.where(cdlimg == i)
                VerboseOut('Class %s: %s -> %s samples' % (i, num, len(locs[0])), 2)
        gippy.GeoRaster_byte(cdlout[0]).Write(cdlimg)
        VerboseOut('CDL processed in %s' % (datetime.now() - start))        

    if args.command == 'train':
        start = datetime.now()
        inv = project_inventory(args.indir)
        VerboseOut('Extracting pixel temporal vectors from %s' % os.path.basename(args.indir), 2)

        import pdb
        pdb.set_trace()

        # Create day of year array TODO - auto from file list
        tmp = map(int, (args.days.split(',')))
        days = [d for d in range(tmp[0], tmp[1]+1)]                
        
        for p in products:
            filenames = {}
            ### REPLACE WITH INVENTORY
            for filename in glob.glob('training_data/*%s*tif' % p):
                basename = os.path.basename(filename)
                #year = int(basename[0:4])
                day = int(basename[4:7])
                filenames[day] = filename
                # what to do if more than one day?....average?
            day_index = sorted(filenames.keys())
            ########

            # Replace with gippy temporal read/open class (or factory function)
            img = gippy.GeoImage(filenames[day_index[0]])
            filenames.pop(day_index[0],None)
            for day in sorted(filenames.keys()):
                img.AddBand(gippy.GeoImage(filenames[day])[0])
            data = gippy.GeoImage_float(img).Extract(cdlout[0])

            truth = data[0,:]
            data = data[1:,:]

            df = pandas.DataFrame(data, index=day_index)
            df[df == -32768] = numpy.nan
            df = df.reindex(index=days)
            for col in df.columns: df[col] = df[col].interpolate()
            df.fillna(method='bfill',inplace=True)
            df.fillna(method='ffill',inplace=True)
            df.T.to_pickle( os.path.join(args.output,'%s.pkl' % p) )
            VerboseOut('Extracted %s training data in %s' % (p,datetime.now()-start),2)

        pandas.DataFrame(truth).to_pickle(os.path.join(args.output,'truth.pkl'))
        VerboseOut('Done training: %s' % (datetime.now()-start))

    elif args.command == 'test':
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

    elif args.command == 'class':
        from gippy.data.landsat import LandsatData
        maxclouds = 39
        if not os.path.exists(args.output): os.makedirs(args.output)
        inv = LandsatData.inventory(site=args.site, dates=args.dates, days=args.days, products=products, maxclouds=maxclouds)

        if args.series:
            date_series = [inv.dates[0].strftime('%Y-%j') + ',' + d.strftime('%Y-%j') for d in inv.dates]
        else:
            date_series = [args.dates]

        updatedir = os.path.join(args.output,'updates')
        if not os.path.exists(updatedir): os.makedirs(updatedir)
            
        for dates in date_series:
            start = datetime.now()
            inv = LandsatData.inventory(site=args.site, dates=dates, days=args.days, products=products, maxclouds=maxclouds)
            inv.project(datadir=os.path.join(args.output,'gipdata'))

            days = numpy.array([int(d.strftime('%j')) for d in inv.dates])
            #VerboseOut('Preparing data from %s dates between %s and %s' % (len(inv.dates), inv.dates[0],inv.dates[-1]))
            
            # extract days from training data
            i = 0
            for key in dfs:
                if i == 0:
                    data = dfs[key][days]
                else:
                    data = data.join(dfs[key][days],rsuffix=key)
                i = i+1

            classifier = RandomForestClassifier()
            classifier.fit(data, truth)

            # extract data from images - TODO - improve assembling time series
            imgarr = []
            for p in products:
                image = gippy.GeoImage_float(inv.get_timeseries(p))
                arr = numpy.squeeze(image.TimeSeries(days.astype('double')))
                if len(days) == 1:
                    dims = arr.shape
                    arr = arr.reshape(1,dims[0],dims[1])
                imgarr.append(arr)
                nodata = image[0].NoDataValue()
            data = numpy.vstack(tuple(imgarr))

            classmap = classify(classifier, data, nodata)

            imgout = gippy.GeoImage(os.path.join(updatedir,'classmap_day%s_update' % days[-1]), cdl, gippy.GDT_Byte, 1)
            imgout.SetNoData(0)
            imgout.CopyColorTable(cdl)
            gippy.GeoRaster_byte(imgout[0]).Write(classmap)
            
            if args.series:
                if len(days) == 1:
                    classmap_todate = classmap
                else:
                    locs = numpy.where(classmap != 0)
                    classmap_todate[locs] = classmap[locs]
                    imgout = None
                    imgout = gippy.GeoImage(os.path.join(args.output,'classmap_day%s' % days[-1]), cdl, gippy.GDT_Byte, 1)
                    imgout.SetNoData(0)
                    imgout.CopyColorTable(cdl)
                    gippy.GeoRaster_byte(imgout[0]).Write(classmap_todate)
                    imgout.CopyColorTable(cdl)
            
            inv = None
            image = None
            imgarr = None
            imgout = None
            data = None
            #gc.collect()

            VerboseOut('Classified in %s' % (datetime.now() - start))

    elif args.command == 'analyze':
        # compare directory of classmaps with truth map
        if args.truth is None:
            header = 'Day, Rice (ha)'
        else:
            cdl = gippy.GeoImage(args.truth)
            truth = gippy.GeoRaster_byte(cdl[0]).Read()
            header = 'Day, Overall Accuracy, Rice Kappa, Rice Accuracy, Rice Reliability, Rice (ha)'

        basename = os.path.basename(args.input)
        f = open(os.path.join(args.input,basename + '_analysis.csv'), 'w')
        f.write(header+'\n')
        classmaps = glob.glob(os.path.join(args.input, 'classmap*tif'))
        for classmap in classmaps:
            gimg = gippy.GeoImage(classmap)
            day = gimg.Basename()[12:15]
            pred = gippy.GeoRaster_byte(gimg[0]).Read()          
            # get total rice, in hectares
            ricetotal = len(numpy.where(pred == 3)[0]) * 0.09

            if args.truth is None:
                line = '%s, %4.2f' % (day, ricetotal)
            else:
                result = cmatrix(truth, pred, nodata=0)
                loc = numpy.where(result['cm'].columns == 3)

                cm = result['cm'].values
                ricecm = numpy.zeros((2,2))
                ricecm[1,1] = cm[loc,loc]
                ricecm[0,1] = (cm[loc,:].sum() - ricecm[1,1])
                ricecm[1,0] = (cm[:,loc].sum() - ricecm[1,1])
                ricecm[0,0] = cm.sum().sum() - ricecm.sum()
                rowsums = cm.sum(axis=0)
                colsums = cm.sum(axis=1)
                total = rowsums.sum().astype('float')
                err = (rowsums/total * colsums/total).sum()
                kappa = (numpy.diag(ricecm).sum()/total - err)/(1. - err)
                
                line = '%s, %4.2f, %4.2f, %4.2f, %4.2f, %8.2f' % \
                    (day, result['accuracy'], kappa, result['prod_accuracy'][loc], result['user_accuracy'][loc], ricetotal)
                f.write(line+'\n')
            print line

        f.close()

    elif args.command == 'show':
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


if __name__ == "__main__":
    main()
