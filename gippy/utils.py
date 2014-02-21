#!/usr/bin/env python

import os
import errno
import gippy
import datetime
import calendar

def VerboseOut(obj, level=1):
    if gippy.Options.Verbose() >= level: 
        #pprint.PrettyPrinter().pprint(obj)
        if not isinstance(obj,(list,tuple)): obj = [obj]
        for o in obj: print o

def File2List(filename):
    f = open(filename)
    txt = f.readlines()
    txt2 = []
    for t in txt: txt2.append( t.rstrip('\n') )
    return txt2

def List2File(lst,filename):
    f = open(filename,'w')
    f.write('\n'.join(lst)+'\n')
    f.close()

def RemoveFiles(filenames):
    for f in filenames:
        try:
            os.remove(f)
        except OSError as e:
            if e.errno != errno.ENOENT: raise
            continue

def _parse_date(dstring, last=False):
    """ Parses string of YYYY or YYYY-MM or YYYY-MM-DD or YYYY-DOY and returns date object """
    d = dstring.split('-')
    if len(d) == 2 and len(d[1]) == 3:
        dttmp = datetime.datetime(int(d[0]),1,1) + datetime.timedelta(days=int(d[1])-1)
        d[1] = dttmp.month
        d.append(dttmp.day)
    if ( not last ):
        if (len(d) == 1): d.append('1')
        if (len(d) == 2): d.append('1')
    else:
        if (len(d) == 1): d.append('12')
        if (len(d) == 2): d.append(calendar.monthrange(int(d[0]),int(d[1]))[1] )
    return datetime.date(int(d[0]),int(d[1]),int(d[2]))

def parse_dates(dstring):
    """ Parses string of 1 or 2 dates separated by a comma.  Valid format for each of the 2 dates:
        YYYY
        YYYY-MM
        YYYY-MM-DD
        YYYY-DOY
    """
    try:
        (d1,d2) = dstring.replace(',',' ').split()
        return (_parse_date(d1),_parse_date(d2,True))
    except:
        return (_parse_date(dstring),_parse_date(dstring,True))


"""
def http_fetch(obj, tile, date, dataset):

    httploc = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005' # only some datasets available here

    STAGE = '/titan/data/modis/stage'

    os.chdir(STAGE)

    year, month, day = date.timetuple()[:3]
    doy = date.timetuple[7]

    pattern = ''.join(['(', MDS, '.A', year, doy, '.', tile, '.005.\d{13}.hdf)'])
    mainurl = ''.join([httploc, '/', year, '.', '%02d'%month, '.', '%02d'%day])

    try:
        listing = urllib.urlopen(mainurl).readlines()
    except Exception, e:
        listing = None
        print 'unable to access %s' % mainurl
        continue

    cpattern = re.compile(pattern)
    name = None
    for item in listing:
        if cpattern.search(item):
            if 'xml' in item:
                continue
            print 'found in', item.strip()
            name = cpattern.findall(item)[0]
            print 'it is', name
            url = ''.join([mainurl, '/', name])
            print 'the url is', url
            try:
                urllib.urlretrieve(url, name)
                retrieved.append(name)
            except Exception, e:
                print 'unable to retrieve %s from %s' % (name, url)

    time.sleep(2)
"""
