#!/usr/bin/env python

import os
import sys
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
import gippy
import numpy


class SLD(object):

    def str(self):
        return minidom.parseString(ET.tostring(self.top)).toprettyxml(indent="   ")

    def __init__(self, name):
        self.name = name
        self.top = ET.Element('StyledLayerDescriptor')
        self.top.set('version', "1.0.0")
        self.top.set('xsi:schemaLocation', "http://www.opengis.net/sld StyledLayerDescriptor.xsd")
        self.top.set('xmlns', "http://www.opengis.net/sld")
        self.top.set('xmlns:ogc', "http://www.opengis.net/ogc")
        self.top.set('xmlns:xlink', "http://www.w3.org/1999/xlink")
        self.top.set('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance")

        l = ET.SubElement(self.top, 'NamedLayer')
        ET.SubElement(l, 'Name').text = name
        s = ET.SubElement(l, 'UserStyle')
        ET.SubElement(s, 'Title').text = name
        fts = ET.SubElement(s, 'FeatureTypeStyle')
        rule = ET.SubElement(fts, 'Rule')
        self.raster_symbolizer = ET.SubElement(rule, 'RasterSymbolizer')
        self.Channels()

    def Channels(self):
        ET.SubElement(self.raster_symbolizer, 'Opacity').text = "1.0"
        a = ET.SubElement(self.raster_symbolizer, 'ChannelSelection')
        # Hard coded for gray only now
        b = ET.SubElement(a, 'GrayChannel')
        ET.SubElement(b, 'SourceChannelName').text = "1"

    def AddColorMap(self, mode, entries=[]):
        self.colormap = ET.SubElement(self.raster_symbolizer, 'ColorMap', type=mode)
        for e in entries:
            ET.SubElement(self.colormap, 'ColorMapEntry', 
                color='#'+str(e['color']), opacity=str(e['opacity']), quantity=str(e['value']), label=str(e['label']))

    def AddNoData(self, val):
        ET.SubElement(self.colormap, 'ColorMapEntry', color='#000000', opacity="0", quantity=str(val), label='No Data')
        
    def write(self, filename):
        f = open(filename, 'w')
        f.write(self.str())
        f.close()

if __name__ == "__main__":
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='Create SLD file for a raster', formatter_class=dhf)
    parser0.add_argument('file', help='Image file to create SLD for')
    parser0.add_argument('-r', '--ramp', nargs='*', help='Create ramp color map')
    args = parser0.parse_args()

    img = gippy.GeoImage(args.file)
    stats = numpy.squeeze(img[0].Stats())

    if args.ramp is not None:
        if len(args.ramp) == 0:
            args.ramp = ['000000', 'FFFFFF']
        numsteps = len(args.ramp)
        interval = (numsteps-1)/(stats[1]-stats[0])
        entries = [{'color': args.ramp[s], 'value': s*interval+stats[0], 'opacity':1, 'label':s} for s in range(0, numsteps)]

    bname = os.path.basename(os.path.splitext(args.file)[0])

    sld = SLD(bname)
    sld.AddColorMap('ramp', entries)
    if img[0].NoData():
        sld.AddNoData(img[0].NoDataValue())

    sld.write(bname + '.sld')
    print 'Style written to %s.sld' % bname

"""
import argparse
import xmltodict
import struct

from pdb import set_trace

if __name__ == "__main__":
    hformat = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='Raster Attribute Table 2 SLD', formatter_class=hformat)

    parser0.add_argument('file', help='XML File with Raster Attribute Table')

    args = parser0.parse_args()

    #xml = ElementTree.parse(args.file)
    doc = open(args.file,"r")
    xml = xmltodict.parse(doc.read())
    doc.close()

    att = xml['PAMDataset']['PAMRasterBand']['GDALRasterAttributeTable']  

    flds = [str(a['Name']) for a in att['FieldDefn']]
    
    reds = [int(255*float(a['F'][3])) for a in att['Row']]
    greens = [int(255*float(a['F'][4])) for a in att['Row']]
    blues = [int(255*float(a['F'][5])) for a in att['Row']]
    values = [str(a['F'][1]) for a in att['Row']]
    labels = [str(a['F'][7]) for a in att['Row']]
    labels2 = [str(a['F'][8]) for a in att['Row']]
    labels3 = [str(a['F'][9]) for a in att['Row']]
    labels4 = [str(a['F'][10]) for a in att['Row']]

    for i in range(0,len(labels)):
        rgb = (reds[i], greens[i], blues[i])
        col = struct.pack('BBB', *rgb).encode('hex')
        #print labels[i]
        #print labels2[i]
        #print labels3[i]
        #print labels4[i]
        print '<sld:ColorMapEntry color="#%s" opacity="1.0" quantity="%s" label="%s" />' % (col, values[i], labels[i])
"""