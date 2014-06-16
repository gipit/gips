#!/usr/bin/env python

import argparse
import numpy
import gippy
from gipif.inventory import project_inventory
from gipif.utils import VerboseOut

__version__ = '0.7.0'


def main():
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(description='GIPIF Flood Detect', formatter_class=dhf)

    #group = parser0.add_argument_group('inventory arguments')
    parser0.add_argument('datadir', help='GIPIF Project directory', default='./')
    parser0.add_argument('-v', '--verbose', help='Verbosity - 0: quiet, 1: normal, 2: debug', default=1, type=int)
    parser0.add_argument('-p', '--product', help='Product to operate on', required=True)

    group = parser0.add_argument_group('algorithm options')
    group.add_argument('-o', '--output', help='Output file name', default='flood.tif')
    group.add_argument('--water', help='Threshold for water', default=-12.0, type=float)
    group.add_argument('--land', help='Threshold for land', default=-10.0, type=float)

    args = parser0.parse_args()
    gippy.Options.SetVerbose(args.verbose)

    # inventory project directory
    inv = project_inventory(args.datadir)

    VerboseOut('Flood Detect v%s' % __version__)
    VerboseOut('Processing %s dates' % len(inv))
    th0 = args.water
    th1 = args.land

    # Input image(s)
    dates = sorted(inv.keys())
    filenames = [inv[date][args.product] for date in dates]
    img = gippy.GeoImage(filenames)
    nodata = img[0].NoDataValue()

    # Flood detect algorithm
    imgout = gippy.GeoImage(args.output, img, gippy.GDT_Byte, img.NumBands()+2)
    imgout[0].SetDescription('hits')
    imgout[1].SetDescription('num_observations')
    for b in range(0, img.NumBands()):
        imgout[b+2].SetDescription(str(dates[b]))
    imgout.SetNoData(0)
    VerboseOut('Processing %s' % dates[0], 2)
    if th0 > th1:
        th0 = -th0
        th1 = -th1
        for b in range(0, img.NumBands()):
            img[b] = img[b] * -1

    last_water = (img[0] < th0).Read()
    last_water[numpy.where(last_water == nodata)] = 0
    hits = numpy.zeros(last_water.shape)
    numobs = numpy.zeros(last_water.shape)
    days = last_water
    imgout[2].Write(days)
    for b in range(1, img.NumBands()):
        VerboseOut('Processing %s' % dates[b], 2)
        days = days + (last_water * (dates[b]-dates[0]).days)
        # Increment # of observations
        numobs = numobs + img[b].DataMask()
        # Current dry land mask
        dried = (img[b] > th1).Read()
        dried[numpy.where(dried == nodata)] = 0
        # Find transitional pixels (water->dry) and increment
        inds = numpy.where(numpy.logical_and(dried, last_water))
        hits[inds] = hits[inds] + 1
        #imgout[b+1].Write(hits)
        # update water mask
        water = (img[b] < th0).Read()
        water[numpy.where(water == nodata)] = 0
        # add new water regions (1 day)
        nowater_inds = numpy.where(days == 0)
        days[nowater_inds] = water[nowater_inds]
        last_water = numpy.minimum(last_water + water, 1)
        # Reset water mask if it was dry
        dry_inds = numpy.where(dried > 0)
        last_water[dry_inds] = 0
        days[dry_inds] = 0
        imgout[b+2].Write(days)
    imgout[0].Write(hits)
    imgout[1].Write(numobs)

    imgout = None
    img = None

if __name__ == "__main__":
    main()

"""
    //! Rice detection algorithm
    GeoImage RiceDetect(const GeoImage& image, string filename, vector<int> days, float th0, float th1, int dth0, int dth1) {
        if (Options::Verbose() > 1) cout << "RiceDetect(" << image.Basename() << ") -> " << filename << endl;

        GeoImage imgout(filename, image, GDT_Byte, image.NumBands());
        imgout.SetNoData(0);
        imgout[0].SetDescription("rice");
        for (unsigned int b=1;b<image.NumBands();b++) {
            imgout[b].SetDescription("day"+to_string(days[b]));
        }

        CImg<float> cimg;
        CImg<unsigned char> cimg_datamask, cimg_dmask;
        CImg<int> cimg_th0, cimg_th0_prev, cimg_flood, cimg_flood_prev;
        int delta_day;

        for (unsigned int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
            if (Options::Verbose() > 3) cout << "Chunk " << iChunk << " of " << image[0].NumChunks() << endl;
            cimg = image[0].Read<float>(iChunk);
            cimg_datamask = image[0].DataMask(iChunk);
            CImg<int> DOY(cimg.width(), cimg.height(), 1, 1, 0);
            CImg<int> cimg_rice(cimg.width(), cimg.height(), 1, 1, 0);
            cimg_th0_prev = cimg.get_threshold(th0);
            cimg_flood = (cimg_th0_prev^1).mul(cimg_datamask);

            for (unsigned int b=1;b<image.NumBands();b++) {
                if (Options::Verbose() > 3) cout << "Day " << days[b] << endl;
                delta_day = days[b]-days[b-1];
                cimg = image[b].Read<float>(iChunk);
                cimg_datamask = image[b].DataMask(iChunk);
                cimg_th0 = cimg.get_threshold(th0);
                // Replace any nodata values with the last value
                cimg_forXY(cimg_datamask, x, y) {
                    if (cimg_datamask(x,y) == 0) { cimg_th0(x,y) = cimg_th0_prev(x,y); }
                }

                DOY += delta_day;                                       // running total of days
                DOY.mul(cimg_flood);                                    // reset if it hasn't been flooded yet
                DOY.mul(cimg_th0);                                      // reset if in hydroperiod

                cimg_dmask = DOY.get_threshold(dth1,false,true)^=1;      // mask of where past high date
                DOY.mul(cimg_dmask);

                // locate (and count) where rice criteria met
                CImg<unsigned char> newrice = cimg.threshold(th1,false,true) & DOY.get_threshold(dth0,false,true);
                cimg_rice = cimg_rice + newrice;

                // update flood map
                cimg_flood |= (cimg_th0^1);
                // remove new found rice pixels, and past high date
                cimg_flood.mul(newrice^=1).mul(cimg_dmask);

                //imgout[b].Write(DOY, iChunk);
                imgout[b].Write(DOY, iChunk);
                cimg_th0_prev = cimg_th0;
            }
            imgout[0].Write(cimg_rice,iChunk);              // rice map count
        }
        return imgout;
    }
"""
