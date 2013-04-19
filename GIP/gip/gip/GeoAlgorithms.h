/*
 * gip_GeoAlgorithms.h
 *
 *  Created on: Aug 26, 2011
 *      Author: mhanson
 */

#ifndef GIP_GEOALGORITHMS_H
#define GIP_GEOALGORITHMS_H

#include <gip/GeoImage.h>
#include <gip/GeoRaster.h>

namespace gip {

    //! Copy input raster band into output raster band
    GeoRaster Copy(const GeoRaster& Input, GeoRaster& Output, UNITS units);

    //! Copy input file into output file
    GeoImage Copy(const GeoImage& Input, GeoImage& Output, UNITS units=RAW);

	//! Copy input file into new output file
	GeoImage Copy(const GeoImage&, std::string, UNITS units=RAW, GDALDataType=GDT_Unknown);

	//! Create new file of standard indices: NDVI, EVI, LSWI, NDSI, BI
	GeoImage Indices(const GeoImage&, std::string, bool=true, bool=true, bool=true, bool=true, bool=true);

	//! Create new file with AutoCloud algorithm
	GeoImage AutoCloud(const GeoImage&, std::string, int=4000, float=0.2, float=14, float=0.2, int=20);

	//! Create new file with a basic cloud mask
	GeoImage Fmask(const GeoImage&, std::string);

	//! Spectral Matched Filter
	//GeoImage SMF(const GeoImage& image, std::string, CImg<double>);

    //! Create a mask of NoData values
	//GeoRaster NoDataMask(const GeoImage& image, std::string filename="");

	//! Replace all 'Inf' results with the bands NoData value
	GeoImage InfReplace(GeoImage&);

    // Replace all 'NaN' results with bands NoData value
	//GeoImage NoDataReplace(GeoImage&);

	//! Calculate spectral covariance
	//CImg<double> SpectralCovariance(const GeoImage&);

    //! Calculate spectral correlation
	//CImg<double> SpectralCorrelation(const GeoImage&, CImg<double> covariance=CImg<double>() );

	//int test(int =1, char** =defaultargv("test") );
	//int hdf2tiff(std::string filename);

} // namespace gip

#endif /* GIP_GEOALGORITHMS_H_ */
