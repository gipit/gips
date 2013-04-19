/*
 * gip_GeoAlgorithms.cpp
 *
 *  Created on: Aug 26, 2011
 *      Author: mhanson
 */

#include <gip/Options.h>
#include <gip/GeoAlgorithms.h>
#include <gip/GeoImageIO.h>
#include <gip/gip_CImg.h>

namespace gip {
    using std::string;
	using std::vector;
	using std::cout;
	using std::endl;
	typedef boost::geometry::model::d2::point_xy<float> point;
	typedef boost::geometry::model::box<point> bbox;

    //! Copy input raster band into output raster band
	GeoRaster Copy(const GeoRaster& Input, GeoRaster& Output, UNITS units) {
        switch (Output.DataType()) {
            case GDT_Byte: GeoRasterIO<unsigned char>(Output).Copy(Input, units);
                break;
            case GDT_UInt16: GeoRasterIO<unsigned short>(Output).Copy(Input, units);
                break;
            case GDT_Int16: GeoRasterIO<short>(Output).Copy(Input, units);
                break;
            case GDT_UInt32: GeoRasterIO<unsigned int>(Output).Copy(Input, units);
                break;
            case GDT_Int32: GeoRasterIO<int>(Output).Copy(Input, units);
                break;
            case GDT_Float32: GeoRasterIO<float>(Output).Copy(Input, units);
                break;
            case GDT_Float64: GeoRasterIO<double>(Output).Copy(Input, units);
                break;
            default: GeoRasterIO<unsigned char>(Output).Copy(Input, units);
        }
        return Output;
	}

    //! Copy input file into output file
	GeoImage Copy(const GeoImage& Input, GeoImage& Output, UNITS units) {
	    if (Input.NumBands() != Output.NumBands()) {
            throw;
	    }
	    for (unsigned int i=0; i<Output.NumBands(); i++) Copy(Input[i], Output[i], units);
	    // Set colors
		Colors colors = Input.GetColors();
		for (unsigned int i=0;i<Input.NumBands();i++) {
			//std::cout << "Setting band " << i+1 << " to " << colors[i+1] << std::endl;
			Output.SetColor(colors[i+1], i+1);
		}
		return Output;
	}

	//! Copy input file into new output file
	GeoImage Copy(const GeoImage& image, string filename, UNITS units, GDALDataType datatype) {
	    // TODO: if not supplied base output datatype on units (REFLECTIVITY should be float, etc)
	    if (datatype == GDT_Unknown) {
	        datatype = image.DataType();
	    }
		GeoImage ImageOut(filename, image, datatype);
        return Copy(image,ImageOut,units);
	}

	//! Create multi-band image of various indices calculated from input
	GeoImage Indices(const GeoImage& ImageIn, string filename, bool ndvi, bool evi, bool lswi, bool ndsi, bool bi) {
		int numbands(0);
		if (ndvi) numbands++;
		if (evi) numbands++;
		if (lswi) numbands++;
		if (ndsi) numbands++;
		if (bi) numbands++;
		//if (satvi) numbands++;
		if (numbands == 0) {
			std::cout << "No indices selected for calculation!" << std::endl;
			return ImageIn;
		}
		GeoImageIO<float> imgin(ImageIn);
		GeoImageIO<float> imgout(GeoImage(filename, imgin, GDT_Float32, numbands));
		// Assume NoData is set!
		float nodatain = imgin[0].NoDataValue();
		float nodataout = -32768;
		//if (imgin[0].NoData())
		imgout.SetNoData(nodataout);
		// Main algorithm
		CImg<float> red, nir, blue, swir1, green, swir2, out;

        // need to add overlap
        std::vector<bbox> Chunks = imgout.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            red = imgin["Red"].Read(*iChunk, REFLECTIVITY);
            green = imgin["Green"].Read(*iChunk, REFLECTIVITY);
            blue = imgin["Blue"].Read(*iChunk, REFLECTIVITY);
            nir = imgin["NIR"].Read(*iChunk, REFLECTIVITY);
            swir1 = imgin["SWIR1"].Read(*iChunk, REFLECTIVITY);
            swir2 = imgin["SWIR2"].Read(*iChunk, REFLECTIVITY);
            int currentband(0);
            if (ndvi) {
                out = (nir-red).div(nir+red);
                cimg_forXY(out,x,y) if (red(x,y) == nodatain || nir(x,y) == nodatain) out(x,y) = nodataout;
                imgout[currentband++].Write(out,*iChunk);
            }
            if (evi) {
                out = 2.5*(nir-red).div(nir + 6*red - 7.5*blue + 1);
                cimg_forXY(out,x,y) if (red(x,y) == nodatain || nir(x,y) == nodatain || blue(x,y) == nodatain) out(x,y) = nodataout;
                imgout[currentband++].Write(out,*iChunk);
            }
            if (lswi) {
                out = (nir-swir1).div(nir+swir1);
                cimg_forXY(out,x,y) if (red(x,y) == nodatain || nir(x,y) == nodatain) out(x,y) = nodataout;
                imgout[currentband++].Write(out,*iChunk);
            }
            if (ndsi) {
                out = (swir1-green).div(swir1+green);
                cimg_forXY(out,x,y) if (red(x,y) == nodatain || nir(x,y) == nodatain) out(x,y) = nodataout;
                imgout[currentband++].Write(out,*iChunk);
            }
            if (bi) {
                out = 0.5*(blue+nir);
                cimg_forXY(out,x,y) if (blue(x,y) == nodatain || nir(x,y) == nodatain) out(x,y) = nodataout;
                imgout[currentband++].Write(out,*iChunk);
            }
            /*if (satvi) {
                out = (1.1 * ;
                cimg_forXY(out,x,y) if (red(x,y) == nodatain || nir(x,y) == nodatain) out(x,y) = nodataval;
                imgout[currentband++].Write(out,*iChunk);
            }*/
        }
        // Set descriptions
        int currentband(0);
        if (ndvi) imgout[currentband++].SetDescription("NDVI");
        if (evi) imgout[currentband++].SetDescription("EVI");
        if (lswi) imgout[currentband++].SetDescription("LSWI");
        if (ndsi) imgout[currentband++].SetDescription("NDSI");
        //if (satvi) imgout[currentband++].SetDescription("SATVI");
        if (bi) imgout[currentband++].SetDescription("BI");
        return imgout;
	}

	//! Auto cloud mask
	GeoImage AutoCloud(const GeoImage& image, string filename, int cheight, float minred, float maxtemp, float maxndvi, int morph) {
	    typedef float outtype;
        GeoImageIO<float> imgin(image);
        GeoImageIO<outtype> imgout(GeoImage(filename, image, GDT_Byte, 1));
        imgout.SetNoData(0);

        CImg<float> red, temp, ndvi;
        CImg<outtype> mask;

        // need to add overlap
        std::vector<bbox> Chunks = image.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {

            red = imgin["Red"].Read(*iChunk, REFLECTIVITY);
            temp = imgin["LWIR"].Read(*iChunk, REFLECTIVITY);
            ndvi = imgin.NDVI(*iChunk);

            mask =
                temp.threshold(maxtemp,false,true)^=1 &
                ndvi.threshold(maxndvi,false,true)^=1 &
                red.get_threshold(minred);

            if (morph != 0) mask.dilate(morph);

            imgout[0].Write(mask,*iChunk);
            /*CImg<double> stats = img.get_stats();
            cout << "stats " << endl;
            for (int i=0;i<12;i++) cout << stats(i) << " ";
            cout << endl;*/
        }
        return imgout;
	}

    //! Basic cloud mask
    GeoImage BasicCloudMask(const GeoImage& image, string filename) {
        GeoImageIO<float> imgin(image);
        GeoImageIO<unsigned char> imgout(GeoImage(filename,image,GDT_Byte,1));
        //imgout.SetNoData(0);

        CImg<unsigned char> cloudmask;
        std::vector<bbox> Chunks = image.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            cloudmask = imgin.BasicCloudMask(*iChunk);
            //nodatamask = imgin.NoDataMask(*iChunk);
            //cimg_forXY(cloudmask,x,y) if (nodatamask(x,y))

            imgout[0].Write(cloudmask,*iChunk);
            /*CImg<double> stats = img.get_stats();
            cout << "stats " << endl;
            for (int i=0;i<12;i++) cout << stats(i) << " ";
            cout << endl;*/
        }
        return imgout;
    }

	//! Cloud masking algorithm
	/*GeoImage fmask(const GeoImage& image) {
        //GeoRaster& thermal, const GeoRaster& red, const GeoRaster& ndvi) {
        GeoImageIO<unsigned char> mask(filename,thermal,GDT_Byte,1);

        std::vector<bbox> Chunks = image.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            CImg<float> img = src.Read(*iChunk);

            CImg<unsigned char> snow_mask, water_mask;

            CImg<T> imgout;
            //CImg<unsigned char> mask(img.thresh)
            if (Gain() != 1.0 || Offset() != 0.0)
                imgout = (img-Offset()) / Gain();
            else imgout = img;
            if (src.NoDataValue() != NoDataValue()) {
                cimg_forXY(img,x,y) { if (img(x,y) == src.NoDataValue()) imgout(x,y) = NoDataValue(); }
            }
            Write(imgout,*iChunk);
            //CImg<double> stats = img.get_stats();
            //cout << "stats " << endl;
            //for (int i=0;i<12;i++) cout << stats(i) << " ";
            //cout << endl;
        }
	}*/

	//! Spectral Matched Filter, with missing data
	/*GeoImage SMF(const GeoImage& image, string filename, CImg<double> Signature) {
		GeoImage output(filename, image, GDT_Float32, 1);

		// Band Means
		CImg<double> means(image.NumBands());
		for (unsigned int b=0;b<image.NumBands();b++) means(b) = image[b].Mean();

		//vector< box<point> > Chunks = ImageIn.Chunk();
		return output;
	}*/

	/*CImg<double> SpectralCovariance(const GeoImage& image) {
		typedef double T;

		unsigned int NumBands(image.NumBands());
		CImg<double> Covariance(NumBands, NumBands);

		unsigned int b;
		vector<GeoRasterIO<T> > bands;
		for (unsigned int b=0;b<image.NumBands();b++) {
			bands.push_back( GeoRasterIO<T>(image[b]) );
		}

		GeoRasterIO<unsigned char> mask( NoDataMask(image) );

		// Calculate Covariance
		//double TotalPixels(image.Size());
		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		//cout << "chunks = " << Chunks.size() << endl;
		CImg<T> bandchunk;
		CImg<unsigned char> maskchunk;
		for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
			int chunksize = boost::geometry::area(*iChunk);
			//cout << "chunksz " << chunksize << endl;
			CImg<T> matrixchunk(NumBands, chunksize);
			maskchunk = mask.Read(*iChunk);

			int p(0);
			for (b=0;b<NumBands;b++) {
				//cout << "band" << b << endl;
				CImg<T> bandchunk( bands[b].Read(*iChunk) );
				p = 0;
				cimg_forXY(bandchunk,x,y) {
					if (maskchunk(x,y)) matrixchunk(b,p++) = bandchunk(x,y);
				}
				//cout << "p = " << matrixchunk[p-1] << endl;
			}
			if (p != (int)image.Size()) matrixchunk.crop(0,0,NumBands-1,p-1);
			Covariance += (matrixchunk.get_transpose() * matrixchunk)/(mask.ValidSize()-1);
		}
		//cout << "done cov" << endl;
		// Subtract Mean
		CImg<double> means(NumBands);
		for (b=0; b<NumBands; b++) means(b) = image[b].Mean(); //cout << "Mean b" << b << " = " << means(b) << endl; }
		Covariance -= (means.get_transpose() * means);

		if (Options::Verbose() > 0) {
			cout << image.Basename() << " Spectral Covariance Matrix:" << endl;
			cimg_forY(Covariance,y) {
				cout << "\t";
				cimg_forX(Covariance,x) {
					cout << std::setw(18) << Covariance(x,y);
				}
				cout << endl;
			}
		}
		return Covariance;
	}

	CImg<double> SpectralCorrelation(const GeoImage& image, CImg<double> covariance) {
		// Correlation matrix
		if (covariance.size() == 0) covariance = SpectralCovariance(image);

		unsigned int NumBands = image.NumBands();
		unsigned int b;

		// Subtract Mean
		//CImg<double> means(NumBands);
		//for (b=0; b<NumBands; b++) means(b) = image[b].Mean();
		//covariance -= (means.get_transpose() * means);

		CImg<double> stddev(NumBands);
		for (b=0; b<NumBands; b++) stddev(b) = image[b].StdDev();
		CImg<double> Correlation = covariance.div(stddev.get_transpose() * stddev);

		if (Options::Verbose() > 0) {
			cout << image.Basename() << " Spectral Correlation Matrix:" << endl;
			cimg_forY(Correlation,y) {
				cout << "\t";
				cimg_forX(Correlation,x) {
					cout << std::setw(18) << Correlation(x,y);
				}
				cout << endl;
			}
		}

		return Correlation;
	}*/

	/*GeoRaster NoDataMask(const GeoImage& image, string filename) {
		typedef double T;
		if (filename == "") filename = image.Filename() + "_NoDataMask";
		GeoImage mask(filename, image, GDT_Byte, 1);
		GeoRasterIO<unsigned char> mask0(mask[0]);

		vector<GeoRasterIO<T> > Pbands;
		for (unsigned int b=0;b<image.NumBands();b++)
			Pbands.push_back( GeoRasterIO<T>(image[b]) );

		// Iterate through chunks
		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		CImg<T> imgchunk;
		CImg<unsigned char> imgout;
		unsigned int validpixels(0);
		for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
			// Read in all bands and make mask
			imgchunk = Pbands[0].Read(*iChunk);
			imgout = Pbands[0].NoDataMask(imgchunk);
			for (unsigned int b=1;b<image.NumBands();b++) {
				imgout &= Pbands[b].NoDataMask( Pbands[b].Read(*iChunk) );
			}
			validpixels += imgout.sum();
			mask0.Write(imgout, *iChunk);
		}
		mask[0].SetValidSize(validpixels);
		return mask[0];
	}*/

	//CImg<double> SpectralCorrelation(const GeoImage& image) {
	//}

	//! Replaces all NoData with NaN
	/*GeoImage NoDataReplace(GeoImage& image) {
		typedef float T;
		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		for (unsigned int b=0;b<image.NumBands();b++) {
			if (!image[b].NoData()) continue;
			double val = image[b].NoDataValue();
			GeoRasterIO<T> band(image[b]);
			string origunits = band.Units();
            // Retrieve Raw units
            band.SetUnits("Raw");
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> img = band.Read(*iChunk);
				cimg_for(img,ptr,T) if (*ptr == val) *ptr = NAN;
				band.Write(img,*iChunk);
			}
			band.ClearNoData();
			band.SetUnits(origunits);
		}
		return image;
	}*/

	//! Replaces all Inf/-Inf pixels with NoDataValue
	GeoImage InfReplace(GeoImage& image) {
		typedef float T;
		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		for (unsigned int b=0;b<image.NumBands();b++) {
			GeoRasterIO<T> band(image[b]);
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> img = band.Read(*iChunk, RAW);
				T nodata = band.NoDataValue();
				cimg_forXY(img,x,y)	if (std::isinf(img(x,y))) img(x,y) = nodata;
				band.Write(img,*iChunk);
			}
		}
		return image;
	}

	/*char** defaultargv(const char* ="");
	char** defaultargv(const char* name) {
		char** argv = new char*[2];
		argv[0] = new char[strlen(name)];
		strcpy(argv[0],name);
		argv[1] = NULL;
		return argv;
	}*/

	/*int test(int ac, char **argv) {
		cout << "ac = " << ac << endl;
		for (int i=0;i<ac;i++) {
			if (argv[i] != NULL) cout << argv[i] << endl; else ac=i;
		}
		cout << "ac = " << ac << endl;
		return 1;
	}*/

} // namespace gip

