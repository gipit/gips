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
	//using namespace boost::accumulators;

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

    //! Fmask cloud mask
    GeoImage Fmask(const GeoImage& image, string filename, int tolerance, int dilate) {
        GeoImageIO<float> imgin(image);

        // Output masks
        GeoImageIO<unsigned char> imgout(GeoImage(filename,image,GDT_Byte,2));
        imgout.SetNoData(0);

        // Output probabilties (for debugging/analysis)
        GeoImageIO<float> probout(GeoImage(filename + "_prob", image, GDT_Float32, 1));
        float nodata = -32768;
        probout.SetNoData(nodata);

        CImg<unsigned char> mask, wmask, lmask, nodatamask, redsatmask, greensatmask;
        CImg<float> red, nir, green, swir1, swir2, BT, ndvi, ndsi, white, vprob;
        float _ndvi, _ndsi;
        //int cloudpixels(0);
        CImg<double> wstats(image.Size()), lstats(image.Size());
        int wloc(0), lloc(0);

        vector<bbox> Chunks = image.Chunk();
        vector<bbox>::const_iterator iChunk;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            red = imgin["Red"].Read(*iChunk, REFLECTIVITY);
            green = imgin["Green"].Read(*iChunk, REFLECTIVITY);
            nir = imgin["NIR"].Read(*iChunk, REFLECTIVITY);
            swir1 = imgin["SWIR1"].Read(*iChunk, REFLECTIVITY);
            swir2 = imgin["SWIR2"].Read(*iChunk, REFLECTIVITY);
            BT = imgin["LWIR"].Read(*iChunk, REFLECTIVITY);
            nodatamask = imgin.NoDataMask(*iChunk)^=1;
            // floodfill....seems bad way
            //shadowmask = nir.draw_fill(nir.width()/2,nir.height()/2,)
            ndvi = (nir-red).div(nir+red);
            ndsi = (green-swir1).div(green+swir1);
            redsatmask = imgin["Red"].SaturationMask(*iChunk);
            greensatmask = imgin["Green"].SaturationMask(*iChunk);
            white = imgin.Whiteness(*iChunk);
            vprob = red;

            // Calculate "variability probability"
            cimg_forXY(vprob,x,y) {
                _ndvi = (redsatmask(x,y) && nir(x,y) > red(x,y)) ? 0 : abs(ndvi(x,y));
                _ndsi = (greensatmask(x,y) && swir1(x,y) > green(x,y)) ? 0 : abs(ndsi(x,y));
                vprob(x,y) = 1 - std::max(white(x,y), std::max(_ndsi, _ndvi));
            }
            probout[0].Write(vprob, *iChunk);

            // Potential cloud layer
            mask =
                swir2.threshold(0.03)
                & BT.get_threshold(27,false,true)^=1
                // NDVI
                & ndvi.threshold(0.8,false,true)^=1
                // NDSI
                & ndsi.threshold(0.8,false,true)^=1
                & imgin.HazeMask(*iChunk)
                & imgin.Whiteness(*iChunk).threshold(0.7,false,true)^=1
                & nir.div(swir1).threshold(0.75);

            //cloudpixels += mask.sum();

            // Water and land masks
            wmask = imgin.WaterMask(*iChunk);
            imgout[0].Write(mask, *iChunk);
            imgout[1].Write(wmask, *iChunk);

            //lmask = wmask^=1 & mask^=1;
            lmask = (wmask^1) & (mask^=1) & nodatamask;
            wmask = wmask & swir2^=1 & nodatamask;

            cimg_forXY(BT,x,y) if (wmask(x,y)) wstats[wloc++] = BT(x,y);
            cimg_forXY(BT,x,y) if (lmask(x,y)) lstats[lloc++] = BT(x,y);

            /*lBT = BT;
            wBT = BT;
            cimg_forXY(BT,x,y) if (!wmask(x,y)) wBT(x,y) = nodata;
            cimg_forXY(BT,x,y) if (!mask(x,y)) lBT(x,y) = nodata;
            probout[0].Write(wBT,*iChunk);
            probout[1].Write(lBT,*iChunk);*/
        }

        // If not enough non-cloud pixels then return existing mask
        //if (cloudpixels >= (0.999*imgout[0].Size())) return imgout;

        // Calculate some thresholds
        double zlo(-0.9345);
        double zhi(0.9345);

        CImg<double> _wstats( wstats.crop(0,0,0,0, wloc-1,0,0,0).stats() );
        CImg<double> _lstats( lstats.crop(0,0,0,0, lloc-1,0,0,0).stats() );
        //CImg<double> _wstats = probout[0].ComputeStats();
        //CImg<double> _lstats = probout[1].ComputeStats();
        double wmean(_wstats(2)), wstddev(sqrt(_wstats(3)));
        double lmean(_lstats(2)), lstddev(sqrt(_lstats(3)));
        double Twater(zhi*wstddev + wmean);
        double Tlo(zlo*lstddev + lmean);
        double Thi(zhi*lstddev + lmean);

        if (Options::Verbose() > 0) {
            cout << "Running fmask in " << Chunks.size() << " chunks with tolerance of " << tolerance << endl;
            cout << "Temperature stats:" << endl;
            cout << "  Water: " << wmean << " s = " << wstddev << endl;
            cout << "  Water: " << lmean << " s = " << lstddev << endl;
            cout << "  Twater (82.5%) = " << Twater << endl;
            cout << "  Tlo (17.5%) = " << Tlo << endl;
            cout << "  Thi (82.5%) = " << Thi << endl;
        }

        // 2nd pass cloud probabilities over land
        CImg<float> wprob, lprob;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            // Cloud over Water prob
            BT = imgin["LWIR"].Read(*iChunk,REFLECTIVITY);
            nodatamask = imgin.NoDataMask(*iChunk);

            vprob = probout[0].Read(*iChunk);

            // temp probability x variability probability
            lprob = ((Thi + 4-BT)/=(Thi+4-(Tlo-4))).mul( vprob );
            //1 - imgin.NDVI(*iChunk).abs().max(imgin.NDSI(*iChunk).abs()).max(imgin.Whiteness(*iChunk).abs()) );

            cimg_forXY(nodatamask,x,y) if (nodatamask(x,y) == 1) lprob(x,y) = nodata;

            // Cloud probability over water
            probout[0].Write( lprob, *iChunk);
        }

        // Apply thresholds to get final max
        float tol = (tolerance-3)*0.1;
        float wthresh = 0.5 + tol;
        CImg<double> stats = probout[0].ComputeStats();
        float lthresh( zhi*stats(3) + stats(2) + 0.2 + tol );

        if (Options::Verbose() > 0) {
            std::cout << "Cloud probability thresholds:" << std::endl;
            std::cout << "  Over Water = " << wthresh << std::endl;
            std::cout << "  Over Land = " << lthresh << std::endl;
        }

        // 3x3 filter of 1's for majority filter
        CImg<int> filter(3,3,1,1, 1);
        int esize = 5;
        CImg<int> erode_elem(esize,esize,1,1,1);
        CImg<int> dilate_elem(esize+dilate,esize+dilate,1,1,1);

        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            mask = imgout[0].Read(*iChunk);
            wmask = imgout[1].Read(*iChunk);

            nodatamask = imgin.NoDataMask(*iChunk);
            BT = imgin["LWIR"].Read(*iChunk,REFLECTIVITY);
            swir1 = imgin["SWIR1"].Read(*iChunk,REFLECTIVITY);

            // temp probability x brightness probability
            wprob = ((Twater - BT)/=4.0).mul( swir1.min(0.11)/=0.11 );
            lprob = probout[0].Read(*iChunk);

            // Combine probabilities of land and water
            mask &= ( (wmask & wprob.threshold(wthresh)) |= ((wmask != 1) & lprob.get_threshold(lthresh)));
            // Add in obvious outliers
            mask |=
                (lprob.threshold(0.99) & (wmask != 1)) |=
                (BT.threshold(Tlo-25,false,true)^=1);

            // Majority filter
            mask.convolve(filter).threshold(5);

            // Erode, then dilate twice
            mask.erode(erode_elem).dilate(dilate_elem);
            mask.dilate(dilate_elem);

            cimg_forXY(nodatamask,x,y) if (nodatamask(x,y) == 1) mask(x,y) = 0;
            imgout[0].Write(mask, *iChunk);
        }

        // Add shadow mask: bitmask

        // NoData regions
        return imgout;
    }

	//! Spectral Matched Filter, with missing data
	/*GeoImage SMF(const GeoImage& image, string filename, CImg<double> Signature) {
		GeoImage output(filename, image, GDT_Float32, 1);

		// Band Means
		CImg<double> means(image.NumBands());
		for (unsigned int b=0;b<image.NumBands();b++) means(b) = image[b].Mean();

		//vector< box<point> > Chunks = ImageIn.Chunk();
		return output;
	}*/

	CImg<double> SpectralCovariance(const GeoImage& image) {
		typedef double T;

		GeoImageIO<T> img(image);

		unsigned int NumBands(image.NumBands());
		CImg<double> Covariance(NumBands, NumBands);

		// Calculate Covariance
		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		CImg<T> bandchunk;
		CImg<unsigned char> mask;
		for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
			int chunksize = boost::geometry::area(*iChunk);
			CImg<T> matrixchunk(NumBands, chunksize);
			mask = img.NoDataMask(*iChunk);
			int validsize = mask.size() - mask.sum();

			int p(0);
			for (unsigned int b=0;b<NumBands;b++) {
				cout << "band" << b << endl;
				CImg<T> bandchunk( img[b].Read(*iChunk) );
				p = 0;
				cimg_forXY(bandchunk,x,y) {
					if (mask(x,y)==0) matrixchunk(b,p++) = bandchunk(x,y);
				}
				//cout << "p = " << matrixchunk[p-1] << endl;
			}
			if (p != (int)image.Size()) matrixchunk.crop(0,0,NumBands-1,p-1);
			Covariance += (matrixchunk.get_transpose() * matrixchunk)/(validsize-1);
		}
		cout << "done cov" << endl;
		// Subtract Mean
		CImg<double> means(NumBands);
		for (unsigned int b=0; b<NumBands; b++) means(b) = image[b].Mean(); //cout << "Mean b" << b << " = " << means(b) << endl; }
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
/*
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


    //! k-means unsupervised classifier
    /*GeoImage kmeans( const GeoImage& image, string filename, int classes, int iterations, float threshold ) {
        //if (Image.NumBands() < 2) throw GIP::Gexceptions::errInvalidParams("At least two bands must be supplied");
        if (Options::Verbose()) {
            cout << image.Basename() << ": " << "k-means Unsupervised Classifier:" << endl
                << "\tClasses = " << classes << endl
                << "\tIterations = " << iterations << endl
                << "\tPixel Change Threshold = " << threshold << "%" << endl;
        }
        // Calculate threshold in # of pixels
        threshold = threshold/100.0 * image.XSize() * image.YSize();

        // Create new output image
        GeoImage kimg(filename, image, GDT_Byte);

        // Get initial class estimates (uses random pixels)
        CImg<float> ClassMeans = Images.GetPixelClasses(classes);

        CImg<float> Pixel, img, DistanceToClass(classes), NumSamples(classes,1,1,1,0);
        CImg<int> imgout;
        CImg<double> RunningTotal(classes+1,Images.NumBands(),1,1,0);

        int NumPixelChange, iteration=0, PercentComplete=0;
        if (Images.GetOptions().verbose) cout << "0% complete" << StatusNewLine;
        do {
            NumPixelChange = 0;
            vector<GIP::Rect<> > Chunks = NewImage.Chunk();
            vector<GIP::Rect<> >::iterator iChunk;
            int chunknum = 0;
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                img = Images.GetRegion(*iChunk);
                imgout = NewImage.GetRegion(*iChunk);
                CImg<T> stats;
                cimg_forXY(img,x,y) {
                    // Calculate distance between this pixel and all classes
                    Pixel = img.get_crop(x,y,0,0,x,y,0,img.spectrum()-1).unroll('x');
                    cimg_forY(ClassMeans,C) {
                        CImg<T> ThisClass = ClassMeans.get_line(C);
                        DistanceToClass(C) = (Pixel - ThisClass).dot(Pixel - ThisClass);
                    }
                    // Get closest distance
                    stats = DistanceToClass.get_stats();
                    // Check if new class is same as old class
                    if (imgout(x,y) != (stats(4)+1)) NumPixelChange++;
                    imgout(x,y) = stats(4)+1;
                    NumSamples(stats(4))++;
                    cimg_forY(RunningTotal,yband) RunningTotal(stats(4)+1,yband) += Pixel(yband);
                }
                NewImage.WriteRegion(imgout,*iChunk);
                PercentComplete += (100/Chunks.size())/iterations; //*(iteration+1)/iterations) + (++chunknum*100/Chunks.size());
                if (Images.GetOptions().verbose) cout << "\r" << PercentComplete << "% complete" << StatusNewLine;
            }
            //if (Images.GetOptions().verbose) cout << endl;
            // Calculate new Mean class vectors
            for (int i=0; i<classes; i++) {
                if (NumSamples(i) > 0) {
                    cimg_forX(ClassMeans,x) {
                        ClassMeans(x,i) = RunningTotal(i+1,x)/NumSamples(i);
                        RunningTotal(i+1,x) = 0;
                    }
                    NumSamples(i) = 0;
                }
                // Output Class Vectors
                //if (Images.GetOptions().verbose) {
                 //   cout << "Class " << i << " vector: ";
                  //  cimg_forX(ClassMeans,x) cout << ClassMeans(x,i) << "  ";
                   // cout << endl;
                //}
            }
            PercentComplete = (100*(iteration+1)/iterations);
            if (Images.GetOptions().verbose) {
                cout << " Iteration " << iteration+1 << ": " << 100.0*NumPixelChange/Images.Size() << "% pixels changed classes" << endl;
                cout << PercentComplete << "% complete" << StatusNewLine;
            }
        } while ( (++iteration < iterations) && (NumPixelChange > threshold) );

        NewImage.SetDescription("k-means");
        NewImage.SetMetadata("color", "gray");
        NewImage.GetGDALDataset()->FlushCache();
        if (Images.GetOptions().verbose) cout << endl << Images.Basename() << ": k-means end" << endl;
        return NewImage;
    }*/

	//! Apply a mask to existing file (where mask>0 change to NoDataValue)
	GeoImage ApplyMask(const GeoImage& image, GeoRaster& mask) {
        for (unsigned int i=0; i<image.NumBands(); i++) {
            switch (image.DataType()) {
                case GDT_Byte: GeoRasterIO<unsigned char>(image[i]).ApplyMask(mask);
                    break;
                case GDT_UInt16: GeoRasterIO<unsigned short>(image[i]).ApplyMask(mask);
                    break;
                case GDT_Int16: GeoRasterIO<short>(image[i]).ApplyMask(mask);
                    break;
                case GDT_UInt32: GeoRasterIO<unsigned int>(image[i]).ApplyMask(mask);
                    break;
                case GDT_Int32: GeoRasterIO<int>(image[i]).ApplyMask(mask);
                    break;
                case GDT_Float32: GeoRasterIO<float>(image[i]).ApplyMask(mask);
                    break;
                case GDT_Float64: GeoRasterIO<double>(image[i]).ApplyMask(mask);
                    break;
                default: GeoRasterIO<unsigned char>(image[i]).ApplyMask(mask);
            }
        }
        return image;
	}

    //! Create mask based on NoData values for all bands
	GeoRaster CreateMask(const GeoImage& image, string filename) {
		typedef float T;

		GeoImageIO<T> imageIO(image);
		CImg<T> imgchunk;

		if (filename == "") filename = image.Basename() + "_mask";
		GeoImage mask(filename, image, GDT_Byte, 1);
		GeoRasterIO<unsigned char> mask0(mask[0]);
		CImg<unsigned char> imgout;
		mask.SetNoData(0);
		//unsigned int validpixels(0);

		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
			//imgchunk = imageIO[0].Read(*iChunk);
			//imgout = Pbands[0].NoDataMask(imgchunk);
			//for (unsigned int b=1;b<image.NumBands();b++) {
            //        imgout &= Pbands[b].NoDataMask( Pbands[b].Read(*iChunk) );
			//}
			//validpixels += imgout.sum();
			imgout = imageIO.NoDataMask(*iChunk)^=1;
			mask0.Write(imgout,*iChunk);
		}
		//mask[0].SetValidSize(validpixels);
		return mask[0];
	}

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

