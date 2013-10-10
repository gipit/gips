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

    //! Copy input file into output file
	GeoImage Copy(const GeoImage& Input, GeoImage& Output, UNITS units) {
	    if (Input.NumBands() != Output.NumBands()) throw;
        for (unsigned int i=0; i<Output.NumBands(); i++) Output[i].Copy(Input[i], units);
	    // Set colors
		Colors colors = Input.GetColors();
		for (unsigned int i=0;i<Input.NumBands();i++) {
			//std::cout << "Setting band " << i+1 << " to " << colors[i+1] << std::endl;
			Output.SetColor(colors[i+1], i+1);
		}
		// Copies if only one band
		Output.CopyColorTable(Input);
		return Output;
	}

	//! Copy input file into new output file
	GeoImage Copy(const GeoImage& image, string filename, UNITS units, GDALDataType datatype) {
	    // TODO: if not supplied base output datatype on units (REFLECTIVITY should be float, etc)
	    if (datatype == GDT_Unknown) datatype = image.DataType();
		GeoImage ImageOut(filename, image, datatype);
        return Copy(image,ImageOut,units);
	}

	//! Calculate Radiance
	GeoImage Rad(const GeoImage& img, string filename) {
	    if (img[0].Units() != "radiance") {
	        throw std::runtime_error("image not in radiance units");
	    }
        typedef float T;
        GeoImageIO<T> imgIO(img);
        GeoImageIO<T> imgoutIO(GeoImage(filename, img, GDT_Int16));
        imgoutIO.SetNoData(-32768); // TODO - set nodata option
        imgoutIO.SetGain(0.1);
        std::vector<bbox> Chunks = img.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        CImg<T> cimg;
        CImg<unsigned char> nodata;
        for (unsigned int b=0;b<img.NumBands();b++) {
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                cimg = imgIO[b].Read(*iChunk);
                nodata = imgIO[b].NoDataMask(*iChunk);
                // only if nodata not same between input and output images
                cimg_forXY(cimg,x,y) { if (nodata(x,y)) cimg(x,y) = imgoutIO[b].NoDataValue(); }
                imgoutIO[b].Write(cimg,*iChunk);
            }
        }
        return imgoutIO;
	}

    //! Calculate reflectance assuming read image is radiance
	GeoImage Ref(const GeoImage& img, string filename) {
	    if (img[0].Units() != "radiance") {
	        throw std::runtime_error("image not in radiance units (required to calculate reflectance)");
	    }
	    typedef float t;
        GeoImageIO<t> imgIO(img);
        GeoImageIO<t> imgoutIO(GeoImage(filename, img, GDT_Int16));
        imgoutIO.SetNoData(-32768); // TODO - set nodata option
        std::vector<bbox> Chunks = img.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        CImg<float> cimg;
        CImg<unsigned char> nodata;
        for (unsigned int b=0;b<img.NumBands();b++) {
            if (img[b].Thermal()) imgoutIO[b].SetGain(0.01); else imgoutIO[b].SetGain(0.0001);
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                cimg = imgIO[b].Ref(*iChunk);
                nodata = imgIO[b].NoDataMask(*iChunk);
                cimg_forXY(cimg,x,y) { if (nodata(x,y)) cimg(x,y) = imgoutIO[b].NoDataValue(); }
                imgoutIO[b].Write(cimg,*iChunk);
            }
        }
        return imgoutIO;
	}

    GeoImage RGB(const GeoImage& image, string filename) {
        GeoImageIO<unsigned char> ImageOut(GeoImage(filename, image, GDT_Byte));
        GeoImageIO<float> imageIO(image);
        ImageOut.SetNoData(0);
        CImg<float> stats, cimg;
        CImg<unsigned char> mask;
        std::vector<bbox> Chunks = ImageOut.Chunk();
        std::vector<bbox>::const_iterator iChunk;
        for (unsigned int b=0;b<image.NumBands();b++) {
            stats = imageIO[b].ComputeStats();
            float lo = std::max(stats(2) - 3*stats(3), stats(0)-1);
            float hi = std::min(stats(2) + 3*stats(3), stats(1));
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                cimg = imageIO[b].Read(*iChunk);
                mask = imageIO[b].NoDataMask(*iChunk);
                ((cimg-=lo)*=(255.0/(hi-lo))).max(0.0).min(255.0);
                cimg_forXY(cimg,x,y) { if (mask(x,y)) cimg(x,y) = ImageOut[b].NoDataValue(); }
                ImageOut[b].Write(CImg<unsigned char>().assign(cimg.round()),*iChunk);
            }
        }
        return ImageOut;
    }

    // Multiply together all permutations
    /*GeoImage Permutations(const GeoImage& image1, const GeoImage& image2, string filename) {
        GeoImageIO<float> imageout(GeoImage(filename, image1, GDT_Float32, image1.NumBands()*image2.NumBands()));
        CImg<float> cimgout;
        int bandout(0);
        CImg<unsigned char> mask;
        float nodataout = -32768;
        imageout.SetNoData(nodataout);

        std::vector<bbox>::const_iterator iChunk;
        for (unsigned int i1=0; i1<image1.NumBands(); i1++) {
            GeoRasterIO<float> img1(image1[i1]);
            std::vector<bbox> Chunks = img1.Chunk();
            for (unsigned int i2=0; i2<image2.NumBands(); i2++) {
                GeoRasterIO<float> img2(image2[i2]);
                imageout[bandout].SetDescription(img1.Description()+'-'+img2.Description());
                for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                    cimgout = img1.Read(*iChunk).mul(img2.Read(*iChunk));
                    mask = img1.NoDataMask(*iChunk)|=(img2.NoDataMask(*iChunk));
                    cimg_forXY(mask,x,y) if (mask(x,y)) cimgout(x,y) = nodataout;
                    imageout[bandout].Write(cimgout, *iChunk);
                }
            bandout++;
            }
        }
        return imageout;
    }*/

    //! Convert lo-high of index into probability
    GeoImage Index2Probability(const GeoImage& image, string filename, float min, float max) {
        // Need method of generating new GeoImage with GeoRaster template in
        int bandnum = 1;
        GeoImageIO<float> imagein(image);
        GeoImageIO<float> imageout(GeoImage(filename, image, GDT_Float32, 2));
        float nodatain = imagein[0].NoDataValue();
        float nodataout = -32768;
        imageout.SetNoData(nodataout);

        std::vector<bbox>::const_iterator iChunk;
        std::vector<bbox> Chunks = image.Chunk();
        CImg<float> cimgin, cimgout;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            cimgin = imagein[bandnum-1].Read(*iChunk);
            cimgout = (cimgin - min)/(max-min);
            cimgout.min(1.0).max(0.0);
            cimg_forXY(cimgin,x,y) if (cimgin(x,y) == nodatain) cimgout(x,y) = nodataout;
            imageout[0].Write(cimgout, *iChunk);
            cimg_for(cimgout,ptr,float) if (*ptr != nodataout) *ptr = 1.0 - *ptr;
            imageout[1].Write(cimgout, *iChunk);
        }
        return imageout;
    }

    //! Perform band math (hard coded subtraction)
    GeoImage BandMath(const GeoImage& image, string filename, int band1, int band2) {
        GeoImageIO<float> imagein(image);
        GeoImageIO<float> imageout(GeoImage(filename, image, GDT_Float32, 1));
        float nodataout = -32768;
        imageout.SetNoData(nodataout);
        std::vector<bbox>::const_iterator iChunk;
        std::vector<bbox> Chunks = image.Chunk();
        CImg<float> cimgout;
        CImg<unsigned char> mask;
        for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
            mask = imagein[band1-1].NoDataMask(*iChunk)|=(imagein[band2-1].NoDataMask(*iChunk));
            cimgout = imagein[band1-1].Read(*iChunk) - imagein[band2-1].Read(*iChunk);
            cimg_forXY(mask,x,y) if (mask(x,y)) cimgout(x,y) = nodataout;
            imageout[0].Write(cimgout,*iChunk);
        }
        return imageout;
    }

	//! Rewrite file (applying processing, masks, etc)
	/* GeoImage Process(const GeoImage& image) {
	    for (unsigned int i=0l i<Output.NumBands(); i++) {

	    }
	}
	// Apply a mask to existing file (where mask>0 change to NoDataValue)
	GeoImage Process(const GeoImage& image, GeoRaster& mask) {
	    image.AddMask(mask);
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
	}	*/


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

	/*CImg<double> SpectralCovariance(const GeoImage& image) {
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
	}*/
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
    GeoImage kmeans( const GeoImage& image, string filename, int classes, int iterations, float threshold ) {
        //if (Image.NumBands() < 2) throw GIP::Gexceptions::errInvalidParams("At least two bands must be supplied");
        if (Options::Verbose()) {
            cout << image.Basename() << " - k-means unsupervised classifier:" << endl
                << "  Classes = " << classes << endl
                << "  Iterations = " << iterations << endl
                << "  Pixel Change Threshold = " << threshold << "%" << endl;
        }
        // Calculate threshold in # of pixels
        threshold = threshold/100.0 * image.Size();

        GeoImageIO<float> img(image);
        // Create new output image
        GeoImageIO<unsigned char> imgout(GeoImage(filename, image, GDT_Byte, 1));

        // Get initial class estimates (uses random pixels)
        CImg<float> ClassMeans = img.GetPixelClasses(classes);

        int i;
        CImg<double> Pixel, C_img, DistanceToClass(classes), NumSamples(classes), ThisClass;
        CImg<unsigned char> C_imgout, C_mask;
        CImg<double> RunningTotal(classes,image.NumBands(),1,1,0);

        vector<bbox> Chunks = image.Chunk();
        vector<bbox>::const_iterator iChunk;
        int NumPixelChange, iteration=0;
        do {
            NumPixelChange = 0;
            for (i=0; i<classes; i++) NumSamples(i) = 0;
            if (Options::Verbose()) cout << "  Iteration " << iteration+1 << std::flush;

            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                C_img = img.Read(*iChunk);
                C_mask = img.NoDataMask(*iChunk);
                C_imgout = imgout[0].Read(*iChunk);

                CImg<double> stats;
                cimg_forXY(C_img,x,y) { // Loop through image
                    // Calculate distance between this pixel and all classes
                    if (!C_mask(x,y)) {
                        Pixel = C_img.get_crop(x,y,0,0,x,y,0,C_img.spectrum()-1).unroll('x');
                        cimg_forY(ClassMeans,cls) {
                            ThisClass = ClassMeans.get_row(cls);
                            DistanceToClass(cls) = (Pixel - ThisClass).dot(Pixel - ThisClass);
                        }
                        // Get closest distance and see if it's changed since last time
                        stats = DistanceToClass.get_stats();
                        if (C_imgout(x,y) != (stats(4)+1)) {
                            NumPixelChange++;
                            C_imgout(x,y) = stats(4)+1;
                        }
                        NumSamples(stats(4))++;
                        cimg_forY(RunningTotal,iband) RunningTotal(stats(4),iband) += Pixel(iband);
                    } else C_imgout(x,y) = 0;
                }
                imgout[0].Write(C_imgout,*iChunk);
                if (Options::Verbose()) cout << "." << std::flush;
            }

            // Calculate new Mean class vectors
            for (i=0; i<classes; i++) {
                if (NumSamples(i) > 0) {
                    cimg_forX(ClassMeans,x) {
                        ClassMeans(x,i) = RunningTotal(i,x)/NumSamples(i);
                        RunningTotal(i,x) = 0;
                    }
                    NumSamples(i) = 0;
                }
            }
            if (Options::Verbose()) cout << 100.0*((double)NumPixelChange/image.Size()) << "% pixels changed class" << endl;
            if (Options::Verbose()>1) cimg_printclasses(ClassMeans);
        } while ( (++iteration < iterations) && (NumPixelChange > threshold) );

        imgout[0].SetDescription("k-means");
        //imgout.GetGDALDataset()->FlushCache();
        return imgout;
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
	GeoImage FixBadPixels(GeoImage& image) {
		typedef float T;
		vector<bbox> Chunks = image.Chunk();
		vector<bbox>::const_iterator iChunk;
		for (unsigned int b=0;b<image.NumBands();b++) {
			GeoRasterIO<T> band(image[b]);
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> img = band.Read(*iChunk, RAW);
				T nodata = band.NoDataValue();
				cimg_forXY(img,x,y)	if ( std::isinf(img(x,y)) || std::isnan(img(x,y)) ) img(x,y) = nodata;
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

