/*
 * gip_GeoAlgorithms.cpp
 *
 *  Created on: Aug 26, 2011
 *      Author: mhanson
 */

#include <gip/GeoAlgorithms.h>
#include <gip/GeoImageIO.h>
#include <gip/gip_CImg.h>

#include <gdal/ogrsf_frmts.h>
#include <gdal/gdalwarper.h>

class CutlineTransformer : public OGRCoordinateTransformation
{
public:

    void         *hSrcImageTransformer;

    virtual OGRSpatialReference *GetSourceCS() { return NULL; }
    virtual OGRSpatialReference *GetTargetCS() { return NULL; }

    virtual int Transform( int nCount, double *x, double *y, double *z = NULL ) {
        int nResult;

        int *pabSuccess = (int *) CPLCalloc(sizeof(int),nCount);
        nResult = TransformEx( nCount, x, y, z, pabSuccess );
        CPLFree( pabSuccess );

        return nResult;
    }

    virtual int TransformEx( int nCount, double *x, double *y, double *z = NULL, int *pabSuccess = NULL ) {
        return GDALGenImgProjTransform( hSrcImageTransformer, TRUE, nCount, x, y, z, pabSuccess );
    }
};

namespace gip {
    using std::string;
	using std::vector;
	using std::cout;
	using std::endl;

	void test(const GeoImage& img) {
	    cout << img.Info() << endl;
	    cout << img[0].Info() << endl;

        GeoImageIO<float> img0(img);

	    img0[0] > 5;
	    cout << img0.Info() << endl;
	    cout << img.Info() << endl;
	    img0[0] = img0[0] > 5;
	    cout << img0.Info() << endl;
	    cout << img.Info() << endl;

	}

    //! Create mask based on NoData values for all bands
	GeoRaster CreateMask(const GeoImage& image, string filename) {
		typedef float T;

		GeoImageIO<T> img(image);
		CImg<T> imgchunk;

		//if (filename == "") filename = image.Basename() + "_mask";
		GeoImage mask(filename, img, GDT_Byte, 1);
		GeoRasterIO<unsigned char> mask0(mask[0]);
		CImg<unsigned char> imgout;
		mask.SetNoData(0);
		//unsigned int validpixels(0);
        for (unsigned int iChunk=1; iChunk<=img[0].NumChunks(); iChunk++) {
        	//imgchunk = imageIO[0].Read(*iChunk);
			//imgout = Pbands[0].NoDataMask(imgchunk);
			//for (unsigned int b=1;b<image.NumBands();b++) {
            //        imgout &= Pbands[b].NoDataMask( Pbands[b].Read(*iChunk) );
			//}
			//validpixels += imgout.sum();
			imgout = img.NoDataMask(iChunk);
			mask0.Write(imgout,iChunk);
		}
		//mask[0].SetValidSize(validpixels);
		return mask[0];
	}

	//! Calculate Radiance
	// TODO - combine with Copy/Process
	GeoImage Rad(const GeoImage& image, string filename) {
	    if (image[0].Units() != "radiance") {
	        throw std::runtime_error("image not in radiance units");
	    }
	    image.SetUnitsOut("radiance");
        GeoImageIO<float> img(image);
        GeoImageIO<float> imgout(GeoImage(filename, img, GDT_Int16));
        imgout.SetNoData(-32768); // TODO - set nodata option
        imgout.SetGain(0.1);
        imgout.SetUnits("radiance");
        CImg<float> cimg;
        Colors colors = img.GetColors();
        CImg<unsigned char> nodata;
        for (unsigned int b=0;b<img.NumBands();b++) {
            imgout.SetColor(colors[b+1], b+1);
            for (unsigned int iChunk=1; iChunk<=img[b].NumChunks(); iChunk++) {
                cimg = img[b].Read(iChunk);
                nodata = img[b].NoDataMask(iChunk);
                // only if nodata not same between input and output images
                cimg_forXY(cimg,x,y) { if (!nodata(x,y)) cimg(x,y) = imgout[b].NoDataValue(); }
                imgout[b].Write(cimg,iChunk);
            }
        }
        return imgout;
	}

    //! Calculate reflectance assuming read image is radiance
	GeoImage Ref(const GeoImage& image, string filename) {
        if (Options::Verbose() > 1)
            std::cout << "Reflectance(" << image.Basename() << ") -> " << filename << std::endl;
	    if ((image[0].Units() != "radiance") && (image[0].Units() != "reflectance")) {
	        throw std::runtime_error("image not in compatible units for reflectance");
	    }
        image.SetUnitsOut("reflectance");
        GeoImageIO<float> img(image);
        GeoImageIO<float> imgout(GeoImage(filename, img, GDT_Int16));
        imgout.SetNoData(-32768); // TODO - set nodata option
        imgout.SetUnits("reflectance");
        CImg<float> cimg;
        CImg<unsigned char> nodata;
        Colors colors = img.GetColors();
        for (unsigned int b=0;b<img.NumBands();b++) {
            if (img[b].Thermal()) imgout[b].SetGain(0.01); else imgout[b].SetGain(0.0001);
            imgout.SetColor(colors[b+1], b+1);
            for (unsigned int iChunk=1; iChunk<=img[b].NumChunks(); iChunk++) {
                cimg = img[b].Read(iChunk);
                nodata = img[b].NoDataMask(iChunk);
                cimg_forXY(cimg,x,y) { if (!nodata(x,y)) cimg(x,y) = imgout[b].NoDataValue(); }
                imgout[b].Write(cimg,iChunk);
            }
        }
        return imgout;
	}

	//! Calculate radar backscatter for all bands
	GeoImage SigmaNought(const GeoImage& image, string filename, float CF) {
	    GeoImageIO<float> img(image);
	    GeoImageIO<float> imgout(GeoImage(filename, img, GDT_Float32));
	    float nodataval = -32768;
	    imgout.SetNoData(nodataval);
	    img.SetUnitsOut("other");
	    imgout.SetUnits("other");
	    CImg<float> cimg;
	    CImg<unsigned char> nodata;
	    Colors colors = img.GetColors();
        for (unsigned int b=0;b<img.NumBands();b++) {
            imgout.SetColor(colors[b+1], b+1);
            for (unsigned int iChunk=1; iChunk<=img[b].NumChunks(); iChunk++) {
                cimg = img[b].Read(iChunk);
                cimg = cimg.pow(2).log10() * 10 + CF;
                nodata = img[b].NoDataMask(iChunk);
                cimg_forXY(cimg,x,y) { if (!nodata(x,y)) cimg(x,y) = nodataval; }
                imgout[b].Write(cimg,iChunk);
                //imgoutIO[b].Write(nodata,iChunk);
            }
        }
        return imgout;
	}

    //! Generate 3-band RGB image scaled to 1 byte for easy viewing
    GeoImage RGB(const GeoImage& image, string filename) {
        GeoImageIO<float> img(image);
        img.SetUnitsOut("reflectance");
        img.PruneToRGB();
        GeoImageIO<unsigned char> imgout(GeoImage(filename, img, GDT_Byte));
        imgout.SetNoData(0);
        imgout.SetUnits("other");
        CImg<float> stats, cimg;
        CImg<unsigned char> mask;
        for (unsigned int b=0;b<img.NumBands();b++) {
            stats = img[b].ComputeStats();
            float lo = std::max(stats(2) - 3*stats(3), stats(0)-1);
            float hi = std::min(stats(2) + 3*stats(3), stats(1));
            for (unsigned int iChunk=1; iChunk<=img[b].NumChunks(); iChunk++) {
                cimg = img[b].Read(iChunk);
                mask = img[b].NoDataMask(iChunk);
                ((cimg-=lo)*=(255.0/(hi-lo))).max(0.0).min(255.0);
                //cimg_printstats(cimg,"after stretch");
                cimg_forXY(cimg,x,y) { if (!mask(x,y)) cimg(x,y) = imgout[b].NoDataValue(); }
                imgout[b].Write(CImg<unsigned char>().assign(cimg.round()),iChunk);
            }
        }
        return imgout;
    }

	//! Merge images into one file and crop to vector
	GeoImage CookieCutter(vector<std::string> imgnames, string filename, string vectorname, float xres, float yres) {
	    // TODO - pass in vector of GeoRaster's instead
        if (Options::Verbose() > 2) {
            std::cout << filename << ": CookieCutter" << std::endl;
        }

        // Open input images
        vector<GeoImage> imgs;
        vector<std::string>::const_iterator iimgs;
        for (iimgs=imgnames.begin();iimgs!=imgnames.end();iimgs++) imgs.push_back(GeoImage(*iimgs));
        unsigned int bsz = imgs[0].NumBands();
        GDALDataType dtype = imgs[0].DataType();

	    // Create output file based on input vector
        OGRDataSource *poDS = OGRSFDriverRegistrar::Open(vectorname.c_str());
        OGRLayer *poLayer = poDS->GetLayer(0);
        OGREnvelope extent;
        poLayer->GetExtent(&extent, true);
        // Need to convert extent to resolution units
        int xsize = (int)(0.5 + (extent.MaxX - extent.MinX) / xres);
        int ysize = (int)(0.5 + (extent.MaxY - extent.MinY) / yres);
        GeoImage imgout(filename, xsize, ysize, bsz, dtype);
        imgout.CopyMeta(imgs[0]);
        imgout.CopyColorTable(imgs[0]);
        for (unsigned int b=0;b<bsz;b++) imgout[b].CopyMeta(imgs[0][b]);

        double affine[6];
        affine[0] = extent.MinX;
        affine[1] = xres;
        affine[2] = 0;
        affine[3] = extent.MaxY;
        affine[4] = 0;
        affine[5] = -yres;
        char* wkt = NULL;
        poLayer->GetSpatialRef()->exportToWkt(&wkt);
        imgout.GetGDALDataset()->SetProjection(wkt);
        imgout.GetGDALDataset()->SetGeoTransform(affine);
        // Compute union
        /*OGRPolygon* site; // = new OGRPolygon();
        OGRFeature *poFeature;
        poLayer->ResetReading();
        while( (poFeature = poLayer->GetNextFeature()) != NULL )
        {
            OGRGeometry *poGeometry;
            poGeometry = poFeature->GetGeometryRef();

            site = (OGRPolygon*)site->Union(poGeometry);
            OGRFeature::DestroyFeature( poFeature );
        }*/

        // Combine shape geoemtries into single geometry cutline
        OGRGeometry* site = OGRGeometryFactory::createGeometry( wkbMultiPolygon );
        OGRGeometry* poGeometry;
        OGRFeature *poFeature;
        poLayer->ResetReading();
        poFeature = poLayer->GetNextFeature();
        site = poFeature->GetGeometryRef();
        while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
            poGeometry = poFeature->GetGeometryRef();

            if( poGeometry == NULL ) fprintf( stderr, "ERROR: Cutline feature without a geometry.\n" );

            //OGRwkbGeometryType eType = wkbFlatten(poGeometry->getGeometryType());
            site = site->Union(poGeometry);
            /*if( eType == wkbPolygon )
                site->addGeometry(poGeometry);
            else if( eType == wkbMultiPolygon ) {
                for(int iGeom = 0; iGeom < OGR_G_GetGeometryCount( poGeometry ); iGeom++ ) {
                    site->addGeometry( poGeometry->getGeometryRef(iGeom)  );
                }
            }
            else fprintf( stderr, "ERROR: Cutline not of polygon type.\n" );*/

            OGRFeature::DestroyFeature( poFeature );
        }
        OGRDataSource::DestroyDataSource( poDS );

        // Cutline transform to pixel coordinates
        char **papszOptionsCutline = NULL;
        //papszOptionsCutline = CSLSetNameValue( papszOptionsCutline, "DST_SRS", wkt );
        //papszOptionsCutline = CSLSetNameValue( papszOptionsCutline, "SRC_SRS", wkt );
        //papszOptionsCutline = CSLSetNameValue( papszOptionsCutline, "INSERT_CENTER_LONG", "FALSE" );
        CutlineTransformer oTransformer;

        /* The cutline transformer will *invert* the hSrcImageTransformer */
        /* so it will convert from the cutline SRS to the source pixel/line */
        /* coordinates */
        oTransformer.hSrcImageTransformer = GDALCreateGenImgProjTransformer2( imgout.GetGDALDataset(), NULL, papszOptionsCutline );
        site->transform(&oTransformer);

        GDALDestroyGenImgProjTransformer( oTransformer.hSrcImageTransformer );
        CSLDestroy( papszOptionsCutline );

        // Warp options
        GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
        psWarpOptions->hDstDS = imgout.GetGDALDataset();
        psWarpOptions->nBandCount = bsz;
        psWarpOptions->panSrcBands = (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
        psWarpOptions->panDstBands = (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
        psWarpOptions->padfSrcNoDataReal = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
        psWarpOptions->padfSrcNoDataImag = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
        psWarpOptions->padfDstNoDataReal = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
        psWarpOptions->padfDstNoDataImag = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount );
        for (unsigned int b=0;b<bsz;b++) {
            psWarpOptions->panSrcBands[b] = b+1;
            psWarpOptions->panDstBands[b] = b+1;
            psWarpOptions->padfSrcNoDataReal[b] = imgs[0][b].NoDataValue();
            psWarpOptions->padfDstNoDataReal[b] = imgout[b].NoDataValue();
            psWarpOptions->padfSrcNoDataImag[b] = 0.0;
            psWarpOptions->padfDstNoDataImag[b] = 0.0;
        }
        if (Options::Verbose() > 2)
            psWarpOptions->pfnProgress = GDALTermProgress;
        else psWarpOptions->pfnProgress = GDALDummyProgress;
        char **papszOptions = NULL;
        //psWarpOptions->hCutline = site;
        //papszOptions = CSLSetNameValue(papszOptions,"SKIP_NOSOURCE","YES");
        //papszOptions = CSLSetNameValue(papszOptions,"INIT_DEST",NULL);
        //site->exportToWkt(&wkt);
        //papszOptions = CSLSetNameValue(papszOptions,"CUTLINE",wkt);
        psWarpOptions->papszWarpOptions = CSLDuplicate(papszOptions);
        //psWarpOptions->hCutline = site;

        GDALWarpOperation oOperation;
        // Perform warp for each input file
        vector<GeoImage>::iterator iimg;
        for (iimg=imgs.begin();iimg!=imgs.end();iimg++) {
            if (Options::Verbose() > 2) std::cout << "Warping file " << iimg->Basename() << std::endl;
            psWarpOptions->hSrcDS = iimg->GetGDALDataset();
            psWarpOptions->pTransformerArg =
                GDALCreateGenImgProjTransformer( iimg->GetGDALDataset(), iimg->GetGDALDataset()->GetProjectionRef(),
                                                imgout.GetGDALDataset(), imgout.GetGDALDataset()->GetProjectionRef(), TRUE, 0.0, 0 );
            psWarpOptions->pfnTransformer = GDALGenImgProjTransform;
            oOperation.Initialize( psWarpOptions );
            //std::cout << CPLGetLastErrorMsg() << std::endl;
            oOperation.ChunkAndWarpMulti( 0, 0, imgout.XSize(), imgout.YSize() );

            GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
        }
        GDALDestroyWarpOptions( psWarpOptions );

        return imgout;
	}

	void NDVI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, std::vector<std::string>({"NDVI"})); }
	void EVI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"EVI"}); }
    void LSWI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"LSWI"}); }
    void NDSI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"NDSI"}); }
    void BI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"BI"}); }
    void SATVI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"SATVI"}); }
    void NDTI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"NDTI"}); }
    void CRC(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"CRC"}); }
    void CRCm(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"CRCM"}); }
    void iSTI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"ISTI"}); }
    void STI(const GeoImage& ImageIn, std::string filename) { return Indices(ImageIn, filename, {"STI"}); }

	//! Create multi-band image of various indices calculated from input
	//GeoImage Indices(const GeoImage& ImageIn, string filename, bool ndvi, bool evi, bool lswi, bool ndsi, bool bi) {
	//void Indices(const GeoImage& ImageIn, string basename, std::initializer_list<std::string> list) {
    //    Indices(ImageIn, basename, std::vector<std::string>(list));
    //}

    void Indices(const GeoImage& ImageIn, string basename, std::vector<std::string> products) {
        ImageIn.SetUnitsOut("reflectance");
        GeoImageIO<float> imgin(ImageIn);
		float nodataout = -32768;

        std::map< string, GeoImageIO<float> > imagesout;
        vector<string>::const_iterator iprod;
        for (iprod=products.begin(); iprod!=products.end(); iprod++) {
            //imagesout[*iprod] = GeoImageIO<float>(GeoImage(basename + '_' + *iprod, imgin, GDT_Int16));
            imagesout[*iprod] = GeoImageIO<float>(GeoImage(basename, imgin, GDT_Int16, 1));
            imagesout[*iprod].SetNoData(nodataout);
            imagesout[*iprod].SetGain(0.0001);
            imagesout[*iprod].SetUnits("other");
            imagesout[*iprod][0].SetDescription(*iprod);
        }
        if (imagesout.size() == 0) throw std::runtime_error("No indices selected for calculation!");

        std::map< string, std::vector<string> > colors;
        colors["NDVI"] = {"NIR","RED"};
        colors["EVI"] = {"NIR","RED","BLUE"};
        colors["LSWI"] = {"NIR","SWIR1"};
        colors["NDSI"] = {"SWIR1","GREEN"};
        colors["BI"] = {"BLUE","NIR"};
        colors["SATVI"] = {"SWIR1","RED"};
        // Tillage indices
        colors["NDTI"] = {"SWIR2","SWIR1"};
        colors["CRC"] = {"SWIR1","SWIR2","BLUE"};
        colors["CRCM"] = {"SWIR1","SWIR2","GREEN"};
        colors["ISTI"] = {"SWIR1","SWIR2"};
        colors["STI"] = {"SWIR1","SWIR2"};

        // Figure out what colors are needed
        std::set< string > used_colors;
        std::set< string >::const_iterator isstr;
        std::vector< string >::const_iterator ivstr;
        for (iprod=products.begin(); iprod!=products.end(); iprod++) {
            for (ivstr=colors[*iprod].begin();ivstr!=colors[*iprod].end();ivstr++) {
                used_colors.insert(*ivstr);
            }
            if (Options::Verbose() > 2) std::cout << "Product " << *iprod << std::endl;
        }

		CImg<float> red, green, blue, nir, swir1, swir2, cimgout, cimgmask;

        // need to add overlap
        for (unsigned int iChunk=1; iChunk<=ImageIn[0].NumChunks(); iChunk++) {
            for (isstr=used_colors.begin();isstr!=used_colors.end();isstr++) {
                if (*isstr == "RED") red = imgin["RED"].Read(iChunk);
                else if (*isstr == "GREEN") green = imgin["GREEN"].Read(iChunk);
                else if (*isstr == "BLUE") blue = imgin["BLUE"].Read(iChunk);
                else if (*isstr == "NIR") nir = imgin["NIR"].Read(iChunk);
                else if (*isstr == "SWIR1") swir1 = imgin["SWIR1"].Read(iChunk);
                else if (*isstr == "SWIR2") swir2 = imgin["SWIR2"].Read(iChunk);
            }
            if (Options::Verbose() > 2) {
                std::cout << "Colors used: ";
                for (isstr=used_colors.begin();isstr!=used_colors.end();isstr++) std::cout << " " << *isstr;
                std::cout << std::endl;
            }

            for (iprod=products.begin(); iprod!=products.end(); iprod++) {
                std::cout << "product = " << *iprod << std::endl;
                //string p = iprod->toupper();
                if (*iprod == "NDVI") {
                    cimgout = (nir-red).div(nir+red);
                } else if (*iprod == "EVI") {
                    cimgout = 2.5*(nir-red).div(nir + 6*red - 7.5*blue + 1);
                } else if (*iprod == "LSWI") {
                    cimgout = (nir-swir1).div(nir+swir1);
                } else if (*iprod == "NDSI") {
                    cimgout = (green-swir1).div(green+swir1);
                } else if (*iprod == "BI") {
                    cimgout = 0.5*(blue+nir);
                } else if (*iprod == "SATVI") {
                    float L(0.5);
                    cimgout = (((1.0+L)*(swir1 - red)).div(swir1+red+L)) - (0.5*swir2);
                // Tillage indices
                } else if (*iprod == "NDTI") {
                    cimgout = (swir1-swir2).div(swir1+swir2);
                } else if (*iprod == "CRC") {
                    cimgout = (swir1-blue).div(swir2+blue);
                } else if (*iprod == "CRCM") {
                    cimgout = (swir1-green).div(swir2+green);
                } else if (*iprod == "ISTI") {
                    cimgout = swir2.div(swir1);
                } else if (*iprod == "STI") {
                    cimgout = swir1.div(swir2);
                }
                if (Options::Verbose() > 2) std::cout << "Getting mask" << std::endl;
                // TODO don't read mask again...create here
                cimgmask = imgin.NoDataMask(iChunk, colors[*iprod]);
                cimg_forXY(cimgout,x,y) if (!cimgmask(x,y)) cimgout(x,y) = nodataout;
                imagesout[*iprod].Write(cimgout,iChunk);
            }
        }
        //return imagesout;
	}

	//! Auto cloud mask - toaref input
	/*GeoImage AutoCloud(const GeoImage& image, string filename, int cheight, float minred, float maxtemp, float maxndvi, int morph) {
	    typedef float outtype;
        GeoImageIO<float> imgin(image);
        GeoImageIO<outtype> imgout(GeoImage(filename, image, GDT_Byte, 1));
        imgout.SetNoData(0);

        CImg<float> red, nir, temp, ndvi;
        CImg<outtype> mask;

        // need to add overlap
        for (int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {

            red = imgin["RED"].Ref(iChunk);
            temp = imgin["LWIR"].Ref(iChunk);
            nir = imgin["NIR"].Ref(iChunk);
            ndvi = (nir-red).div(nir+red);

            mask =
                temp.threshold(maxtemp,false,true)^=1 &
                ndvi.threshold(maxndvi,false,true)^=1 &
                red.get_threshold(minred);

            if (morph != 0) mask.dilate(morph);

            imgout[0].Write(mask,iChunk);
            //CImg<double> stats = img.get_stats();
            //cout << "stats " << endl;
            //for (int i=0;i<12;i++) cout << stats(i) << " ";
            //cout << endl;
        }
        return imgout;
	}*/

	//! ACCA (Automatic Cloud Cover Assessment) takes in TOA Reflectance and temperature
	GeoImage ACCA(const GeoImage& img, string filename) {
	    img.SetUnitsOut("reflectance");
        GeoImageIO<float> imgin(img);

        float th_red(0.08);
        float th_ndsi(0.7);
        float th_temp(27);
        float th_comp(225);
        float th_nirred(2.0);
        float th_nirgreen(2.0);
        float th_nirswir1(1.0);
        //float th_warm(210);

        GeoImageIO<unsigned char> imgout(GeoImage(filename, imgin, GDT_Byte, 4));
        imgout.SetNoData(0);
        imgout.SetUnits("other");
        imgout[0].SetDescription("finalmask");
        imgout[1].SetDescription("datamask");
        imgout[2].SetDescription("pass2");
        imgout[3].SetDescription("pass1");

        CImg<float> red, green, nir, swir1, temp, ndsi, b56comp;
        CImg<unsigned char> nonclouds, ambclouds, clouds, mask;
        float cloudsum(0), scenesize(0);

        if (Options::Verbose()) cout << img.Basename() << " - ACCA" << endl;
        for (unsigned int iChunk=1; iChunk<=imgin[0].NumChunks(); iChunk++) {
            red = imgin["RED"].Read(iChunk);
            green = imgin["GREEN"].Read(iChunk);
            nir = imgin["NIR"].Read(iChunk);
            swir1 = imgin["SWIR1"].Read(iChunk);
            temp = imgin["LWIR"].Read(iChunk);

            mask = imgin.NoDataMask(iChunk, {"RED","GREEN","NIR","SWIR1","LWIR"});

            ndsi = (green - swir1).div(green + swir1);
            b56comp = (1.0 - swir1).mul(temp + 273.15);

            // Pass one
            nonclouds = // 1's where they are non-clouds
                // Filter1
                (red.get_threshold(th_red)^=1) |=
                // Filter2
                ndsi.get_threshold(th_ndsi) |=
                // Filter3
                temp.get_threshold(th_temp);

            ambclouds =
                (nonclouds^1).mul(
                // Filter4
                b56comp.get_threshold(th_comp) |=
                // Filter5
                nir.get_div(red).threshold(th_nirred) |=
                // Filter6
                nir.get_div(green).threshold(th_nirgreen) |=
                // Filter7
                (nir.get_div(swir1).threshold(th_nirswir1)^=1) );

            clouds =
                (nonclouds + ambclouds)^=1;

                // Filter8 - warm/cold
                //b56comp.threshold(th_warm) + 1);

            //nonclouds.mul(mask);
            clouds.mul(mask);
            ambclouds.mul(mask);

            cloudsum += clouds.sum();
            scenesize += mask.sum();

            imgout[3].Write(clouds,iChunk);
            imgout[2].Write(ambclouds,iChunk);
            //imgout[0].Write(nonclouds,iChunk);
            if (Options::Verbose() > 3) std::cout << "Processed chunk " << iChunk << " of " << imgin[0].NumChunks() << std::endl;
        }
        // Cloud statistics
        float cloudcover = cloudsum / scenesize;
        CImg<float> tstats = imgin["LWIR"].AddMask(imgout[3]).ComputeStats();
        if (Options::Verbose() > 1) {
            cout.precision(4);
            cout << "   Cloud Cover = " << cloudcover*100 << "%" << endl;
            cimg_print(tstats, "Cloud stats(min,max,mean,sd,skew,count)");
        }

        // Pass 2 (thermal processing)
        bool addclouds(false);
        if ((cloudcover > 0.004) && (tstats(2) < 22.0)) {
            float th0 = imgin["LWIR"].Percentile(83.5);
            float th1 = imgin["LWIR"].Percentile(97.5);
            if (tstats[4] > 0) {
                float th2 = imgin["LWIR"].Percentile(98.75);
                float shift(0);
                shift = tstats[3] * ((tstats[4] > 1.0) ? 1.0 : tstats[4]);
                //cout << "Percentiles = " << th0 << ", " << th1 << ", " << th2 << ", " << shift << endl;
                if (th2-th1 < shift) shift = th2-th1;
                th0 += shift;
                th1 += shift;
            }
            imgin["LWIR"].ClearMasks();
            CImg<float> warm_stats = imgin["LWIR"].AddMask(imgout[2]).AddMask(imgin["LWIR"] < th1).AddMask(imgin["LWIR"] > th0).ComputeStats();
            if (Options::Verbose() > 1) cimg_print(warm_stats, "Warm Cloud stats(min,max,mean,sd,skew,count)");
            imgin["LWIR"].ClearMasks();
            if (((warm_stats(5)/scenesize) < 0.4) && (warm_stats(2) < 22)) {
                if (Options::Verbose() > 2) cout << "Accepting warm clouds" << endl;
                imgout[2].AddMask(imgin["LWIR"] < th1).AddMask(imgin["LWIR"] > th0);
                addclouds = true;
            } else {
                // Cold clouds
                CImg<float> cold_stats = imgin["LWIR"].AddMask(imgout[2]).AddMask(imgin["LWIR"] < th0).ComputeStats();
                if (Options::Verbose() > 1) cimg_print(cold_stats, "Cold Cloud stats(min,max,mean,sd,skew,count)");
                imgin["LWIR"].ClearMasks();
                if (((cold_stats(5)/scenesize) < 0.4) && (cold_stats(2) < 22)) {
                    if (Options::Verbose() > 2) cout << "Accepting cold clouds" << endl;
                    imgout[2].AddMask(imgin["LWIR"] < th0);
                    addclouds = true;
                } else
                    if (Options::Verbose() > 2) cout << "Rejecting all ambiguous clouds" << endl;
            }
        } else imgin["LWIR"].ClearMasks();

        int ksize(3);
        //CImg<int> selem(esize,esize,1,1,1);
        // 3x3 filter of 1's for majority filter
        CImg<int> filter(ksize,ksize,1,1, 1);
        int th_majority(((ksize*ksize)+1)/2);
        int padding((ksize+1)/2);

        //int esize = 5;
        //CImg<int> erode_elem(esize,esize,1,1,1);
        //CImg<int> dilate_elem(esize+dilate,esize+dilate,1,1,1);

        for (int b=0;b<4;b++) imgout[b].Chunk(ksize);
        for (int b=0;b<imgin.NumBands();b++) imgin[b].Chunk(ksize);

        for (unsigned int iChunk=1; iChunk<=imgout[0].NumChunks(); iChunk++) {
            if (Options::Verbose() > 3) std::cout << "Chunk " << iChunk << " of " << imgout[0].NumChunks() << std::endl;
            clouds = imgout[3].Read(iChunk);
            // should this be a |= ?
            if (addclouds) clouds += imgout[2].Read(iChunk);
            // Majority filter
            clouds|=(imgin.SaturationMask(iChunk)^=1);
            clouds = (clouds + clouds.get_convolve(filter).threshold(th_majority)).threshold(1);
            //clouds.erode(selem).dilate(selem);
            imgout[2].Write(clouds,iChunk);
            // Inverse and multiply by nodata mask to get good data mask
            imgout[1].Write((clouds^=1).mul(imgin.NoDataMask(iChunk)), iChunk);
            //imgout[0].Write(clouds|=imgin.SaturationMask(iChunk), iChunk);
            imgout[0].Write(clouds, iChunk);
            //imgout[0].Write(clouds^=1, iChunk);

            // Erode, then dilate twice
            //mask.erode(erode_elem).dilate(dilate_elem);
            //mask.dilate(dilate_elem);
        }

        return imgout;
	}

    //! Fmask cloud mask
    /*GeoImage Fmask(const GeoImage& image, string filename, int tolerance, int dilate) {
        image.SetUnitsOut("reflectance");
        GeoImageIO<float> imgin(image);

        // Output masks
        GeoImageIO<unsigned char> imgout(GeoImage(filename,image,GDT_Byte,2));
        imgout.SetNoData(0);

        // Output probabilties (for debugging/analysis)
        GeoImageIO<float> probout(GeoImage(filename + "_prob", image, GDT_Float32, 1));
        float nodata = -32768;
        probout.SetNoData(nodata);

        CImg<unsigned char> mask, wmask, lmask, nodatamask, redsatmask, greensatmask;
        CImg<float> red, nir, green, blue, swir1, swir2, BT, ndvi, ndsi, white, vprob;
        float _ndvi, _ndsi;
        //int cloudpixels(0);
        CImg<double> wstats(image.Size()), lstats(image.Size());
        int wloc(0), lloc(0);

        for (unsigned int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
            blue = imgin["BLUE"].Read(iChunk);
            red = imgin["RED"].Read(iChunk);
            green = imgin["GREEN"].Read(iChunk);
            nir = imgin["NIR"].Read(iChunk);
            swir1 = imgin["SWIR1"].Read(iChunk);
            swir2 = imgin["SWIR2"].Read(iChunk);
            BT = imgin["LWIR"].Read(iChunk);
            nodatamask = imgin.NoDataMask(iChunk);
            // floodfill....seems bad way
            //shadowmask = nir.draw_fill(nir.width()/2,nir.height()/2,)
            ndvi = (nir-red).div(nir+red);
            ndsi = (green-swir1).div(green+swir1);
            redsatmask = imgin["RED"].SaturationMask(iChunk);
            greensatmask = imgin["GREEN"].SaturationMask(iChunk);
            white = imgin.Whiteness(iChunk);
            vprob = red;

            // Calculate "variability probability"
            cimg_forXY(vprob,x,y) {
                _ndvi = (redsatmask(x,y) && nir(x,y) > red(x,y)) ? 0 : abs(ndvi(x,y));
                _ndsi = (greensatmask(x,y) && swir1(x,y) > green(x,y)) ? 0 : abs(ndsi(x,y));
                vprob(x,y) = 1 - std::max(white(x,y), std::max(_ndsi, _ndvi));
            }
            probout[0].Write(vprob, iChunk);

            // Potential cloud layer
            mask =
                swir2.threshold(0.03)
                & BT.get_threshold(27,false,true)^=1
                // NDVI
                & ndvi.threshold(0.8,false,true)^=1
                // NDSI
                & ndsi.threshold(0.8,false,true)^=1
                // HazeMask
                & (blue - 0.5*red - 0.08).threshold(0.0)
                & imgin.Whiteness(iChunk).threshold(0.7,false,true)^=1
                & nir.div(swir1).threshold(0.75);

            //cloudpixels += mask.sum();

            // Water and land masks
            CImg<bool> wmask(red.width(),red.height(),1,1,false);
            cimg_forXY(wmask,x,y) {
                if ( ((ndvi(x,y) < 0.01) && (nir(x,y) < 0.11)) || ((ndvi(x,y) < 0.1) && (nir(x,y) < 0.05)) ) wmask(x,y) = true;
            }
            imgout[0].Write(mask, iChunk);
            imgout[1].Write(wmask, iChunk);

            //lmask = wmask^=1 & mask^=1;
            lmask = (wmask^1) & (mask^=1) & nodatamask;
            wmask = wmask & swir2^=1 & nodatamask;

            cimg_forXY(BT,x,y) if (wmask(x,y)) wstats[wloc++] = BT(x,y);
            cimg_forXY(BT,x,y) if (lmask(x,y)) lstats[lloc++] = BT(x,y);

            //lBT = BT;
            //wBT = BT;
            //cimg_forXY(BT,x,y) if (!wmask(x,y)) wBT(x,y) = nodata;
            //cimg_forXY(BT,x,y) if (!mask(x,y)) lBT(x,y) = nodata;
            //probout[0].WriteChunk(wBT,*iChunk);
            //probout[1].WriteChunk(lBT,*iChunk);
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

        if (Options::Verbose() > 1) {
            cout << "Running fmask in " << image[0].NumChunks() << " chunks with tolerance of " << tolerance << endl;
            cout << "Temperature stats:" << endl;
            cout << "  Water: " << wmean << " s = " << wstddev << endl;
            cout << "  Water: " << lmean << " s = " << lstddev << endl;
            cout << "  Twater (82.5%) = " << Twater << endl;
            cout << "  Tlo (17.5%) = " << Tlo << endl;
            cout << "  Thi (82.5%) = " << Thi << endl;
        }

        // 2nd pass cloud probabilities over land
        CImg<float> wprob, lprob;
        for (unsigned int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
            // Cloud over Water prob
            BT = imgin["LWIR"].Read(iChunk);
            nodatamask = imgin.NoDataMask(iChunk);

            vprob = probout[0].Read(iChunk);

            // temp probability x variability probability
            lprob = ((Thi + 4-BT)/=(Thi+4-(Tlo-4))).mul( vprob );
            //1 - imgin.NDVI(*iChunk).abs().max(imgin.NDSI(*iChunk).abs()).max(imgin.Whiteness(*iChunk).abs()) );

            cimg_forXY(nodatamask,x,y) if (!nodatamask(x,y)) lprob(x,y) = nodata;

            // Cloud probability over water
            probout[0].Write( lprob, iChunk);
        }

        // Apply thresholds to get final max
        float tol = (tolerance-3)*0.1;
        float wthresh = 0.5 + tol;
        CImg<double> stats = probout[0].ComputeStats();
        float lthresh( zhi*stats(3) + stats(2) + 0.2 + tol );

        if (Options::Verbose() > 1) {
            std::cout << "Cloud probability thresholds:" << std::endl;
            std::cout << "  Over Water = " << wthresh << std::endl;
            std::cout << "  Over Land = " << lthresh << std::endl;
        }

        // 3x3 filter of 1's for majority filter
        CImg<int> filter(3,3,1,1, 1);
        int esize = 5;
        CImg<int> erode_elem(esize,esize,1,1,1);
        CImg<int> dilate_elem(esize+dilate,esize+dilate,1,1,1);

        for (unsigned int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
            mask = imgout[0].Read(iChunk);
            wmask = imgout[1].Read(iChunk);

            nodatamask = imgin.NoDataMask(iChunk);
            BT = imgin["LWIR"].Read(iChunk);
            swir1 = imgin["SWIR1"].Read(iChunk);

            // temp probability x brightness probability
            wprob = ((Twater - BT)/=4.0).mul( swir1.min(0.11)/=0.11 );
            lprob = probout[0].Read(iChunk);

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

            cimg_forXY(nodatamask,x,y) if (!nodatamask(x,y)) mask(x,y) = 0;
            imgout[0].Write(mask, iChunk);
        }

        // Add shadow mask: bitmask

        // NoData regions
        return imgout;
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

        CImg<float> cimgin, cimgout;
        for (unsigned int iChunk=1; iChunk<=image[bandnum-1].NumChunks(); iChunk++) {
            cimgin = imagein[bandnum-1].Read(iChunk);
            cimgout = (cimgin - min)/(max-min);
            cimgout.min(1.0).max(0.0);
            cimg_forXY(cimgin,x,y) if (cimgin(x,y) == nodatain) cimgout(x,y) = nodataout;
            imageout[0].Write(cimgout, iChunk);
            cimg_for(cimgout,ptr,float) if (*ptr != nodataout) *ptr = 1.0 - *ptr;
            imageout[1].Write(cimgout, iChunk);
        }
        return imageout;
    }

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

        int NumPixelChange, iteration=0;
        do {
            NumPixelChange = 0;
            for (i=0; i<classes; i++) NumSamples(i) = 0;
            if (Options::Verbose()) cout << "  Iteration " << iteration+1 << std::flush;

            for (unsigned int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
                C_img = img.Read(iChunk);
                C_mask = img.NoDataMask(iChunk);
                C_imgout = imgout[0].Read(iChunk);

                CImg<double> stats;
                cimg_forXY(C_img,x,y) { // Loop through image
                    // Calculate distance between this pixel and all classes
                    if (C_mask(x,y)) {
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
                imgout[0].Write(C_imgout,iChunk);
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
/*
    //! Rice detection algorithm
    GeoImage RiceDetect(const GeoImage& img, string filename, vector<int> days, int maxcrops, float th0, float th1, int dth0, int dth1) {

        GeoImageIO<float> imgio(img);
        int numbands = 2 * maxcrops + 3;
        GeoImage imgout(filename, img, GDT_Byte, numband);
        GeoImage<unsigned char> imgoutio(imgout);

        if (Options::Verbose() > 1) {
            vector<int>::const_iterator v;
            std::cout << "Days ";
            for (v=days.begin(); v!=days.end(); v++) std::cout << *v << " ";
            std::cout << endl;
        }

        int chunknum(0);
        int peaknum(0);
        for (int iChunk=1; iChunk<=image[0].NumChunks(); iChunk++) {
            CImgList<T> CImage( Image.ReadAsList(*iChunk) );
            CImgList<T> CImageOut( ImageOut.ReadAsList(*iChunk) );

            CImg<int> cimgout1, cimgout2;

            // get chunk size?
            p1 = iChunk->min_corner();
            p2 = iChunk->max_corner();
            width = p2.x()-p1.x()+1;
            height = p2.y()-p1.y()+1;
            // Reset running DOY to all zero
            CImg<int> DOY( width, height, 1, 1, 0 );

            CImg<int> Clow, Chigh, matched;

            // Seed with first band
            CImg<float> cimg = imgio[0].Read(iChunk);
            CImg<int> ClowEver( cimg.get_threshold(th0)^1 );
            CImg<int> ChighEver( cimg.get_threshold(th1,false,true) );
            cimgout1 = ClowEver;
            cimgout2 = ChighEver;
            for (unsigned int b=1;b<Image.NumBands();b++) {
                cimg = imgio[b].Read(iChunk);
                // Get low and high thresholds for this band
                Clow = cimg.get_threshold(th0);
                Chigh = cimg.get_threshold(th1,false,true);

                // Temporal processing
                // Where <= low threshold, set DOY to 0
                DOY.mul(Clow);
                // Where > low threshold, add to previous DOY
                DOY += Clow * doy[b];

                // Update flags if it was Ever low/high
                ClowEver |= (Clow^=1); // Clow now represents <= th0
                ChighEver |= Chigh;

                // Add to number of troughs/peaks
                cimgout1 += Clow;
                cimgout2 += Chigh;

                // More temporal processing
                // If low threshold was never met, change DOY to zero
                DOY.mul(ClowEver^1);
                Chigh = Chigh & DOY.get_threshold(dth0,true) & (DOY.get_threshold(dth1)^1);
                matched += Chigh;
                // Reset if high date has passed
                DOY.mul(DOY.get_threshold(dth1));
                // Loop through
                if (maxcrops != 0) {
                    cimg_forXY(CImageOut[0],x,y) {
                        if (Clow(x,y)) {
                            peaknum = CImageOut[1](x,y);
                            if (peaknum > 0 && peaknum <= maxcrops) {
                                CImageOut[(peaknum-1)*2+3](x,y) = doy[b];
                            }
                        }
                        if (Chigh(x,y)) {
                            peaknum = CImageOut[2](x,y);
                            if (peaknum > 0 && peaknum <= maxcrops) {
                                CImageOut[(peaknum-1)*2+4](x,y) = doy[b];
                                // Day of Year
                                CImageOut[(peaknum-1)*2+3](x,y) = DOY(x,y); // * Chigh(x,y);
                                // Yield (value)
                                CImageOut[(peaknum-1)*2+4](x,y) = CImage[b](x,y);
                            }
                        }

                    }
                }
            }
            CImageOut[0] = ClowEver & ChighEver;
            CImg<T> tmp = CImageOut[0].get_threshold(0);
            for (unsigned int c=0;c<CImageOut.size();c++) {
                CImageOut[c].mul(tmp);
            }

            // Write out images
            ImageOut.Write(CImageOut.get_append('c'), *iChunk, true);
        }

    }
*/

    // Perform band math (hard coded subtraction)
    /*GeoImage BandMath(const GeoImage& image, string filename, int band1, int band2) {
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
            imageout[0].WriteChunk(cimgout,*iChunk);
        }
        return imageout;
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

} // namespace gip

