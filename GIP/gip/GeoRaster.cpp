#include <gip/GeoRaster.h>
#include <gip/GeoImage.h>

#include <gip/GeoRasterIO.h>

using namespace std;

namespace gip {

	// Copy constructor
	GeoRaster::GeoRaster(const GeoRaster& image, GeoFunction func)
		: GeoData(image), _GDALRasterBand(image._GDALRasterBand), _Masks(image._Masks), _NoData(image._NoData), //_ValidSize(image._ValidSize),
            _minDC(image._minDC), _maxDC(image._maxDC), _K1(image._K1), _K2(image._K2), _Esun(image._Esun),
            _Atmosphere(image._Atmosphere), _Functions(image._Functions) {
        if (func.Function() != "") AddFunction(func);
		//std::cout << Basename() << ": GeoRaster copy (" << this << ")" << std::endl;
	}

	// Assignment
	GeoRaster& GeoRaster::operator=(const GeoRaster& image) {
		// Check for self assignment
		if (this == &image) return *this;
		//_GeoData = image._GeoData;
		GeoData::operator=(image);
		_GDALRasterBand = image._GDALRasterBand;
		_Masks = image._Masks;
		_NoData = image._NoData;
		//_ValidSize = image._ValidSize;
		_minDC = image._minDC;
		_maxDC = image._maxDC;
		_K1 = image._K1;
		_K2 = image._K2;
		_Esun = image._Esun;
		_Functions = image._Functions;
		_Atmosphere = image._Atmosphere;
		//cout << _GeoImage->Basename() << ": " << ref << " references (GeoRaster Assignment)" << endl;
		return *this;
	}

    //! Process passed raster band into this raster band
    GeoRaster& GeoRaster::Process(const GeoRaster& img, bool RAW) {
        switch (DataType()) {
            case GDT_Byte: return GeoRasterIO<unsigned char>(*this).Process(img, RAW);
            case GDT_UInt16: return GeoRasterIO<unsigned short>(*this).Process(img, RAW);
            case GDT_Int16: return GeoRasterIO<short>(*this).Process(img, RAW);
            case GDT_UInt32: return GeoRasterIO<unsigned int>(*this).Process(img, RAW);
            case GDT_Int32: return GeoRasterIO<int>(*this).Process(img, RAW);
            case GDT_Float32: return GeoRasterIO<float>(*this).Process(img, RAW);
            case GDT_Float64: return GeoRasterIO<double>(*this).Process(img, RAW);
            default: return GeoRasterIO<unsigned char>(*this).Process(img, RAW);
            // TODO - remove default. This should throw exception
        }
    }

	string GeoRaster::Info(bool showstats) const {
		std::stringstream info;
		//info << _GeoImage->Basename() << " - b" << _GDALRasterBand->GetBand() << ":" << endl;
		info << XSize() << " x " << YSize() << " " << DataType() << ": " << Description();
		info << " (GeoData: " << _GDALDataset.use_count() << " " << _GDALDataset << ")";
		info << " RasterBand &" << _GDALRasterBand << endl;
        info << "\t\tGain = " << Gain() << ", Offset = " << Offset(); //<< ", Units = " << Units();
        if (_NoData)
			info << ", NoData = " << NoDataValue() << endl;
        else info << endl;
        if (showstats) {
            cimg_library::CImg<float> stats = this->ComputeStats();
        	info << "\t\tMin = " << stats(0) << ", Max = " << stats(1) << ", Mean = " << stats(2) << " =/- " << stats(3) << endl;
        }
        if (!_Functions.empty()) info << "\t\tFunctions:" << endl;
        for (unsigned int i=0;i<_Functions.size();i++) {
        	info << "\t\t\t" << _Functions[i].Function() << " " << _Functions[i].Operand() << endl;
        }
		//_GeoImage->GetGDALDataset()->Reference(); int ref = _GeoImage->GetGDALDataset()->Dereference();
		//info << "  GDALDataset: " << _GDALDataset.use_count() << " (&" << _GDALDataset << ")" << endl;
        return info.str();
	}

    //! Compute stats
    cimg_library::CImg<float> GeoRaster::ComputeStats(bool RAW) const {
        CImg<float> stats;
        /*switch (DataType()) {
            case GDT_Byte: stats = GeoRasterIO<unsigned char>(*this).ComputeStats(RAW);
            case GDT_UInt16: stats = GeoRasterIO<unsigned short>(*this).ComputeStats(RAW);
            case GDT_Int16: stats = GeoRasterIO<short>(*this).ComputeStats(RAW);
            case GDT_UInt32: stats = GeoRasterIO<unsigned int>(*this).ComputeStats(RAW);
            case GDT_Int32: stats = GeoRasterIO<int>(*this).ComputeStats(RAW);
            case GDT_Float32: stats = GeoRasterIO<float>(*this).ComputeStats(RAW);
            case GDT_Float64: stats = GeoRasterIO<double>(*this).ComputeStats(RAW);
            default: stats = GeoRasterIO<unsigned char>(*this).ComputeStats(RAW);
            // TODO - remove default. This should throw exception
        }*/
        stats = GeoRasterIO<double>(*this).ComputeStats(RAW);
        if (Options::Verbose() > 1) {
            cout << Basename() << " stats: " << stats(0);
            for (int i=1;i<4;i++) cout << ", " << stats(i);
        }
        return stats;
    }

    //! Compute histogram
    cimg_library::CImg<long> GeoRaster::Histogram(int bins) const {
        CImg<double> cimg;
        CImg<float> stats = ComputeStats();
        CImg<long> hist(bins,1,1,1,0);
        float nodata = NoData();
        GeoRasterIO<double> img(GeoRasterIO<double>(*this));
        //CImg<long> hist = GeoRasterIO<double>(*this).Read(chunknum).get_histogram(bins); //, stats(0), stats(1));
        for (int iChunk=1; iChunk<=NumChunks(); iChunk++) {
            cimg = img.Read(iChunk);
            cimg_for(cimg,ptr,double) {
                if (*ptr != nodata) hist[(unsigned int)((*ptr-stats(0)*bins/(stats(1)-stats(0))))]++;
            }
        }
        if (Options::Verbose() > 3) hist.display_graph(0,3);
        return hist;
    }

} // namespace gip
