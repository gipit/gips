#include <gip/GeoRaster.h>
#include <gip/GeoImage.h>

#include <gip/GeoRasterIO.h>

using namespace std;

namespace gip {
	// Copy constructor
	GeoRaster::GeoRaster(const GeoRaster& image)
		: GeoData(image), _GDALRasterBand(image._GDALRasterBand), _Masks(image._Masks), _NoData(image._NoData), _ValidSize(image._ValidSize),
            _minDC(image._minDC), _maxDC(image._maxDC), _K1(image._K1), _K2(image._K2), _Esun(image._Esun),
            _Atmosphere(image._Atmosphere), _Functions(image._Functions) {
		//std::cout << Basename() << ": GeoRaster copy (" << this << ")" << std::endl;
	}

	// Copy and add processing function
	GeoRaster::GeoRaster(const GeoRaster& image, GeoFunction func)
		: GeoData(image), _GDALRasterBand(image._GDALRasterBand), _Masks(image._Masks), _NoData(image._NoData), _ValidSize(image._ValidSize),
            _minDC(image._minDC), _maxDC(image._maxDC), _K1(image._K1), _K2(image._K2), _Esun(image._Esun),
            _Atmosphere(image._Atmosphere), _Functions(image._Functions) {
		_Functions.push_back(func);
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
		_ValidSize = image._ValidSize;
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

    //! Copy passed raster band into this raster band
    GeoRaster& GeoRaster::Copy(const GeoRaster& img, UNITS units) {
        switch (DataType()) {
            case GDT_Byte: GeoRasterIO<unsigned char>(*this).Copy(img, units);
                break;
            case GDT_UInt16: GeoRasterIO<unsigned short>(*this).Copy(img, units);
                break;
            case GDT_Int16: GeoRasterIO<short>(*this).Copy(img, units);
                break;
            case GDT_UInt32: GeoRasterIO<unsigned int>(*this).Copy(img, units);
                break;
            case GDT_Int32: GeoRasterIO<int>(*this).Copy(img, units);
                break;
            case GDT_Float32: GeoRasterIO<float>(*this).Copy(img, units);
                break;
            case GDT_Float64: GeoRasterIO<double>(*this).Copy(img, units);
                break;
            default: GeoRasterIO<unsigned char>(*this).Copy(img, units);
                break;
        }
        return *this;
    }

    //! Compute stats
    cimg_library::CImg<float> GeoRaster::ComputeStats(UNITS units) const {
        switch (DataType()) {
            case GDT_Byte: return GeoRasterIO<unsigned char>(*this).ComputeStats(units);
            case GDT_UInt16: return GeoRasterIO<unsigned short>(*this).ComputeStats(units);
            case GDT_Int16: return GeoRasterIO<short>(*this).ComputeStats(units);
            case GDT_UInt32: return GeoRasterIO<unsigned int>(*this).ComputeStats(units);
            case GDT_Int32: return GeoRasterIO<int>(*this).ComputeStats(units);
            case GDT_Float32: return GeoRasterIO<float>(*this).ComputeStats(units);
            case GDT_Float64: return GeoRasterIO<double>(*this).ComputeStats(units);
            default: return GeoRasterIO<unsigned char>(*this).ComputeStats(units);
            // TODO - remove default. This should throw exception
        }
    }

	string GeoRaster::Info(bool stats) const {
		std::stringstream info;
		//info << _GeoImage->Basename() << " - b" << _GDALRasterBand->GetBand() << ":" << endl;
		info << XSize() << " x " << YSize() << " " << DataType() << ": " << Description();
		info << " (GeoData: " << _GDALDataset.use_count() << " " << _GDALDataset << ")";
		info << " RasterBand &" << _GDALRasterBand << endl;
        info << "\t\tGain = " << Gain() << ", Offset = " << Offset(); //<< ", Units = " << Units();
        if (_NoData)
			info << ", NoData = " << NoDataValue() << endl;
        else info << endl;
        if (stats) {
        	info << "\t\tMin = " << Min() << ", Max = " << Max() << ", Mean = " << Mean() << " =/- " << StdDev() << endl;
        }
        if (!_Functions.empty()) info << "\t\tFunctions:" << endl;
        for (unsigned int i=0;i<_Functions.size();i++) {
        	info << "\t\t\t" << _Functions[i].Function() << " " << _Functions[i].Operand() << endl;
        }
		//_GeoImage->GetGDALDataset()->Reference(); int ref = _GeoImage->GetGDALDataset()->Dereference();
		//info << "  GDALDataset: " << _GDALDataset.use_count() << " (&" << _GDALDataset << ")" << endl;
        return info.str();
	}

} // namespace gip
