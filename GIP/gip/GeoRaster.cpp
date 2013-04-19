#include <gip/GeoRaster.h>
#include <gip/GeoImage.h>

using namespace std;

namespace gip {
	// Copy constructor
	GeoRaster::GeoRaster(const GeoRaster& image)
		: GeoData(image), _GDALRasterBand(image._GDALRasterBand), _NoData(image._NoData), _ValidSize(image._ValidSize),
            _Sensor(image._Sensor), _Atmosphere(image._Atmosphere), _Functions(image._Functions) {
		//std::cout << Basename() << ": GeoRaster copy (" << this << ")" << std::endl;
	}

	// Copy and add processing function
	GeoRaster::GeoRaster(const GeoRaster& image, GeoFunction func)
		: GeoData(image), _GDALRasterBand(image._GDALRasterBand), _NoData(image._NoData), _ValidSize(image._ValidSize),
            _Sensor(image._Sensor), _Atmosphere(image._Atmosphere), _Functions(image._Functions) {
		_Functions.push_back(func);
	}

	// Assignment
	GeoRaster& GeoRaster::operator=(const GeoRaster& image) {
		// Check for self assignment
		if (this == &image) return *this;
		//_GeoData = image._GeoData;
		GeoData::operator=(image);
		_GDALRasterBand = image._GDALRasterBand;
		_NoData = image._NoData;
		_ValidSize = image._ValidSize;
		_Functions = image._Functions;
		_Sensor = image._Sensor;
		_Atmosphere = image._Atmosphere;
		//cout << _GeoImage->Basename() << ": " << ref << " references (GeoRaster Assignment)" << endl;
		return *this;
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
