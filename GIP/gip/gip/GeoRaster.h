#ifndef GIP_GEORASTER_H
#define GIP_GEORASTER_H

#include <boost/filesystem.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <gip/gip_CImg.h>
#include <gip/GeoData.h>
#include <gip/Atmosphere.h>

#include <iostream>
#include <iomanip>

// for tolowercase
#include <boost/algorithm/string.hpp>

namespace gip {
	//typedef bgeo::model::d2::point_xy<float> point;

    //enum UNITS {RAW, RADIANCE, REFLECTIVITY};
    //enum READMODE {NORMAL, RAW, ATMOSPHERE};
	//template<class T> class GeoRasterIO;

	class GeoFunction {
	public:
        GeoFunction()
            : _Function("") {}
		GeoFunction(std::string f, double v)
			: _Function(f), _dOperand(v) {
		}
		std::string Function() const { return _Function; }
		//string Operand() const { return _Operand; }
		double Operand() const { return _dOperand; }
	private:
		std::string _Function;
		//string _Operand;
		double _dOperand;
	};

	//! Extended GDALRasterBand class
	/*!
		The GeoRaster class wraps the GDALRasterBand class
	*/
	class GeoRaster : public GeoData {
		friend class GeoImage;
	public:
		//! \name Constructors/Destructors
		//! Constructor for new band
		GeoRaster(const GeoData& geodata, int bandnum=1)
            : GeoData(geodata), _NoData(false), //_ValidSize(0),
            _minDC(1), _maxDC(255), _K1(0), _K2(0), _Esun(0),
            _Atmosphere() {
			LoadBand(bandnum);
		}
		//! Copy constructor
		GeoRaster(const GeoRaster& image, GeoFunction func=GeoFunction());
		// Copy and add new function
		//GeoRaster(const GeoRaster& image, GeoFunction func);
		//! Assignment Operator
		GeoRaster& operator=(const GeoRaster& image);
		//! Destructor
		~GeoRaster() {}

		//! \name File Information
		//! Band X Size
		unsigned int XSize() const { return _GDALRasterBand->GetXSize(); }
		//! Band Y Size
		unsigned int YSize() const { return _GDALRasterBand->GetYSize(); }
		//! Get GDALDatatype
		GDALDataType DataType() const { return _GDALRasterBand->GetRasterDataType(); }
		//! Output file info
		std::string Info(bool showstats=false) const;

		// Total Valid Pixels
		//unsigned int ValidSize() const { return _ValidSize; }
		// Get datatype
		//std::string DataTypeStr() const { return GDALGetDataTypeName(DataType()); }

		//! Get GDALRasterBand object - use cautiously
		GDALRasterBand* GetGDALRasterBand() const { return _GDALRasterBand; }

		//! \name Metadata functions
		//! Get Band Description
		std::string Description() const {
			return _GDALRasterBand->GetDescription();
		}
		//! Set Band Description
		void SetDescription(std::string desc) {
			_GDALRasterBand->SetDescription(desc.c_str());
			// Also set description in dataset metadata since band desc doesn't work at least in GTiff
			this->SetMeta("Band "+to_string(_GDALRasterBand->GetBand()), desc);
		}
		//! Set Color Interp
		void SetColor(std::string col) {
		    // Is this used in other GDAL aware programs?
		    GDALColorInterp gdalcol;
		    if (col == "Red")
                gdalcol = GCI_RedBand;
            else if (col == "Green")
                gdalcol = GCI_GreenBand;
            else if (col == "Blue")
                gdalcol = GCI_BlueBand;
            else gdalcol = GCI_GrayIndex;
			_GDALRasterBand->SetColorInterpretation(gdalcol);
		}
        //! Copy category names from another band
		void CopyCategoryNames(const GeoRaster& raster) {
            _GDALRasterBand->SetCategoryNames(raster.GetGDALRasterBand()->GetCategoryNames());
		}
		//! Get GDAL Unit type
		std::string Units() const {
		    std::string units( _GDALRasterBand->GetUnitType() );
		    boost::algorithm::to_lower(units);
            return units;
        }
		//! Get gain
		float Gain() const { return _GDALRasterBand->GetScale(); }
		//! Get offset
		float Offset() const { return _GDALRasterBand->GetOffset(); }
		//! Set Unit type
		GeoRaster& SetUnits(std::string units) { _GDALRasterBand->SetUnitType(units.c_str()); return *this; }
		//! Set gain
		GeoRaster& SetGain(float gain) { _GDALRasterBand->SetScale(gain); return *this; }
		//! Set offset
		GeoRaster& SetOffset(float offset) { _GDALRasterBand->SetOffset(offset); return *this; }
		//! Flag indicating if NoData value is used or not
		bool NoData() const { return _NoData; }
		//! Get NoDataValue
		double NoDataValue() const {
		    //std::cout << "NoDataValue" << std::endl;
			//if (_NoData) return _GDALRasterBand->GetNoDataValue(); else return 0;
			return _GDALRasterBand->GetNoDataValue();
		}
		//! Set No Data value
		GeoRaster& SetNoData(double val) {
		    //std::cout << "SetNoData " << val << std::endl;
			_GDALRasterBand->SetNoDataValue(val);
			_NoData = true;
			return *this;
		}
		//! Clear NoData
		void ClearNoData() {
			_GDALRasterBand->SetNoDataValue( MaxValue() + 1 );
			_NoData = false;
		}
		//! Return maximum value based on datatype
		double MaxValue() const {
		    // TODO - base this on platform, not hard-coded
			switch ( DataType() ) {
				case GDT_Byte: return 255;
				case GDT_UInt16: return 65535;
				case GDT_Int16: return 32767;
				case GDT_UInt32: return 4294967295;
				case GDT_Int32: return 2147183647;
				case GDT_Float32: return 3.4E38;
				default: return 1.79E308;
			}
		}
		//! Return minimum value based on datatype (TODO - get from limits?)
		double MinValue() const {
            switch (DataType()) {
				case GDT_Byte: return 0;
				case GDT_UInt16: return 0;
				case GDT_Int16: return -32768;
				case GDT_UInt32: return 0;
				case GDT_Int32: return -2147183648;
				case GDT_Float32: return -3.4E38;
				default: return -1.79E308;
            }
		}
		// Set NoDataValue to a default based on datatype of band
		/*GeoRaster& SetNoData() {
		    //std::cout << "SetNoData" << std::endl;
			// only set if not already set
			if (_NoData) return;
			double val = -32766;
			GDALDataType dt = DataType();
			if (dt == GDT_Byte) val = 255;
			else if (dt == GDT_UInt16 || dt == GDT_UInt32) val = -val;
			SetNoData(val);
			return *this;
		}*/

        //! \name Calibration and atmospheric functions
		//! Is this a thermal sensor band?
		bool Thermal() const { if (_K1*_K2 == 0) return false; else return true; }
		//! Set thermal band, calling wth no arguments will clear thermal band status
		void SetThermal(float k1=0, float k2=0) { _K1=k1; _K2=k2; }
        //! Sets dyanmic range of sensor (min to max digital counts)
        void SetDynamicRange(int min, int max) {
            _minDC = min;
            _maxDC = max;
        }
        //! Set exo-atmospheric solar irradiance
        void SetEsun(float E) { _Esun = E; }
        // TODO - Does there need to be an atmospheric class???
		//! Is there an atmospheric correction supplied?
		bool Atmosphere() const { return _Atmosphere.Valid(); }
		//! Set atmospheric correction parameters
		GeoRaster& SetAtmosphere(gip::Atmosphere atm) { _Atmosphere = atm; return *this; }
        //! Clear atmospheric correction
        GeoRaster& ClearAtmosphere() { _Atmosphere = gip::Atmosphere(); return *this; }

        //! \name Processing functions
        //! Copy raster band into this raster
        GeoRaster& Process(const GeoRaster&, bool RAW=false);

        //! Adds a mask band (1 for valid), applied on read
        GeoRaster& AddMask(const GeoRaster& band) { _Masks.push_back(band); return *this; }
        //! Remove all masks from band
		GeoRaster& ClearMasks() { _Masks.clear(); return *this; }

        GeoRaster& AddFunction(GeoFunction func) { _Functions.push_back(func); return *this; }
        GeoRaster& ClearFunctions() { _Functions.clear(); return *this; }

        // Logical functions
		//! Greater than
		GeoRaster operator>(double val) { return GeoRaster(*this, GeoFunction(">",val)); }
		//! Greater than or equal to
		GeoRaster operator>=(double val) { return GeoRaster(*this, GeoFunction(">=",val)); }
		//! Less than
		GeoRaster operator<(double val) { return GeoRaster(*this, GeoFunction("<",val)); }
		//! Less than or equal to
		GeoRaster operator<=(double val) { return GeoRaster(*this, GeoFunction("<=",val)); }
		//! Equal to
		GeoRaster operator==(double val) { return GeoRaster(*this, GeoFunction("==",val)); }
		// Basic math
		//! Add constant to every pixel
		GeoRaster operator+(double val) { return GeoRaster(*this, GeoFunction("+",val)); }
		//! Subtract constant from every pixel
		GeoRaster operator-(double val) { return GeoRaster(*this, GeoFunction("-",val)); }

        // Statistics - should these be stored?
		//double Min() const { return (GetGDALStats())[0]; }
		//double Max() const { return (GetGDALStats())[1]; }
		//double Mean() const { return (GetGDALStats())[2]; }
		//double StdDev() const { return (GetGDALStats())[3]; }
        cimg_library::CImg<float> ComputeStats(bool RAW=false) const;

        // TODO - If RAW then can use GDAL Statistics, but compare speeds
		// Retrieve Statistics
		/*cimg_library::CImg<double> GetGDALStats() const {
			double min, max, mean, stddev;
			_GDALRasterBand->GetStatistics(false, true, &min, &max, &mean, &stddev);
			cimg_library::CImg<double> stats(4);
			stats(0) = min;
			stats(1) = max;
			stats(2) = mean;
			stats(3) = stddev;
			return stats;
		}*/
		// Compute Statistics
		/*cimg_library::CImg<double> ComputeGDALStats() const {
			double min, max, mean, stddev;
			_GDALRasterBand->ComputeStatistics(false, &min, &max, &mean, &stddev, NULL, NULL);
			cimg_library::CImg<double> stats(4);
			stats(0) = min;
			stats(1) = max;
			stats(2) = mean;
			stats(3) = stddev;
			return stats;
		}*/
	protected:
        // TODO - examine why not shared pointer? (I think because it's managed by GDALDataset class)
		//! GDALRasterBand
		GDALRasterBand* _GDALRasterBand;

        //! Vector of masks to apply
		std::vector< GeoRaster > _Masks;

		//! Bool if nodata value is used
		bool _NoData;

		//! Number of valid pixels
		//long _ValidSize;

        // Constants
        int _minDC;
        int _maxDC;
		double _K1;
		double _K2;
		//! in-band exo-atmospheric solar irradiance
		double _Esun;

		gip::Atmosphere _Atmosphere;

		//! List of processing functions to apply on reads (in class GeoProcess)
		std::vector<GeoFunction> _Functions;

	private:
		//! Default constructor - private so not callable
		explicit GeoRaster() {}

		//! Load band from GDALDataset
		void LoadBand(int bandnum=1) {
		    // TODO - Do i need to reset GDALDataset?   Maybe it needs to be done here...
		    // In practice this is called right after GDALDataset is set, so not needed
			_GDALRasterBand = _GDALDataset->GetRasterBand(bandnum);
			int pbSuccess(0);
			_NoData = false;
			_GDALRasterBand->GetNoDataValue(&pbSuccess);
			if (pbSuccess != 0) {
				if (pbSuccess == 1) _NoData = true;
			}
            Chunk();
		}
	};
}
#endif
