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
            : GeoData(geodata), _NoData(false), _ValidSize(0),
            _minDC(1), _maxDC(255), _K1(0), _K2(0), _Esun(0),
            _Atmosphere() {
			LoadBand(bandnum);
		}
		//! Copy constructor
		GeoRaster(const GeoRaster& image);
		//! Copy and add new function
		GeoRaster(const GeoRaster& image, GeoFunction func);
		//! Assignment Operator
		GeoRaster& operator=(const GeoRaster& image);
		//! Destructor
		~GeoRaster() {} //_GDALRasterBand->FlushCache(); }

        //! Copy passed raster band into this raster (is this constructor?)
        GeoRaster& Copy(const GeoRaster&, bool RAW=false);

        GeoRaster& AddFunction(GeoFunction func) { _Functions.push_back(func); return *this; }

		//! \name File Information
		//! Band X Size
		unsigned int XSize() const { return _GDALRasterBand->GetXSize(); }
		//! Band Y Size
		unsigned int YSize() const { return _GDALRasterBand->GetYSize(); }
		//! Total Valid Pixels
		unsigned int ValidSize() const { return _ValidSize; }
		//! Set # of total valid pixels
		//void SetValidSize(unsigned int sz) { _ValidSize = sz; }

		// Average wavelength of spectral band (in microns)
		// TODO - Wavelengths ?

		//! Get datatype
		std::string DataTypeStr() const { return GDALGetDataTypeName(DataType()); }
		//! Get GDALDatatype
		GDALDataType DataType() const { return _GDALRasterBand->GetRasterDataType(); }
		//! Get GDALRasterBand object - use cautiously
		GDALRasterBand* GetGDALRasterBand() const { return _GDALRasterBand; }

		//! Output file info
		std::string Info(bool stats=false) const;

		// \name File I/O
		//! Get Band Description
		std::string Description() const {
			return _GDALRasterBand->GetDescription();
		}
		//! Set Band Description
		void SetDescription(std::string desc) {
			// Set description in band
			_GDALRasterBand->SetDescription(desc.c_str());
			// Also set description in dataset metadata since band desc doesn't work at least in GTiff
			//_GeoImage->SetMeta("Band "+to_string(_GDALRasterBand->GetBand()), desc);
		}
		//! Set Color Interp (used in other GDAL aware programs?)
		void SetColor(std::string col) {
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
		//! Set Unit type
		void SetUnits(std::string units) const { _GDALRasterBand->SetUnitType(units.c_str()); }
		//! Get gain
		float Gain() const { return _GDALRasterBand->GetScale(); }
		//! Get offset
		float Offset() const { return _GDALRasterBand->GetOffset(); }
		//! Set gain
		void SetGain(float gain) { _GDALRasterBand->SetScale(gain); }
		//! Set offset
		void SetOffset(float offset) { _GDALRasterBand->SetOffset(offset); }

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
		void SetAtmosphere(gip::Atmosphere atm) { _Atmosphere = atm; }
        //! Clear atmospheric correction
        void ClearAtmosphere() { _Atmosphere = gip::Atmosphere(); }

		// NoData
		//! Flag indicating if NoData value is used or not
		bool NoData() const { return _NoData;
		    /*int pbSuccess(0);
		    _GDALRasterBand->GetNoDataValue(&pbSuccess);
		    std::cout << "pbSuccess = " << pbSuccess << " " << _NoData << std::endl; //return _NoData;
		    if (pbSuccess == 1) return true; else return false;*/
        }
		//! Get NoDataValue
		double NoDataValue() const {
		    //std::cout << "NoDataValue" << std::endl;
			//if (_NoData) return _GDALRasterBand->GetNoDataValue(); else return 0;
			return _GDALRasterBand->GetNoDataValue();
		}
		//! Set No Data value
		void SetNoData(double val) {
		    //std::cout << "SetNoData " << val << std::endl;
			_GDALRasterBand->SetNoDataValue(val);
			_NoData = true;
		}
		//! Set NoDataValue to a default based on datatype of band
		void SetNoData() {
		    //std::cout << "SetNoData" << std::endl;
			// only set if not already set
			if (_NoData) return;
			double val = -32766;
			GDALDataType dt = DataType();
			if (dt == GDT_Byte) val = 255;
			else if (dt == GDT_UInt16 || dt == GDT_UInt32) val = -val;
			SetNoData(val);
		}
		//! Clear NoData
		void ClearNoData() {
		    //std::cout << "ClearNoData" << std::endl;
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

        //! Adds a mask band (1 for valid), applied on read
		GeoRaster& AddMask(const GeoRaster& band) { _Masks.push_back(band); return *this; }

        //! Statistics - should these be stored?
		double Min() const { return (GetGDALStats())[0]; }
		double Max() const { return (GetGDALStats())[1]; }
		double Mean() const { return (GetGDALStats())[2]; }
		double StdDev() const { return (GetGDALStats())[3]; }

		//! Retrieve Statistics
		cimg_library::CImg<double> GetGDALStats() const {
			double min, max, mean, stddev;
			_GDALRasterBand->GetStatistics(false, true, &min, &max, &mean, &stddev);
			cimg_library::CImg<double> stats(4);
			stats(0) = min;
			stats(1) = max;
			stats(2) = mean;
			stats(3) = stddev;
			return stats;
		}

		//! Compute Statistics
		cimg_library::CImg<double> ComputeGDALStats() const {
			double min, max, mean, stddev;
			_GDALRasterBand->ComputeStatistics(false, &min, &max, &mean, &stddev, NULL, NULL);
			cimg_library::CImg<double> stats(4);
			stats(0) = min;
			stats(1) = max;
			stats(2) = mean;
			stats(3) = stddev;
			return stats;
		}

		cimg_library::CImg<float> ComputeStats(bool RAW=false) const;

		// \name Processing functions
		//! Greater than
		GeoRaster operator>(double val) const { return GeoRaster(*this, GeoFunction(">",val)); }
		GeoRaster operator>=(double val) const { return GeoRaster(*this, GeoFunction(">=",val)); }
		GeoRaster operator<(double val) const { return GeoRaster(*this, GeoFunction("<",val)); }
		GeoRaster operator<=(double val) const { return GeoRaster(*this, GeoFunction("<=",val)); }
		GeoRaster operator==(double val) const { return GeoRaster(*this, GeoFunction("==",val)); }
		GeoRaster operator+(double val) const { return GeoRaster(*this, GeoFunction("+",val)); }
		GeoRaster operator-(double val) const { return GeoRaster(*this, GeoFunction("-",val)); }

	protected:
		//! GDALRasterBand
		GDALRasterBand* _GDALRasterBand;

        //! Vector of masks to apply
		std::vector< GeoRaster > _Masks;

		//! Bool if nodata value is used
		bool _NoData;

		//! Number of valid pixels
		long _ValidSize;

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
			_GDALRasterBand = _GDALDataset->GetRasterBand(bandnum);
			int pbSuccess(0);
			_NoData = false;
			_GDALRasterBand->GetNoDataValue(&pbSuccess);
			if (pbSuccess != 0) {
				if (pbSuccess == 1) _NoData = true;
			}
		}

        // Set Reflectance Coefficient given DOY, solar irradiance and solar zenith
        /*void SetRefCoef(int doy, float E, float zenith) {
            float theta = M_PI * zenith/180.0;
            float sundist = (1.0 - 0.016728 * cos(M_PI * 0.9856 * (doy-4.0)/180.0));
            _RefCoef = (E == 0.0) ? 0.0 : (M_PI * sundist * sundist) / (E * cos(theta));
            //std::cout << "doy " << doy << std::endl;
            //std::cout << "E " << E << std::endl;
            //std::cout << "theta " << theta << std::endl;
            //std::cout << "sundist " << sundist << std::endl;
            //std::cout << "coef " << _RefCoef << std::endl;
        }*/

	};
}
#endif
