#ifndef GIP_GEORASTER_H
#define GIP_GEORASTER_H

#include <boost/filesystem.hpp>
#include <boost/function.hpp>

#include <gip/gip_CImg.h>
#include <gip/GeoData.h>
#include <gip/Function.h>
#include <boost/bind.hpp>
#include <gip/Atmosphere.h>

#include <iostream>
#include <iomanip>

// for tolowercase
#include <boost/algorithm/string.hpp>

#include <typeinfo>

namespace gip {
    //typedef bgeo::model::d2::point_xy<float> point;

    //enum UNITS {RAW, RADIANCE, REFLECTIVITY};
    //enum READMODE {NORMAL, RAW, ATMOSPHERE};

    /*class GeoFunction {
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
    };*/

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
            : GeoData(geodata), _NoData(false), _ValidStats(false), _UnitsOut(""),
            _minDC(1), _maxDC(255), _K1(0), _K2(0), _Esun(0),
            _Atmosphere() {
            LoadBand(bandnum);
        }
        //! Copy constructor
        GeoRaster(const GeoRaster& image, Function func=Function());
        // Copy and add new function
        //GeoRaster(const GeoRaster& image, GeoFunction func);
        //! Assignment Operator
        GeoRaster& operator=(const GeoRaster& image);
        //! Destructor
        ~GeoRaster() {}

        //! \name File Information
        std::string Basename() const { return GeoData::Basename() + "[" + Description() + "]"; }
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
        //! Set units out
        void SetUnitsOut(std::string units) const { _UnitsOut = units; } //return *this; }
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
        //! Copy all meta data from another raster
        GeoRaster& CopyMeta(const GeoRaster& img) {
            SetDescription(img.Description());
            SetUnits(img.Units());
            SetGain(img.Gain());
            SetOffset(img.Offset());
            SetNoData(img.NoDataValue());
            //_GDALRasterBand->SetMetadata(img._GDALRasterBand->GetMetadata());
            return *this;
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

        //! Adds a mask band (1 for valid), applied on read
        const GeoRaster& AddMask(const GeoRaster& band) const {
            _ValidStats = false;
            _Masks.push_back(band);
            return *this;
        }
        //! Remove all masks from band
        const GeoRaster& ClearMasks() const {
            if (!_Masks.empty()) _ValidStats = false;
            _Masks.clear();
            return *this;
        }

        /*GeoRaster& AddFunction(Function func) {
            _ValidStats = false;
            _Functions.push_back(func);
            return *this;
        }
        GeoRaster& ClearFunctions() {
            if (!_Functions.empty()) _ValidStats = false;
            _Functions.clear();
            return *this;
        }*/

        //! \name Processing functions

        GeoRaster operator>(const double &val) const {
            return GeoRaster(*this, Function(boost::bind(&CImg<double>::threshold, _1, val, false, true)));
        }
/*
            if (iFunc->Function() == ">") {
                img.threshold(iFunc->Operand(),false,true);
            } else if (iFunc->Function() == ">=") {
                img.threshold(iFunc->Operand(),false,false);
            } else if (iFunc->Function() == "<") {
                img.threshold(iFunc->Operand(),false,false)^=1;
            } else if (iFunc->Function() == "<=") {
                img.threshold(iFunc->Operand(),false,true)^=1;
            } else if (iFunc->Function() == "==") {
                cimg_for(img,ptr,T) if (*ptr == iFunc->Operand()) *ptr = 1; else *ptr = 0;
                //img = img.get_threshold(iFunc->Operand(),false,false) - img.get_threshold(iFunc->Operand(),false,true);
            } else if (iFunc->Function() == "+") {
                img = img + iFunc->Operand();
            } else if (iFunc->Function() == "-") {
                img = img - iFunc->Operand();
*/

        // Logical functions
        //! Greater than
        //GeoRaster operator>(double val) const { return GeoRaster(*this, GeoFunction(">",val)); }
        //! Greater than or equal to
        //GeoRaster operator>=(double val) const { return GeoRaster(*this, GeoFunction(">=",val)); }
        //! Less than
        //GeoRaster operator<(double val) const { return GeoRaster(*this, GeoFunction("<",val)); }
        //! Less than or equal to
        //GeoRaster operator<=(double val) const { return GeoRaster(*this, GeoFunction("<=",val)); }
        //! Equal to
        //GeoRaster operator==(double val) const { return GeoRaster(*this, GeoFunction("==",val)); }
        // Basic math
        //! Add constant to every pixel
        //GeoRaster operator+(double val) const { return GeoRaster(*this, GeoFunction("+",val)); }
        //! Subtract constant from every pixel
        //GeoRaster operator-(double val) const { return GeoRaster(*this, GeoFunction("-",val)); }

        // Statistics - should these be stored?
        //double Min() const { return (GetGDALStats())[0]; }
        //double Max() const { return (GetGDALStats())[1]; }
        //double Mean() const { return (GetGDALStats())[2]; }
        //double StdDev() const { return (GetGDALStats())[3]; }
        cimg_library::CImg<float> ComputeStats() const;

        cimg_library::CImg<float> Histogram(int bins=100, bool cumulative=false) const;

        float Percentile(float p) const;

        // TODO - If RAW then can use GDAL Statistics, but compare speeds
        // Compute Statistics
        /*cimg_library::CImg<double> ComputeGDALStats() const {
            double min, max, mean, stddev;
            _GDALRasterBand->GetStatistics(false, true, &min, &max, &mean, &stddev);
            _GDALRasterBand->ComputeStatistics(false, &min, &max, &mean, &stddev, NULL, NULL);
            cimg_library::CImg<double> stats(4);
            stats(0) = min;
            stats(1) = max;
            stats(2) = mean;
            stats(3) = stddev;
            return stats;
        }*/

        //! \name File I/O
        template<class T> cimg_library::CImg<T> ReadRaw(int chunknum=0) const;
        template<class T> cimg_library::CImg<T> ReadRaw(iRect chunk) const;
        template<class T> cimg_library::CImg<T> Read(int chunknum=0) const;
        template<class T> cimg_library::CImg<T> Read(iRect chunk) const;
        template<class T> GeoRaster& WriteRaw(cimg_library::CImg<T> img, int chunknum=0);
        template<class T> GeoRaster& WriteRaw(cimg_library::CImg<T> img, iRect chunk);
        template<class T> GeoRaster& Write(cimg_library::CImg<T> img, int chunknum=0);
        template<class T> GeoRaster& Process(const GeoRaster& raster);

         //! Get Saturation mask: 1's where it's saturated
        cimg_library::CImg<unsigned char> SaturationMask(int chunk=0) const {
            switch (DataType()) {
                case GDT_Byte: return _Mask<unsigned char>(_maxDC, chunk);
                case GDT_UInt16: return _Mask<unsigned short>(_maxDC, chunk);
                case GDT_Int16: return _Mask<short>(_maxDC, chunk);
                case GDT_UInt32: return _Mask<unsigned int>(_maxDC, chunk);
                case GDT_Int32: return _Mask<int>(_maxDC, chunk);
                case GDT_Float32: return _Mask<float>(_maxDC, chunk);
                case GDT_Float64: return _Mask<double>(_maxDC, chunk);
                default: return _Mask<double>(_maxDC, chunk);
            }
        }

        //! NoData mask: 1's where it's good data
        cimg_library::CImg<unsigned char> NoDataMask(int chunk=0) const {
            // TODO - if NoData not set, return all 1s
            if (!NoData()) {
                int width, height;
                if (chunk == 0) {
                    width = XSize();
                    height = YSize();
                } else {
                    iRect ch = _PadChunks[chunk-1];
                    width = ch.x1()-ch.x0()+1;
                    height = ch.y1()-ch.y0()+1;
                }
                return CImg<unsigned char>(width,height,1,1,0);
            }
            switch (DataType()) {
                case GDT_Byte: return _Mask<unsigned char>(NoDataValue(), chunk);
                case GDT_UInt16: return _Mask<unsigned short>(NoDataValue(), chunk);
                case GDT_Int16: return _Mask<short>(NoDataValue(), chunk);
                case GDT_UInt32: return _Mask<unsigned int>(NoDataValue(), chunk);
                case GDT_Int32: return _Mask<int>(NoDataValue(), chunk);
                case GDT_Float32: return _Mask<float>(NoDataValue(), chunk);
                case GDT_Float64: return _Mask<double>(NoDataValue(), chunk);
                default: return _Mask<double>(NoDataValue(), chunk);
            }
        }

    protected:
        // TODO - examine why not shared pointer? (I think because it's managed by GDALDataset class)
        //! GDALRasterBand
        GDALRasterBand* _GDALRasterBand;

        //! Vector of masks to apply
        mutable std::vector< GeoRaster > _Masks;

        //! Bool if nodata value is used
        bool _NoData;

        //! Valid Stats Flag
        mutable bool _ValidStats;
        //! Statistics
        mutable CImg<double> _Stats;

        mutable std::string _UnitsOut;

        // Constants
        int _minDC;
        int _maxDC;
        double _K1;
        double _K2;
        //! in-band exo-atmospheric solar irradiance
        double _Esun;

        gip::Atmosphere _Atmosphere;

        //! List of processing functions to apply on reads (in class GeoProcess)
        std::vector< boost::function< CImg<double>& > > _Functions;

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

        template<class T> inline cimg_library::CImg<unsigned char> _Mask(T val, int chunk=0) const {
            using cimg_library::CImg;
            CImg<T> img = ReadRaw<T>(chunk);
            CImg<unsigned char> mask(img.width(),img.height(),1,1,0);
            cimg_forXY(img,x,y) if (img(x,y) == val) mask(x,y) = 1;
            return mask;
        }

    }; //class GeoImage

    //! \name File I/O
    template<class T> cimg_library::CImg<T> GeoRaster::ReadRaw(int chunknum) const {
        if (chunknum == 0)
            return ReadRaw<T>( iRect(iPoint(0,0),iPoint(XSize()-1,YSize()-1)) );
        //if (chunknum == _chunknum) return _cimg;
        //_chunknum = chunknum;
        //_cimg.assign( ReadRaw( _PadChunks[chunknum-1] ) );
        //return _cimg;
        return ReadRaw<T>( _PadChunks[chunknum-1] ) ;
    }

    //! Read raw chunk given bounding box
    template<class T> cimg_library::CImg<T> GeoRaster::ReadRaw(iRect chunk) const {
        // This doesn't check for in bounds, should it?
        int width = chunk.x1()-chunk.x0()+1;
        int height = chunk.y1()-chunk.y0()+1;

        T* ptrPixels = new T[width*height];
        CPLErr err = _GDALRasterBand->RasterIO(GF_Read, chunk.x0(), chunk.y0(), width, height, 
            ptrPixels, width, height, type2GDALtype(typeid(T)), 0, 0);
        if (err != CE_None) {
            std::stringstream err;
            err << "error reading " << CPLGetLastErrorMsg();
            throw std::runtime_error(err.str());
        }
        CImg<T> img(ptrPixels,width,height);

        // Apply all masks TODO - cmask need to be float ?
        if (_Masks.size() > 0) {
            if (Options::Verbose() > 3 && (chunk.p0()==iPoint(0,0)))
                std::cout << Basename() << ": Applying " << _Masks.size() << " masks" << std::endl;
            CImg<float> cmask(_Masks[0].Read<float>(chunk));
            for (unsigned int i=1; i<_Masks.size(); i++) {
                cmask.mul(_Masks[i].Read<float>(chunk));
            }
            cimg_forXY(img,x,y) {
                if (cmask(x,y) != 1) img(x,y) = NoDataValue();
            }
        }

        if (Options::Verbose() > 3) std::cout << Basename() << ": read " << chunk << std::endl;

        delete ptrPixels;
        return img;
    }

    //! Retrieve a piece of the image as a CImg
    template<class T> cimg_library::CImg<T> GeoRaster::Read(int chunknum) const {
        if (chunknum == 0)
            return Read<T>( iRect(iPoint(0,0),iPoint(XSize()-1,YSize()-1)) );
        return Read<T>( this->_PadChunks[chunknum-1]);
    }

    //! Retrieve a piece of the image as a CImg
    template<class T> cimg_library::CImg<T> GeoRaster::Read(iRect chunk) const {
        using cimg_library::CImg;

        CImg<T> img(ReadRaw<T>(chunk));
        CImg<T> imgorig(img);

        bool updatenodata = false;
        // Convert data to radiance (if not raw requested)
        if (Gain() != 1.0 || Offset() != 0.0) {
            img = Gain() * (img-_minDC) + Offset();
            updatenodata = true;
        }
        // apply atmosphere if there is one (which would data is radiance units) TODO - check units
        if (Atmosphere()) {
            if (Options::Verbose() > 3 && (chunk.p0()==iPoint(0,0)))
                std::cout << Basename() << ": applying atmosphere" << std::endl;
            double e = (Thermal()) ? 0.95 : 1;  // For thermal band, currently water only
            img = (img - (_Atmosphere.Lu() + (1-e)*_Atmosphere.Ld())) / (_Atmosphere.t() * e);
            updatenodata = true;
        }

        // Convert to reflectance
        if ((Units() == "radiance") && (_UnitsOut == "reflectance")) {
            if (Options::Verbose() > 3 && (chunk.p0()==iPoint(0,0)))
                std::cout << Basename() << ": converting radiance to reflectance" << std::endl;
            if (Thermal()) {
                cimg_for(img,ptr,T) *ptr = (_K2/log(_K1/(*ptr)+1)) - 273.15;
            } else {
                float normrad = Atmosphere() ? (1.0/_Atmosphere.Ld()) : (1.0/_Esun);
                cimg_for(img,ptr,T) *ptr = *ptr * normrad;
            }
            updatenodata = true;
        }

        std::vector<Function>::const_iterator iFunc;
        for (iFunc=_Functions.begin();iFunc!=_Functions.end();iFunc++) {
            if (Options::Verbose() > 3 && (chunk.p0()==iPoint(0,0)))
                std::cout << Basename() << ": Applying function " << iFunc->Function() << std::endl;
            ifunc->Function()(img);
            updatenodata = true;
        }

        // Apply Processing functions

        // If processing was applied update NoData values where needed
        if (NoData() && updatenodata) {
            cimg_forXY(img,x,y) {
                if (imgorig(x,y) == NoDataValue()) img(x,y) = NoDataValue();
            }
        }
        return img;
    }

    //! Write raw CImg to file
    template<class T> GeoRaster& GeoRaster::WriteRaw(cimg_library::CImg<T> img, int chunknum) {
        if (chunknum == 0)
            return WriteRaw(img, iRect(iPoint(0,0),iPoint(XSize()-1,YSize()-1)) );
        return WriteRaw(img, _Chunks[chunknum-1] );
    }

    //! Write raw CImg to file
    template<class T> GeoRaster& GeoRaster::WriteRaw(cimg_library::CImg<T> img, iRect chunk) {
        if (Options::Verbose() > 3) {
            std::cout << Basename() << ": writing " << img.width() << " x " 
                << img.height() << " image to rect " << chunk << std::endl;
        }
        CPLErr err = _GDALRasterBand->RasterIO(GF_Write, chunk.x0(), chunk.y0(), 
            chunk.width(), chunk.height(), img.data(), img.width(), img.height(), 
            type2GDALtype(typeid(T)), 0, 0);
        if (err != CE_None) {
            std::stringstream err;
            err << "error writing " << CPLGetLastErrorMsg();
            throw std::runtime_error(err.str());
        }
        _ValidStats = false;
        return *this;
    }

    //! Write a Cimg to the file
    template<class T> GeoRaster& GeoRaster::Write(cimg_library::CImg<T> img, int chunknum) {
        iRect chunk;
        if (chunknum == 0) {
            chunk = iRect( iPoint(0,0), iPoint(XSize()-1,YSize()-1) );
        } else {
            chunk = _Chunks[chunknum-1];
            iRect pchunk = _PadChunks[chunknum-1];
            if (chunk != pchunk) {
                iPoint p0(chunk.p0()-pchunk.p0());
                iPoint p1 = p0 + iPoint(chunk.width()-1,chunk.height()-1);
                img.crop(p0.x(),p0.y(),p1.x(),p1.y());
            }
        }
        if (Gain() != 1.0 || Offset() != 0.0) {
            cimg_for(img,ptr,T) if (*ptr != NoDataValue()) *ptr = (*ptr-Offset())/Gain();
        }
        if (Options::Verbose() > 3 && (chunk.p0()==iPoint(0,0)))
            std::cout << Basename() << ": Writing (" << Gain() << "x + " << Offset() << ")" << std::endl;
        /*if (BadValCheck) {
            cimg_for(img,ptr,T) if ( std::isinf(*ptr) || std::isnan(*ptr) ) *ptr = NoDataValue();
        }*/
        return WriteRaw(img,chunk); 
    }

    //! Process input band into this
    template<class T> GeoRaster& GeoRaster::Process(const GeoRaster& raster) {
        using cimg_library::CImg;
        for (unsigned int iChunk=1; iChunk<=NumChunks(); iChunk++) {
                CImg<T> cimg = raster.Read<T>(iChunk);
                //WriteChunk(CImg<T>().assign(cimg.round()),*iChunk, RAW);
                Write(cimg,iChunk); //, RAW);
        }
        GDALRasterBand* band = raster.GetGDALRasterBand();
        CopyCategoryNames(raster);
        _GDALRasterBand->SetDescription(band->GetDescription());
        _GDALRasterBand->SetColorInterpretation(band->GetColorInterpretation());
        _GDALRasterBand->SetMetadata(band->GetMetadata());
        CopyCoordinateSystem(raster);
        return *this;
    }

} // namespace GIP

#endif
