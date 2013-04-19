%module gippylib

%{
//#include "gip/Colors.h"
#include "gip/Options.h"
//#include "gip/GeoData.h"
//#include "gip/GeoRaster.h"
#include "gip/GeoImage.h"
#include "gip/GeoAlgorithms.h"
//#include "gdal/gdal_priv.h"

using namespace gip;

namespace gip {
    void reg() { GDALAllRegister(); }
}
%}

// STL bindings
%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(vectors) std::vector<std::string>;
}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// GIP functions to ignore (suppresses warnings)
//%ignore gip::Options::operator=;
%ignore gip::cimg_library::CImg;
%ignore gip::GeoData::operator=;
%ignore gip::Colors::operator[];
%ignore gip::GeoImage::operator[];
//%ignore boost::program_options::options_description;

// GIP headers and classes to be wrapped
%include "gip/Colors.h"
%include "gip/Sensor.h"
%include "gip/Atmosphere.h"
//%include "gip/Options.h"
%include "gip/GeoData.h"
%include "gip/GeoRaster.h"
%include "gip/GeoImage.h"
%include "gip/GeoAlgorithms.h"

enum GDALDataType { GDT_Unknown, GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32,
    GDT_Float32, GDT_Float64 };
    #GDT_CInt16, GDT_CInt32, GDT_CFloat32, GDT_Float64
// work-around, should support scoped enums for C++0x
enum UNITS { RAW, RADIANCE, REFLECTIVITY };

namespace gip {

    // Register file formats with GDAL
    void reg();
    // Wrap up basic options
    class Options {
    public:
        //static std::string ConfigDir();
        //static void SetConfigDir(std::string dir);
        static std::string DefaultFormat();
        static void SetDefaultFormat(std::string format);
        static float ChunkSize();
        static void SetChunkSize(float sz);
        static int Verbose();
        static void SetVerbose(int v);
    };

    //%rename(__getitem__) operator[];

    %extend Colors {
        int __getitem__(std::string col) {
            return self->Colors::operator[](col);
        }
        std::string __getitem__(int col) {
            return self->Colors::operator[](col);
        }
    }

    //%rename(__assign__) GeoImage::operator=;
    //%rename(__assign__) GeoRaster::operator=;

    %extend GeoImage {
        GeoRaster __getitem__(std::string col) {
            return self->GeoImage::operator[](col);
        }
        GeoRaster __getitem__(int band) {
            return self->GeoImage::operator[](band);
        }
        GeoRaster& __setitem__(int band, const GeoRaster& raster) {
            self->operator[](band) = raster;
            return self->GeoImage::operator[](band);
        }
        GeoRaster& __setitem__(std::string col, const GeoRaster& raster) {
            self->operator[](col) = raster;
            return self->operator[](col);
        }
    }

}
