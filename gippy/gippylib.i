%module gippylib
%{
    #define SWIG_FILE_WITH_INIT
    //#include "gip/Colors.h"
    //#include "gip/GeoData.h"
    //#include "gip/GeoRaster.h"
    #include "gip/GeoImage.h"
    #include "gip/GeoAlgorithms.h"
    //#include "gdal/gdal_priv.h"
    #include <python2.7/Python.h>
    #include <numpy/arrayobject.h>
    #include <iostream>
    #include "gip/gip_CImg.h"
    #include "gip/GeoRasterIO.h"
    #include "gip/GeoImageIO.h"

    using namespace gip;

    namespace gip {
        void reg() { GDALAllRegister(); }
    }

    template<typename T> PyObject* CImgToArr(cimg_library::CImg<T> cimg) {
        npy_intp dims[] = { cimg.height(), cimg.width() };
        int typenum;
        if (typeid(T) == typeid(unsigned char)) typenum = NPY_UINT8;
        else if (typeid(T) == typeid(char)) typenum = NPY_INT8;
        else if (typeid(T) == typeid(unsigned short)) typenum = NPY_UINT16;
        else if (typeid(T) == typeid(short)) typenum = NPY_INT16;
        else if (typeid(T) == typeid(unsigned int)) typenum = NPY_UINT32;
        else if (typeid(T) == typeid(int)) typenum = NPY_INT32;
        else if (typeid(T) == typeid(unsigned long)) typenum = NPY_UINT64;
        else if (typeid(T) == typeid(long)) typenum = NPY_INT64;
        else if (typeid(T) == typeid(float)) typenum = NPY_FLOAT32;
        else if (typeid(T) == typeid(double)) typenum = NPY_FLOAT64;
        else throw(std::exception());
        PyObject* arr;
        if (dims[0] == 1)
            arr = PyArray_SimpleNew(1,&dims[1], typenum);
        else arr = PyArray_SimpleNew(2, dims, typenum);
        void *arr_data = PyArray_DATA((PyArrayObject*)arr);
        memcpy(arr_data, cimg.data(), PyArray_ITEMSIZE((PyArrayObject*) arr) * dims[0] * dims[1]);
        return arr;
    }

    template<typename T> cimg_library::CImg<T> ArrToCImg(PyObject* arr) {
        PyArrayObject* _arr = (PyArrayObject*)arr;
        cimg_library::CImg<T> cimg(_arr->data, _arr->dimensions[0], _arr->dimensions[1]);
        return cimg;
    }

%}

%init %{
    import_array();
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

// CImg -> numpy
%typemap (out) cimg_library::CImg<unsigned char> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<char> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<unsigned short> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<short> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<unsigned int> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<int> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<unsigned long> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<long> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<float> { return CImgToArr($1); }
%typemap (out) cimg_library::CImg<double> { return CImgToArr($1); }

// numpy -> CImg
%typemap (in) cimg_library::CImg<unsigned char> { $1 = ArrToCImg<unsigned char>($input); }

/*%typemap (in) cimg_library::CImg<char> { $1 = ArrToCImg<char>($input); }
%typemap (in) cimg_library::CImg<unsigned short> { 1 = ArrToCImg<unsigned short>($input); }
%typemap (in) cimg_library::CImg<short> { $1 = ArrToCImg<short>($input); }
%typemap (in) cimg_library::CImg<unsigned int> { $1 = ArrToCImg<unsigned int>($input); }
%typemap (in) cimg_library::CImg<int> { $1 = ArrToCImg<int>($input); }
%typemap (in) cimg_library::CImg<unsigned long>& { $1 = ArrToCImg<unsigned long>($input); }
%typemap (in) cimg_library::CImg<long> { $1 = ArrToCImg<long>($input); }
%typemap (in) cimg_library::CImg<float> { $1 = ArrToCImg<float>($input); }
%typemap (in) cimg_library::CImg<double> { $1 = ArrToCImg<double>($input); }
*/
/*
%typemap (in) cimg_library::CImg<unsigned char> {
    PyArrayObject* arr = (PyArrayObject*)$input;
    cimg_library::CImg<unsigned char> cimg(arr->data, arr->dimensions[0], arr->dimensions[1]);
    $1 = cimg;
}
%typemap (in) cimg_library::CImg<unsigned char>& {
    PyArrayObject* arr = (PyArrayObject*)$input;
    cimg_library::CImg<unsigned char> cimg(arr->data, arr->dimensions[0], arr->dimensions[1]);
    $1 = &cimg;
}
%typemap (in) cimg_library::CImg<unsigned char>* {
    cimg_library::CImg<unsigned char> cimg($input->data, $input->dimensions[0], $input->dimensions[1]);
    $1 = &cimg;
}
*/


// GIP functions to ignore (suppresses warnings)
// These operators are redefined below
%ignore gip::GeoData::operator=;
%ignore gip::Colors::operator[];
%ignore gip::GeoImage::operator[];
%ignore gip::GeoRaster::operator==;

// GIP headers and classes to be wrapped
%include "gip/Colors.h"
%include "gip/Atmosphere.h"
%include "gip/GeoData.h"
%include "gip/GeoRaster.h"
%include "gip/GeoImage.h"
%include "gip/GeoAlgorithms.h"
%include "gip/GeoRasterIO.h"
%include "gip/GeoImageIO.h"

enum GDALDataType { GDT_Unknown, GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32,
    GDT_Float32, GDT_Float64 };
    #GDT_CInt16, GDT_CInt32, GDT_CFloat32, GDT_Float64
// work-around, should support scoped enums for C++0x
//enum UNITS { RAW, RADIANCE, REFLECTIVITY };

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
        static std::string WorkDir();
        static void SetWorkDir(std::string workdir);
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

    %extend GeoRaster {
        GeoRaster __eq__(double val) {
            return self->operator==(val);
        }
    }

    %template(GeoRaster_byte) gip::GeoRasterIO<unsigned char>;
    %template(GeoRaster_int16) gip::GeoRasterIO<short int>;
    %template(GeoRaster_int32) gip::GeoRasterIO<int>;
    %template(GeoRaster_int64) gip::GeoRasterIO<long>;
    %template(GeoRaster_float) gip::GeoRasterIO<float>;
    %template(GeoRaster_double) gip::GeoRasterIO<double>;

    %template(GeoImage_byte) gip::GeoImageIO<unsigned char>;
    %template(GeoImage_int16) gip::GeoImageIO<short int>;
    %template(GeoImage_int32) gip::GeoImageIO<int>;
    %template(GeoImage_int64) gip::GeoImageIO<long>;
    %template(GeoImage_float) gip::GeoImageIO<float>;
    %template(GeoImage_double) gip::GeoImageIO<double>;
}




