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

    template<typename T> int numpytype() {
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
    }

    // Convert CImg into numpy array
    template<typename T> PyObject* CImgToArr(cimg_library::CImg<T> cimg) {
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

        npy_intp dims[] = { cimg.spectrum(), cimg.depth(), cimg.height(), cimg.width() };
        PyObject* arr;
        int numdim = 4;
        if (cimg.spectrum() == 1) {
            numdim = 3;
            if (cimg.depth() == 1) numdim=2;
        } 
        arr = PyArray_SimpleNew(numdim, &dims[4-numdim], typenum);
        //if (dims[0] == 1)
        //    arr = PyArray_SimpleNew(numdim-1,&dims[1], typenum);
        //else arr = PyArray_SimpleNew(numdim, dims, typenum);
        void *arr_data = PyArray_DATA((PyArrayObject*)arr);
        memcpy(arr_data, cimg.data(), PyArray_ITEMSIZE((PyArrayObject*) arr) * dims[0] * dims[1] * dims[2] * dims[3]);
        return arr;
    }

    // Convert numpy array into CImg...currently 2-D only
    template<typename T> cimg_library::CImg<T> ArrToCImg(PyObject* arr) {
        PyArrayObject* _arr = (PyArrayObject*)arr;
        cimg_library::CImg<T> cimg((T*)_arr->data, _arr->dimensions[1], _arr->dimensions[0]);
        return cimg;
    }

%}

%init %{
    // Not really sure what this does or why it's needed
    import_array();
%}

// STL bindings
%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(vectors) std::vector<std::string>;
}
// TODO - add in map<string,string> to/from dictionary

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
%typemap (in) cimg_library::CImg<char> { $1 = ArrToCImg<char>($input); }
%typemap (in) cimg_library::CImg<unsigned short> { 1 = ArrToCImg<unsigned short>($input); }
%typemap (in) cimg_library::CImg<short> { $1 = ArrToCImg<short>($input); }
%typemap (in) cimg_library::CImg<unsigned int> { $1 = ArrToCImg<unsigned int>($input); }
%typemap (in) cimg_library::CImg<int> { $1 = ArrToCImg<int>($input); }
%typemap (in) cimg_library::CImg<unsigned long>& { $1 = ArrToCImg<unsigned long>($input); }
%typemap (in) cimg_library::CImg<long> { $1 = ArrToCImg<long>($input); }
%typemap (in) cimg_library::CImg<float> { $1 = ArrToCImg<float>($input); }
%typemap (in) cimg_library::CImg<double> { $1 = ArrToCImg<double>($input); }

// TODO - Was trying to quiet warnings...didn't work
//%typemap(typecheck) PyArrayObject * = cimg_library::CImg<unsigned char> ;

// GIP functions to ignore (suppresses warnings)
// These operators are redefined below
%ignore gip::GeoData::operator=;
%ignore gip::Colors::operator[];
%ignore gip::GeoImage::operator[];
%ignore gip::GeoRaster::operator==;

// GIP headers and classes to be wrapped
%include "gip/GeoData.h"
%include "gip/GeoRaster.h"
%include "gip/GeoImage.h"
%include "gip/GeoAlgorithms.h"
// TODO - Not sure this really needs to be wrapped
%include "gip/Colors.h"
// TODO - Not sure this really needs to be wrapped
%include "gip/Atmosphere.h"
//%include "gip/GeoRasterIO.h"
//%include "gip/GeoImageIO.h"

// TODO - improve enums.  C++0x scoped enums ?
enum GDALDataType { GDT_Unknown, GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32,
    GDT_Float32, GDT_Float64 };
    #GDT_CInt16, GDT_CInt32, GDT_CFloat32, GDT_Float64
//enum UNITS { RAW, RADIANCE, REFLECTIVITY };

namespace gip {

    // Register file formats with GDAL
    void reg();

    // Just wrapping basic options.
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

    %extend Colors {
        int __getitem__(std::string col) {
            return self->Colors::operator[](col);
        }
        std::string __getitem__(int col) {
            return self->Colors::operator[](col);
        }
    }

    // This renaming should have worked, but didn't.  Used extend instead 
    //%rename(__getitem__) operator[];
    //%rename(__assign__) GeoImage::operator=;
    //%rename(__assign__) GeoRaster::operator=;

    %extend GeoRaster {
        // Processing functions
        GeoRaster __eq__(double val) {
            return self->operator==(val);
        }
        PyObject* read() {
            if (self->Gain() == 1.0 && self->Offset() == 0.0) {
                switch(self->DataType()) {
                    case 1: return CImgToArr(GeoRasterIO<unsigned char>(*self).Read());
                    case 2: return CImgToArr(GeoRasterIO<unsigned int>(*self).Read());
                    case 3: return CImgToArr(GeoRasterIO<int>(*self).Read());
                    case 4: return CImgToArr(GeoRasterIO<unsigned long>(*self).Read());
                    case 5: return CImgToArr(GeoRasterIO<long>(*self).Read());
                    case 6: return CImgToArr(GeoRasterIO<float>(*self).Read());
                    case 7: return CImgToArr(GeoRasterIO<double>(*self).Read());
                    default: throw(std::exception());
                }
            }
            return CImgToArr(GeoRasterIO<float>(*self).Read());
        }
        GeoRaster& write(PyObject* arr) {
            switch(((PyArrayObject*)arr)->descr->type_num) {
                case NPY_UINT8: GeoRasterIO<unsigned char>(*self).Write(ArrToCImg<unsigned char>(arr));
                case NPY_UINT16: GeoRasterIO<unsigned int>(*self).Write(ArrToCImg<unsigned int>(arr));
                case NPY_INT16: GeoRasterIO<int>(*self).Write(ArrToCImg<int>(arr));
                case NPY_UINT32: GeoRasterIO<unsigned long>(*self).Write(ArrToCImg<unsigned long>(arr));
                case NPY_INT32: GeoRasterIO<long>(*self).Write(ArrToCImg<long>(arr));
                case NPY_FLOAT32: GeoRasterIO<float>(*self).Write(ArrToCImg<float>(arr));
                case NPY_FLOAT64: GeoRasterIO<double>(*self).Write(ArrToCImg<double>(arr));
                default: throw(std::exception());
            }
            return *self;
        }

    }

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

    // These are no longer needed, since GeoImage and GeoRaster have
    // wrapper functions to the GeoImageIO and GeoRasterIO functions
    /*
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
    */
}




