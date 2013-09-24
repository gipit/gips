#ifndef GIP_GEORASTERIO_H
#define GIP_GEORASTERIO_H

#include <exception>
//#include <gip/gip.h>
#include <gip/GeoRaster.h>
#include <cmath>

// only used for debugging
#include <iostream>

namespace gip {
    using cimg_library::CImg;

    // TODO - Compile with C++0x to use scoped enums
    //enum UNITS {RAW, RADIANCE, REFLECTIVITY};

	template<class T> class GeoRasterIO : public GeoRaster {
	public:
        typedef boost::geometry::model::box<point> bbox;
		// Constructors/Destructor
		GeoRasterIO(GeoRaster& img)
			: GeoRaster(img) {}
		GeoRasterIO(const GeoRaster& img)
			: GeoRaster(img) {}
		//GeoRasterIO(GeoRasterIO& img) {}
		~GeoRasterIO() {}

		//! \name File I/O
		//! Retrieve a piece of the image as a CImg
		CImg<T> Read(bbox chunk, UNITS units=RAW) const {
			point p1 = chunk.min_corner();
			point p2 = chunk.max_corner();
			int width = p2.x()-p1.x()+1;
			int height = p2.y()-p1.y()+1;
			T* ptrPixels = new T[width*height];
			//cout << "Reading " << bgeo::dsv(p1) << bgeo::dsv(p2) << " w=" << width << " h=" << height << endl;
			// casting away const, safe because this is a read-only const_cast<GDALRasterBand*>
			CPLErr err = _GDALRasterBand->RasterIO(GF_Read, p1.x(), p1.y(), width, height, ptrPixels, width, height, GDALType(), 0, 0);
			if (err != CE_None) std::cout << "Error reading " << Filename() << ": " << CPLGetLastErrorMsg() << std::endl;
			CImg<T> img(ptrPixels, width, height);
			CImg<T> imgorig(img);

			// Convert data to radiance
			if (units != RAW) {
                if (Gain() != 1.0 || Offset() != 0.0) {
                    img = Gain() * (img-_minDC) + Offset();
                }
			}
			// Coef image
			//CImg<T> coefimg(img);

            // Atmospheric correction to surface-leaving radiance
            double normrad;
            if (Atmosphere()) {
                // For thermal band, currently water only
                double e = (Thermal()) ? 0.95 : 1;
                img = (img - (_Atmosphere.Lu() + (1-e)*_Atmosphere.Ld())) / (_Atmosphere.t() * e);
                normrad = 1.0/_Atmosphere.Ld();
            } else {
                normrad = 1.0/_Esun;
                //std::cout << "normrad, rad " << normrad << ", " << _Sensor.Radiance() << ", " << _Sensor.RefCoef() << std::endl;
            }

            // Convert to reflectance
			if (units == REFLECTIVITY) {
                // Convert to reflectance or temperature here
                if (Thermal()) {
                    cimg_for(img,ptr,T) *ptr = (_K2/log(_K1/(*ptr) + 1)) - 273.15;
                    //img = _Sensor.K2()/(((_Sensor.K1()/img) + 1).log());
                } else {
                    cimg_for(img,ptr,T) *ptr = *ptr * normrad; //* _Sensor.RefCoef();
                }
			}

			// Testing 6s coefficients
			/*if (Atmosphere() && units == REFLECTIVITY && _Atmosphere.Coef()) {
			    double f;
                cimg_for(coefimg,ptr,T) {
                    f = _Atmosphere.Xa * *ptr - _Atmosphere.Xb;
                    *ptr = f/(1 + _Atmosphere.Xc * f);
                }
                CImg<T> stats = (img - coefimg).get_stats();
                std::cout << "Min/Max = " << stats(0) << ", " << stats(1)
                    << " Mean/StdDev = " << stats(2) << " +/- " << stats(3) << std::endl;
                //img = coefimg;
			}*/

			// Apply Processing functions
			std::vector<GeoFunction>::const_iterator iFunc;
			for (iFunc=_Functions.begin();iFunc!=_Functions.end();iFunc++) {
				//std::cout << Basename() << ": Applying Function " << iFunc->Function() << " " << iFunc->Operand() << std::endl;
				if (iFunc->Function() == ">") {
					img.threshold(iFunc->Operand(),false,true);
				} else if (iFunc->Function() == ">=") {
					img.threshold(iFunc->Operand(),false,false);
				} else if (iFunc->Function() == "<") {
					img.threshold(iFunc->Operand(),false,false)^=1;
				} else if (iFunc->Function() == "<=") {
					img.threshold(iFunc->Operand(),false,true)^=1;
				} else if (iFunc->Function() == "+") {
					img = img + iFunc->Operand();
				} else if (iFunc->Function() == "-") {
					img = img - iFunc->Operand();
				}
			}

			// Apply all masks
			if (_Masks.size() > 0) {
                GeoRasterIO<unsigned char> mask(_Masks[0]);
                CImg<unsigned char> cmask(mask.Read(chunk));
                for (unsigned int i=1; i<_Masks.size(); i++) {
                    mask = GeoRasterIO<unsigned char>(_Masks[i]);
                    cmask.mul(mask.Read(chunk));
                }
                img.mul(cmask);
			}

			// If processing was apply NoData values where needed
			if (NoData() && units != RAW) {
                cimg_forXY(img,x,y) {
                    if (imgorig(x,y) == NoDataValue()) img(x,y) = NoDataValue();
                }
			}

			delete ptrPixels;
			return img;
		}

		//! Write a Cimg to the file
		GeoRasterIO<T>& Write(CImg<T>& img, bbox chunk, bool BadValCheck=false) {
			point p1 = chunk.min_corner();
			point p2 = chunk.max_corner();
			int width = p2.x()-p1.x()+1;
			int height = p2.y()-p1.y()+1;
			// This the right place for this??
			/*if (Gain() != 1.0) {
                cimg_for(img,ptr,T) if (*ptr != NoDataValue()) *ptr = *ptr/Gain() - Offset();
                std::cout << "Writing...gain = " << Gain() << ", " << Offset() << std::endl;
			}*/
			if (BadValCheck) {
				cimg_for(img,ptr,T) if ( std::isinf(*ptr) || std::isnan(*ptr) ) *ptr = NoDataValue();
			}
			CPLErr err = _GDALRasterBand->RasterIO(GF_Write, p1.x(), p1.y(), width, height, img.data(), width, height, GDALType(), 0, 0);
			if (err != CE_None) std::cout << "Error writing " << Filename() << ": " << CPLGetLastErrorMsg() << std::endl;
			return *this;
		}

		/*GeoRasterIO<T>& Copy(const GeoRasterIO<T>& src) {
            return Process(src,false);
		}*/

		//! Copy input band into this
		GeoRasterIO<T>& Copy(const GeoRaster& img, UNITS units=RAW) {
		    GeoRasterIO<float> src(img);
            std::vector<bbox> Chunks = Chunk();
            std::vector<bbox>::const_iterator iChunk;
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                    CImg<float> img = src.Read(*iChunk, units);
                    CImg<T> imgout;
                    //CImg<unsigned char> mask(img.thresh)
                    if (Gain() != 1.0 || Offset() != 0.0) {
                        imgout = (img-Offset()) / Gain();
                    } else imgout.assign(img);
                    if (src.NoDataValue() != NoDataValue()) {
                        cimg_forXY(img,x,y) { if (img(x,y) == src.NoDataValue()) imgout(x,y) = NoDataValue(); }
                    }
                    Write(imgout,*iChunk);
            }
            // Copy relevant metadata
            GDALRasterBand* band = src.GetGDALRasterBand();
            if (img.NoData()) SetNoData(img.NoDataValue());
            CopyCategoryNames(src);
            _GDALRasterBand->SetDescription(band->GetDescription());
            _GDALRasterBand->SetColorInterpretation(band->GetColorInterpretation());
            _GDALRasterBand->SetMetadata(band->GetMetadata());
            CopyCoordinateSystem(src);
            return *this;
		}

        //! Apply mask to (where mask>0 make NoDataValue)
		/*GeoRasterIO& ApplyMask(const GeoRaster& mask) {
            GeoRasterIO<unsigned char> maskio(mask);
            CImg<unsigned char> maskimg;
            CImg<T> img;
            std::vector<bbox> Chunks = Chunk();
            std::vector<bbox>::const_iterator iChunk;
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                img = Read(*iChunk);
                maskimg = maskio.Read(*iChunk);
                cimg_forXY(img,x,y) if (maskimg(x,y) == 0) img(x,y) = NoDataValue();
                Write(img,*iChunk);
            }
            return *this;
		}*/

        //! Get Saturation mask
		CImg<bool> SaturationMask(bbox chunk) const {
		    CImg<float> band(Read(chunk, RAW));
		    return band.threshold(_maxDC);
		}

		//! NoData mask (all bands)
		CImg<unsigned char> NoDataMask(bbox chunk) const {
		    CImg<T> img = Read(chunk);
		    CImg<unsigned char> mask(img.width(),img.height(),1,1,0);
		    if (!NoData()) return mask;
		    T nodataval = NoDataValue();
            cimg_forXY(img,x,y) if (img(x,y) == nodataval) mask(x,y) = 1;
            return mask;
		}

		//! Creates map of indices for each value in image (used for rasterized vector images)
		/*map<T,POI> MapPOIs(GeoMask<T> Mask) const {
			map<T,POI> POIs;
			vector< box<point> > Chunks = _GeoImage->Chunk();
			vector< box<point> >::const_iterator iChunk;
			CImg<T> img, mask;
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
					img = Read(*iChunk);
					mask = Mask.Read(*iChunk);
					// Apply Nodata mask
					//if (Mask.Valid())
					img.mul(mask);
					cimg_forXY(img,x,y) {
							if (img(x,y) != 0) {
									// If first pixel with this ID, add new
									if (POIs.find(img(x,y)) != POIs.end()) {
											POIs.insert( pair<T,POI>(img(x,y),POI() ) );
											//cout << "Inserting " << img(x,y) << " " << x << " " << y << endl;
									}
									// Add pixel to this POI collection
									POIs[img(x,y)].AddPixel( point_2d(iChunk->min_corner().x()+x,iChunk->min_corner().y()+y) );
							}
					}
			}
			typename map<T,POI>::iterator iPOI;
			for (iPOI=POIs.begin();iPOI!=POIs.end();iPOI++) {
					iPOI->second.Extent();
			}
			return POIs;
		}*/

		//! Copy raster band
		/*GeoRasterIO<T>& Copy(const GeoRaster& src) {
			return Copy(GeoRasterIO<T>(src));
		}*/

		// \name Processing
		// Normalized difference algorithm (NDVI, NDWI, etc)
		/*GeoRasterIO& NormDiff(const GeoRaster& band1, const GeoRaster& band2, std::string desc) {
			// This file needs to use a NoDataValue - use band1 val if absent
			//double val( NoDataValue() );
			//if (!NoData())
			//	SetNoData( (band1.NoData()) ? band1.NoDataValue() : DefaultNoDataValue() );
			band1.NoData() ? SetNoData(band1.NoDataValue()): SetNoData();
			std::vector<bbox> Chunks = Chunk();
			std::vector<bbox>::const_iterator iChunk;
			GeoRasterIO<T> b1proc(band1);
			GeoRasterIO<T> b2proc(band2);
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> b1 = b1proc.Read(*iChunk);
				CImg<T> b2 = b2proc.Read(*iChunk);
				CImg<T> imgout = (b2-b1).div(b2+b1);

				// Check for NoData
				if (band1.NoData() || band2.NoData()) {
					CImg<T> mask = b1proc.NoDataMask(b1) & b2proc.NoDataMask(b2);
					cimg_forXY(imgout,x,y) if (!mask(x,y)) imgout(x,y) = NoDataValue();
				}
				Write(imgout,*iChunk, true);
			}
			if (desc != "") SetDescription(desc.c_str());
			return *this;
		}*/

		// EVI Algorithm
		/*GeoRasterIO& EVI(const GeoRaster& blue, const GeoRaster& red, const GeoRaster& nir) {
			GeoRasterIO<T> Pblue(blue);
			GeoRasterIO<T> Pred(red);
			GeoRasterIO<T> Pnir(nir);
			float G = 2.5;
			float C1 = 6;
			float C2 = 7.5;
			float L = 1;
			// If blue NoData set then use same value, otherwise set to a default
			blue.NoData() ? SetNoData(blue.NoDataValue()) : SetNoData();
			// Chunk and process
			std::vector<bbox> Chunks = Chunk();
			std::vector<bbox>::const_iterator iChunk;
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> Cblue = Pblue.Read(*iChunk);
				CImg<T> Cred = Pred.Read(*iChunk);
				CImg<T> Cnir = Pnir.Read(*iChunk);
				CImg<T> imgout = G*(Cnir-Cred).div(Cnir + C1*Cred - C2*Cblue + L);
				// Check for NoData
				if (blue.NoData() || red.NoData() || nir.NoData()) {
					CImg<T> mask = Pblue.NoDataMask(Cblue) & Pred.NoDataMask(Cred) & Pnir.NoDataMask(Cnir);
					cimg_forXY(imgout,x,y) if (!mask(x,y)) imgout(x,y) = NoDataValue();
				}
				Write(imgout,*iChunk, true);
			}
			SetDescription("EVI");
			return *this;
		}*/

		// SATVI Algorithm
		/*GeoRasterIO& SATVI(const GeoRaster& red, const GeoRaster& swir1, const GeoRaster& swir2) {
			GeoRasterIO<T> Pred(red);
			GeoRasterIO<T> Pswir1(swir1);
			GeoRasterIO<T> Pswir2(swir2);
			float L = 0.1;
			red.NoData() ? SetNoData(red.NoDataValue()) : SetNoData();
			std::vector<bbox> Chunks = Chunk();
			std::vector<bbox>::const_iterator iChunk;
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> Cred = Pred.Read(*iChunk);
				CImg<T> Cswir1 = Pswir1.Read(*iChunk);
				CImg<T> Cswir2 = Pswir2.Read(*iChunk);
				CImg<T> imgout = (((1.0+L)*(Cswir1 - Cred)).div(Cswir1+Cred+L)) - (0.5*Cswir2);

				// Check for NoData
				if (red.NoData() || swir1.NoData() || swir2.NoData()) {
					CImg<T> mask = Pred.NoDataMask(Cred) & Pswir1.NoDataMask(Cswir1) & Pswir2.NoDataMask(Cswir2);
					cimg_forXY(imgout,x,y) if (!mask(x,y)) imgout(x,y) = NoDataValue();
				}
				Write(imgout,*iChunk);
			}
			SetDescription("SATVI");
			return *this;
		}*/

		// Return mask of NoData/NaN/Inf values
		/*CImg<bool> NoDataMask(const CImg<T>& img) {
			CImg<bool> imgout(img.width(),img.height(),img.depth(),img.spectrum(),1);
			if (NoData()) {
				double val = NoDataValue();
				cimg_forXY(img,x,y)
					if (img(x,y) == val || std::isinf(img(x,y)) || std::isnan(img(x,y))) imgout(x,y) = 0;
			}
			return imgout;
		}*/

	private:
		// Private default constructor prevents direct creation
		GeoRasterIO() {}

		//! Returns GDAL Type corresponding to template type T
		GDALDataType GDALType() const {
			if (typeid(T) == typeid(unsigned char)) return GDT_Byte;
			else if (typeid(T) == typeid(unsigned short)) return GDT_UInt16;
			else if (typeid(T) == typeid(short)) return GDT_Int16;
			else if (typeid(T) == typeid(unsigned int)) return GDT_UInt32;
			else if (typeid(T) == typeid(int)) return GDT_Int32;
			else if (typeid(T) == typeid(float)) return GDT_Float32;
			else if (typeid(T) == typeid(double)) return GDT_Float64;
			else {
				std::cout << "Data Type " << typeid(T).name() << " not supported" << std::endl;
				throw(std::exception());
			}
		}

	}; //class GeoRasterIO
} // namespace gip

#endif
