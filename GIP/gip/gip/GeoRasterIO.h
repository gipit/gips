#ifndef GIP_GEORASTERIO_H
#define GIP_GEORASTERIO_H

#include <exception>
#include <gip/GeoRaster.h>
#include <cmath>

// only used for debugging
#include <iostream>

namespace gip {

    // TODO - Compile with C++0x to use scoped enums
    //enum UNITS {RAW, RADIANCE, REFLECTIVITY};
    typedef boost::geometry::model::box<point> bbox;

	template<class T> class GeoRasterIO : public GeoRaster {
	public:
        typedef boost::geometry::model::box<point> bbox;
        //! \name Constructors/Destructor
		GeoRasterIO(GeoRaster& img)
			: GeoRaster(img) {}
		GeoRasterIO(const GeoRaster& img)
			: GeoRaster(img) {}
		~GeoRasterIO() {}

		//! \name File I/O
		//! Retrieve a piece of the image as a CImg
		cimg_library::CImg<T> Read(int chunknum=0, bool RAW=false) const {
		    bbox chunk;
		    if (chunknum == 0) {
                chunk = bbox(point(0,0),point(XSize()-1,YSize()-1));
		    } else {
                chunk = _Chunks[chunknum-1];
		    }
		    using cimg_library::CImg;
			point p1 = chunk.min_corner();
			point p2 = chunk.max_corner();
			int width = p2.x()-p1.x()+1;
			int height = p2.y()-p1.y()+1;
			//std::cout << Basename() << " reading " << boost::geometry::dsv(p1) << boost::geometry::dsv(p2) << " w=" << width << " h=" << height << std::endl;
			T* ptrPixels = new T[width*height];
			// casting away const, safe because this is a read-only const_cast<GDALRasterBand*>
			CPLErr err = _GDALRasterBand->RasterIO(GF_Read, p1.x(), p1.y(), width, height, ptrPixels, width, height, this->Type(), 0, 0);
			if (err != CE_None) std::cout << "Error reading " << Filename() << ": " << CPLGetLastErrorMsg() << std::endl;
			CImg<T> img(ptrPixels, width, height);
			CImg<T> imgorig(img);

			bool updatenodata = false;
			// Convert data to radiance (if not raw requested)
			if (!RAW) {
                if (Gain() != 1.0 || Offset() != 0.0) {
                    img = Gain() * (img-_minDC) + Offset();
                    updatenodata = true;
                }
                // apply atmosphere if there is one (which would data is radiance units) TODO - check units
                if (Atmosphere()) {
                    if (Options::Verbose() > 1) std::cout << Basename() << ": applying atmosphere" << std::endl;
                    double e = (Thermal()) ? 0.95 : 1;  // For thermal band, currently water only
                    img = (img - (_Atmosphere.Lu() + (1-e)*_Atmosphere.Ld())) / (_Atmosphere.t() * e);
                    updatenodata = true;
                }
			}

			// Apply Processing functions - TODO - should these apply if RAW ?
			std::vector<GeoFunction>::const_iterator iFunc;
			for (iFunc=_Functions.begin();iFunc!=_Functions.end();iFunc++) {
			    if (Options::Verbose() > 1)
                    std::cout << Basename() << ": Applying function " << iFunc->Function() << " " << iFunc->Operand() << std::endl;
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
				}
				updatenodata = true;
			}
			// Apply all masks
			if (_Masks.size() > 0) {
			    if (Options::Verbose() > 1)
                    std::cout << Basename() << ": Applying " << _Masks.size() << " masks" << std::endl;
                GeoRasterIO<unsigned char> mask(_Masks[0]);
                CImg<unsigned char> cmask(mask.Read(chunknum));
                for (unsigned int i=1; i<_Masks.size(); i++) {
                    mask = GeoRasterIO<unsigned char>(_Masks[i]);
                    cmask.mul(mask.Read(chunknum));
                }
                img.mul(cmask);
                cimg_forXY(img,x,y) {
                    if (cmask(x,y) != 1) img(x,y) = NoDataValue();
                }
                updatenodata = true;
			}
			// If processing was applied update NoData values where needed
			if (NoData() && updatenodata) {
                cimg_forXY(img,x,y) {
                    if (imgorig(x,y) == NoDataValue()) img(x,y) = NoDataValue();
                }
			}
			delete ptrPixels;
			return img;
		}

		//GeoRasterIO<T>& Write(cimg_library::CImg<T> img, int chunknum=0, bool RAW=false) {
        //    return Write(img, chunknum, RAW);
		//}

		//! Write a Cimg to the file
		GeoRasterIO<T>& Write(cimg_library::CImg<T> img, int chunknum=0, bool RAW=false) { //, bool BadValCheck=false) {
		    bbox chunk;
		    if (chunknum == 0) {
                chunk = bbox(point(0,0),point(XSize()-1,YSize()-1));
		    } else {
                chunk = _Chunks[chunknum-1];
		    }
			point p1 = chunk.min_corner();
			point p2 = chunk.max_corner();
			int width = p2.x()-p1.x()+1;
			int height = p2.y()-p1.y()+1;
			// This the right place for this??
			if (!RAW && (Gain() != 1.0 || Offset() != 0.0)) {
                cimg_for(img,ptr,T) if (*ptr != NoDataValue()) *ptr = (*ptr-Offset())/Gain();
			}
            if (Options::Verbose() > 2)
                std::cout << Basename() << ": Writing with gain " << Gain() << " and offset " << Offset() << std::endl;
			/*if (BadValCheck) {
				cimg_for(img,ptr,T) if ( std::isinf(*ptr) || std::isnan(*ptr) ) *ptr = NoDataValue();
			}*/
			CPLErr err = _GDALRasterBand->RasterIO(GF_Write, p1.x(), p1.y(), width, height, img.data(), width, height, this->Type(), 0, 0);
			if (err != CE_None) std::cout << "Error writing " << Filename() << ": " << CPLGetLastErrorMsg() << std::endl;
			return *this;
		}

		//! Process input band into this
		GeoRasterIO<T>& Process(const GeoRaster& raster, bool RAW=false) {
		    using cimg_library::CImg;
		    GeoRasterIO<double> rasterIO(raster);
            for (int iChunk=1; iChunk<=NumChunks(); iChunk++) {
                    CImg<double> cimg = rasterIO.Read(iChunk, RAW);
                    //CImg<unsigned char> mask;
                    //if (Gain() != 1.0 || Offset() != 0.0) {
                    //    (cimg-=Offset())/=Gain();
                    //    mask = rasterIO.NoDataMask(*iChunk);
                    //    cimg_forXY(cimg,x,y) { if (mask(x,y)) cimg(x,y) = NoDataValue(); }
                    //}
                    //WriteChunk(CImg<T>().assign(cimg.round()),*iChunk, RAW);
                    Write(CImg<T>().assign(cimg),iChunk, RAW);
            }
            // Copy relevant metadata
            GDALRasterBand* band = raster.GetGDALRasterBand();
            //if (img.NoData()) SetNoData(img.NoDataValue());
            CopyCategoryNames(raster);
            _GDALRasterBand->SetDescription(band->GetDescription());
            _GDALRasterBand->SetColorInterpretation(band->GetColorInterpretation());
            _GDALRasterBand->SetMetadata(band->GetMetadata());
            CopyCoordinateSystem(raster);
            return *this;
		}

        //! Get Saturation mask
		cimg_library::CImg<bool> SaturationMask(int chunk=0) const {
		    cimg_library::CImg<float> band(Read(chunk, true));
		    return band.threshold(_maxDC);
		}

		//! NoData mask
		cimg_library::CImg<unsigned char> NoDataMask(int chunk=0) const {
		    using cimg_library::CImg;
		    CImg<T> img = Read(chunk, true);  // this reads raw
		    CImg<unsigned char> mask(img.width(),img.height(),1,1,0);
		    if (!NoData()) return mask;
		    T nodataval = NoDataValue();
            cimg_forXY(img,x,y) if (img(x,y) == nodataval) mask(x,y) = 1;
            return mask;
		}

		//! Return reflectance (or temperature if thermal band)  (move to GeoRaster)
		cimg_library::CImg<float> Ref(int chunk=0) const {
            cimg_library::CImg<float> cimg = Read(chunk);
            if (Units() == "reflectance") {
                return cimg;
            } else if (Units() != "radiance") {
                throw std::runtime_error("image not in compatible units for reflectance");
            }
            if (Thermal()) {
                cimg_for(cimg,ptr,float) *ptr = (_K2/log(_K1/(*ptr)+1)) - 273.15;
            } else {
                float normrad = Atmosphere() ? (1.0/_Atmosphere.Ld()) : (1.0/_Esun);
                cimg_for(cimg,ptr,float) *ptr = *ptr * normrad;
            }
            return cimg;
		}

		//! Calculate stats
		cimg_library::CImg<float> ComputeStats(bool RAW=false) {
		    using cimg_library::CImg;
		    CImg<T> cimg;
		    T min(MaxValue()), max(MinValue());
		    long count(0);
		    double total(0);
            for (int iChunk=1; iChunk<=NumChunks(); iChunk++) {
                cimg = Read(iChunk, RAW);
                cimg_for(cimg,ptr,T) {
                    if (*ptr != NoDataValue()) {
                        total += *ptr;
                        count++;
                        if (*ptr > max) max = *ptr;
                        if (*ptr < min) min = *ptr;
                    }
                }
            }
            float mean = total/count;
            total = 0;
            for (int iChunk=1; iChunk<=NumChunks(); iChunk++) {
                cimg = Read(iChunk, RAW);
                cimg_for(cimg,ptr,T) {
                    if (*ptr != NoDataValue()) total += (*ptr - mean)*(*ptr - mean);
                }
            }
            float stdev = sqrt(total/count);
            CImg<float> stats(4,1,1,1,(float)min,(float)max,mean,stdev);
            return stats;
		}

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

		// Creates map of indices for each value in image (used for rasterized vector images)
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

	private:
		// Private default constructor prevents direct creation
		GeoRasterIO() {}

		//! Returns GDAL Type corresponding to template type T
		GDALDataType Type() const {
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
