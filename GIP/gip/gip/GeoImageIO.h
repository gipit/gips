/*
 * GeoImageIO.h
 *
 *  Created on: Dec 19, 2011
 *      Author: mhanson
 */

#ifndef GEOIMAGEIO_H_
#define GEOIMAGEIO_H_

#include <gip/GeoImage.h>
#include <gip/GeoRasterIO.h>

namespace gip {

	template<class T> class GeoImageIO : public GeoImage {
	public:
        typedef boost::geometry::model::box<point> bbox;
		//! \name Constructors/Destructor
		GeoImageIO() {}
		GeoImageIO(GeoImage& img)
			: GeoImage(img) {
			LoadRasterIO();
		}
		GeoImageIO(const GeoImage& img)
			: GeoImage(img) {
			LoadRasterIO();
		}
		~GeoImageIO() {}

        // TODO - remove this?
		void LoadRasterIO() {
			for (unsigned int i=0;i<_RasterBands.size();i++) {
				_RasterIOBands.push_back(GeoImage::_RasterBands[i]);
			}
		}

		//! \name Band Operations
		//! Get raster band (0-based index)
		//GeoRasterIO<T> operator[](int band) { return GeoRasterIO<T>(GeoImage::operator[](band)); }
		//! Get raster band, const version
		//const GeoRasterIO<T> operator[](int band) const { return GeoRasterIO<T>(GeoImage::operator[](band)); }
		//! Get raster band by color
		//GeoRasterIO<T> operator[](std::string col) { return GeoRasterIO<T>(GeoImage::operator[](col)); }
		//! Get raster band by color, const version
		//const GeoRasterIO<T> operator[](std::string col) const { return GeoRasterIO<T>(GeoImage::operator[](col)); }

        //! Get raster band (0-based index)
		GeoRasterIO<T>& operator[](int band) { return _RasterIOBands[band]; }
		//! Get raster band, const version
		const GeoRasterIO<T>& operator[](int band) const { return _RasterIOBands[band]; }
		//! Get raster band by color
		GeoRasterIO<T>& operator[](std::string col) {
			// Call const version
			return const_cast<GeoRasterIO<T>&>(static_cast<const GeoImageIO<T>&>(*this)[col]);
        }
		//! Get raster band by color, const version
		const GeoRasterIO<T>& operator[](std::string col) const {
            int index(_Colors[col]);
            if (index > 0)
                return _RasterIOBands[index-1];
            else {
                // TODO - fix this - can't return NULL reference...?
                //std::cout << "No band of that color, returning band 0" << std::endl;
                throw std::out_of_range ("No band of color "+col);
            }
        }

		//! \name File I/O
		//! Read chunk, across all bands
		cimg_library::CImg<T> Read(int chunk=0, bool RAW=false) const {
			cimg_library::CImgList<T> images;
			typename std::vector< GeoRasterIO<T> >::const_iterator iBand;
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				images.insert( iBand->Read(chunk, RAW) );
			}
			//return images.get_append('c','p');
			return images.get_append('v','p');
		}

        //GeoImageIO& Write(const CImg<T> img, int chunk=0, bool BadValCheck=false) {
        //    return Write(img, chunk, BadValCheck);
        //}

		//! Write cube across all bands
		GeoImageIO& Write(const CImg<T> img, int chunk=0, bool BadValCheck=false) {
			typename std::vector< GeoRasterIO<T> >::iterator iBand;
			int i(0);
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				CImg<T> tmp = img.get_channel(i++);
				iBand->Write(tmp, chunk, BadValCheck);
			}
			return *this;
		}
		// Read Cube as list
		/*cimg_library::CImgList<T> ReadAsList(bbox chunk) const {
			cimg_library::CImgList<T> images;
			typename std::vector< GeoRasterIO<T> >::const_iterator iBand;
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				images.insert( iBand->Read(chunk) );
			}
			return images;
		}*/

		//! NoData mask (all bands)
		CImg<unsigned char> NoDataMask(int chunk=0, std::vector<std::string> bands=std::vector<std::string>()) const {
		    unsigned int c;
		    CImg<T> cube = Read(chunk, true);
		    CImg<unsigned char> mask(cube.width(),cube.height(),1,1,1);
		    std::vector<int> ibands;
		    std::vector<int>::const_iterator b;
		    if (bands.size() == 0) {
		        for (c=0; c<NumBands(); c++) ibands.push_back(c);
            } else {
                for (std::vector<std::string>::const_iterator i=bands.begin(); i!=bands.end(); i++) {
                    ibands.push_back(_Colors[*i]-1);
                }
            }
            /*if (Options::Verbose() > 2) {
                std::cout << "Retrieving nodata mask for bands ";
                for (b=ibands.begin();b!=ibands.end();b++) std::cout << *b << " ";
                std::cout << std::endl;
            }*/
            cimg_forXY(cube,x,y) {
                for (std::vector<int>::const_iterator i=ibands.begin(); i!=ibands.end(); i++) {
                    if ( (*this)[*i].NoData() && (cube(x,y,*i) == (*this)[*i].NoDataValue())) mask(x,y) = 0;
                }
            }
            return mask;
		}

		// TODO - examine that these are even generically needed, and aren't just part of algorithms (all or mostly fmask)

		//! Return a mask of snow
		CImg<bool> SnowMask(int chunk=0) const {
            CImg<float> nir( operator[]("NIR").Ref(chunk) );
            CImg<float> green( operator[]("Green").Ref(chunk) );
            CImg<float> swir1( (*this)["SWIR1"].Ref(chunk) );
            CImg<float> temp( operator[]("LWIR").Ref(chunk) );

            float th_nir = 0.11;
            float th_green = 0.1;
            float th_ndsi = 0.15;
            float th_temp = 3.80;

            CImg<float> ndsi( (green-swir1).div(green+swir1) );
            CImg<bool> mask = ndsi.threshold(th_ndsi) & nir.threshold(th_nir)
                & green.threshold(th_green) & temp.threshold(th_temp,false,true)^=1;
            return mask;
		}

        //! Return a mask of water (and possibley clear-sky pixels)
		CImg<bool> WaterMask(int chunk=0) const {
            CImg<float> red( operator[]("Red").Ref(chunk) );
            CImg<float> nir( operator[]("NIR").Ref(chunk) );
            CImg<float> ndvi = (nir-red).get_div(nir+red);
            CImg<bool> mask(red.width(),red.height(),1,1,false);
            cimg_forXY(mask,x,y) {
                if ( ((ndvi(x,y) < 0.01) && (nir(x,y) < 0.11)) || ((ndvi(x,y) < 0.1) && (nir(x,y) < 0.05)) ) mask(x,y) = true;
            }
            return mask;
		}

        //! Return haze mask
        CImg<bool> HazeMask(int chunk=0) const {
            CImg<float> red( operator[]("Red").Ref(chunk) );
            CImg<float> blue( operator[]("Blue").Ref(chunk) );
            CImg<bool> mask( (blue - 0.5*red - 0.08).threshold(0.0) );
            return mask;
        }

        CImg<float> Whiteness(int chunk=0) const {
            // RAW or RADIANCE ?
            CImg<float> red = operator[]("Red").Read(chunk, true);
            CImg<float> green = operator[]("Green").Read(chunk, true);
            CImg<float> blue = operator[]("Blue").Read(chunk, true);
            CImg<float> white(red.width(),red.height());
            float mu;
            cimg_forXY(white,x,y) {
                mu = (red(x,y) + green(x,y) + blue(x,y))/3;
                white(x,y) = (abs(red(x,y)-mu) + abs(green(x,y)-mu) + abs(blue(x,y)-mu))/mu;
            }
            // Saturation?  If pixel saturated make Whiteness 0 ?
            return white;
        }

        // NDVI
		/*CImg<float> NDVI(bbox chunk) const {
            CImg<float> red( (*this)["Red"].Ref(chunk) );
            CImg<float> nir( (*this)["NIR"].Ref(chunk) );
            CImg<float> ndvi( (nir-red).div(nir+red) );
            // Check for NoData
            //if (*this["Red"].NoData() || *this["NIR"].NoData()) {
            //    double val = *this["Red"].NoDataValue();
            //    cimg_forXY(ndvi,x,y)
            //        if (red(x,y) == val || nir(x,y) == val) ndvi(x,y) = val;
            //}
            return ndvi;
		}*/

        // Get NDSI
		/*CImg<float> NDSI(bbox chunk) const {
            CImg<float> green( (*this)["Green"].Ref(chunk) );
            CImg<float> swir1( (*this)["SWIR1"].Ref(chunk) );
            CImg<float> ndsi( (green-swir1).div(green+swir1) );
            // Check for NoData
            //if (*this["Green"].NoData() || *this["SWIR1"].NoData()) {
            //    double val = *this["Green"].NoDataValue();
            //    cimg_forXY(ndsi,x,y)
            //        if (green(x,y) == val || swir1(x,y) == val) ndsi(x,y) = val;
            //}
            return ndsi;
		}*/

        //! Extract, and interpolate, time series (C is time axis)
        cimg_library::CImg<T> TimeSeries(cimg_library::CImg<double> C) {
            cimg_library::CImg<T> cimg = Read();
            T nodata = _RasterIOBands[0].NoDataValue();
            if (cimg.spectrum() > 2) {
                int lowi, highi;
                float y0, y1, x0, x1;
                for (int c=1; c<cimg.spectrum()-1;c++) {
                    if (Options::Verbose() > 3) cimg_print(C, "days vector");
                    cimg_forXY(cimg,x,y) {
                        if (cimg(x,y,c) == nodata) {
                            // Find next lowest point
                            lowi = highi = 1;
                            while ((cimg(x,y,c-lowi) == nodata) && (lowi<c)) lowi++;
                            while ((cimg(x,y,c+highi) == nodata) && (c+highi < cimg.spectrum()-1) ) highi++;
                            y0 = cimg(x,y,c-lowi);
                            y1 = cimg(x,y,c+highi);
                            x0 = C(c-lowi);
                            x1 = C(c+highi);
                            if ((y0 != nodata) && (y1 != nodata)) {
                                cimg(x,y,c) = y0 + (y1-y0) * ((C(c)-x0)/(x1-x0));
                            }
                        } else if (cimg(x,y,c-1) == nodata) {
                            T val = cimg(x,y,c);
                            for (int i=c-1; i>=0; i--) {
                                if (cimg(x,y,i) == nodata) cimg(x,y,i) = val;
                            }
                        }
                    }
                }
            }
            return cimg;
        }

        //! Extract spectra from select pixels (where mask > 0)
        cimg_library::CImg<T> Extract(const GeoRaster& mask) {
            if (Options::Verbose() > 2 ) std::cout << "Pixel spectral extraction" << std::endl;
            GeoRasterIO<unsigned char> maskio(mask);
            cimg_library::CImg<unsigned char> cmask;
            cimg_library::CImg<T> cimg;
            long count = 0;

            for (int iChunk=1; iChunk<=NumChunks(); iChunk++) {
                cmask = maskio.Read(iChunk);
                cimg_for(cmask,ptr,unsigned char) if (*ptr > 0) count++;
            }
            cimg_library::CImg<T> pixels(count,NumBands()+1,1,1,_RasterIOBands[0].NoDataValue());
            count = 0;
            int ch(0);
            unsigned int c;
            for (int iChunk=1; iChunk<=NumChunks(); iChunk++) {
                cimg = Read(iChunk);
                cmask = maskio.Read(iChunk);
                cimg_forXY(cimg,x,y) {
                    if (cmask(x,y) > 0) {
                        for (c=0;c<NumBands();c++) pixels(count,c+1) = cimg(x,y,c);
                        pixels(count++,0) = cmask(x,y);
                    }
                }
                if (Options::Verbose() > 2) std::cout << "  Chunk " << ch++ << std::endl;
            }
            return pixels;
        }

        //! Get a number of random pixel vectors (spectral vectors)
        CImg<T> GetRandomPixels(int NumPixels) const {
            CImg<T> Pixels(NumBands(), NumPixels);
            srand( time(NULL) );
            bool badpix;
            int p = 0;
            while(p < NumPixels) {
                int col = (double)rand()/RAND_MAX * (XSize()-1);
                int row = (double)rand()/RAND_MAX * (YSize()-1);
                T pix[1];
                badpix = false;
                for (unsigned int j=0; j<NumBands(); j++) {
                    _RasterIOBands[j].GetGDALRasterBand()->RasterIO(GF_Read, col, row, 1, 1, &pix, 1, 1, Type(), 0, 0);
                    if (_RasterIOBands[j].NoData() && pix[0] == _RasterIOBands[j].NoDataValue()) {
                        badpix = true;
                    } else {
                        Pixels(j,p) = pix[0];
                    }
                }
                if (!badpix) p++;
            }
            return Pixels;
        }

        //! Get a number of pixel vectors that are spectrally distant from each other
        CImg<T> GetPixelClasses(int NumClasses) const {
            int RandPixelsPerClass = 500;
            CImg<T> stats;
            CImg<T> ClassMeans(NumBands(), NumClasses);
            // Get Random Pixels
            CImg<T> RandomPixels = GetRandomPixels(NumClasses * RandPixelsPerClass);
            // First pixel becomes first class
            cimg_forX(ClassMeans,x) ClassMeans(x,0) = RandomPixels(x,0);
            for (int i=1; i<NumClasses; i++) {
                CImg<T> ThisClass = ClassMeans.get_row(i-1);
                long validpixels = 0;
                CImg<T> Dist(RandomPixels.height());
                for (long j=0; j<RandomPixels.height(); j++) {
                    // Get current pixel vector
                    CImg<T> ThisPixel = RandomPixels.get_row(j);
                    // Find distance to last class
                    Dist(j) = ThisPixel.sum() ? (ThisPixel-ThisClass).dot( (ThisPixel-ThisClass).transpose() ) : 0;
                    if (Dist(j) != 0) validpixels++;
                }
                stats = Dist.get_stats();
                // The pixel farthest away from last class make the new class
                cimg_forX(ClassMeans,x) ClassMeans(x,i) = RandomPixels(x,stats(8));
                // Toss a bunch of pixels away (make zero)
                CImg<T> DistSort = Dist.get_sort();
                T cutoff = DistSort[RandPixelsPerClass*i]; //(stats.max-stats.min)/10 + stats.min;
                cimg_forX(Dist,x) if (Dist(x) < cutoff) cimg_forX(RandomPixels,x1) RandomPixels(x1,x) = 0;
            }
            // Output Class Vectors
            //if (Options::Verbose()>1) cimg_printclasses(ClassMeans, "Initial Class");
            return ClassMeans;
        }

	protected:
		//! Vector of raster bands
		std::vector< GeoRasterIO<T> > _RasterIOBands;

	private:


		//! Returns GDAL Type corresponding to template type T
		GDALDataType Type() const {
		    // TODO - this is repeated in GeoRaster
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

	}; // class GeoImageIO
} // namespace gip

#endif /* GEOIMAGEIO_H_ */
