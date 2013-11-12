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
		// Constructors/Destructor
		GeoImageIO(GeoImage& img)
			: GeoImage(img) {
			LoadRasterIO();
		}
		GeoImageIO(const GeoImage& img)
			: GeoImage(img) {
			LoadRasterIO();
		}
        GeoImageIO(std::string filename)
            : GeoImage(filename, true) {
            LoadRasterIO();
        }
		~GeoImageIO() {}

		// \name Band Operations
		//! Get raster band (0-based index)
		GeoRasterIO<T> operator[](int band) { return GeoRasterIO<T>(GeoImage::operator[](band)); }
		//! Get raster band, const version
		const GeoRasterIO<T> operator[](int band) const { return GeoRasterIO<T>(GeoImage::operator[](band)); }
		//! Get raster band by color
		GeoRasterIO<T> operator[](std::string col) { return GeoRasterIO<T>(GeoImage::operator[](col)); }
		//! Get raster band by color, const version
		const GeoRasterIO<T> operator[](std::string col) const { return GeoRasterIO<T>(GeoImage::operator[](col)); }

        //! Extract spectra from select pixels (where mask > 0)
        cimg_library::CImg<T> Extract(const GeoRaster& mask) {
            if (Options::Verbose() > 2 ) std::cout << "Pixel spectral extraction" << std::endl;
            GeoRasterIO<unsigned char> maskio(mask);
            cimg_library::CImg<unsigned char> cmask;
            cimg_library::CImg<T> cimg;
            long count = 0;
            float chunksize = Options::ChunkSize();
            Options::SetChunkSize(chunksize/NumBands());
            std::vector<bbox> Chunks = Chunk();
            std::vector<bbox>::const_iterator iChunk;
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                cmask = maskio.Read(*iChunk);
                cimg_for(cmask,ptr,unsigned char) if (*ptr > 0) count++;
            }
            cimg_library::CImg<T> pixels(count,NumBands()+1,1,1,_RasterIOBands[0].NoDataValue());
            count = 0;
            int ch(0);
            unsigned int c;
            for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                cimg = Read(*iChunk);
                cmask = maskio.Read(*iChunk);
                cimg_forXY(cimg,x,y) {
                    if (cmask(x,y) > 0) {
                        for (c=0;c<NumBands();c++) pixels(count,c+1) = cimg(x,y,c);
                        pixels(count++,0) = cmask(x,y);
                    }
                }
                if (Options::Verbose() > 2) std::cout << "  Chunk " << ch++ << std::endl;
            }
            Options::SetChunkSize(chunksize);
            return pixels;
        }

		//! Get raster band by color
		/*GeoRasterIO<T>& operator[](Color::Enum col) {
			// Call const version
			return const_cast<GeoRasterIO<T>&>(static_cast<const GeoImageIO&>(*this)[col]);
		}
		//! Get raster band by color, const version
		const GeoRasterIO<T>& operator[](Color::Enum col) const {
			//std::map<Color::Enum,int>::const_iterator pos;
			//pos = _Colors.find(col);
			int index(_Colors.Band(col));
			if (index > 0)
				return _RasterIOBands[index-1];
			else {
				// TODO - fix this - can't return NULL reference...?
				std::cout << "No band of that color, returning band 0" << std::endl;
				return _RasterIOBands[0];
			}
		}*/

		void LoadRasterIO() {
			for (unsigned int i=0;i<_RasterBands.size();i++) {
				_RasterIOBands.push_back(GeoImage::_RasterBands[i]);
			}
		}

		//! Read Cube
		CImg<T> Read(bbox chunk, bool RAW=false) const {
			cimg_library::CImgList<T> images;
			typename std::vector< GeoRasterIO<T> >::const_iterator iBand;
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				images.insert( iBand->Read(chunk, RAW) );
			}
			//return images.get_append('c','p');
			return images.get_append('v','p');
		}

		//! Read Cube as list
		cimg_library::CImgList<T> ReadAsList(bbox chunk) const {
			cimg_library::CImgList<T> images;
			typename std::vector< GeoRasterIO<T> >::const_iterator iBand;
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				images.insert( iBand->Read(chunk) );
			}
			return images;
		}

		//! Write Cube of data (for all bands)
		GeoImageIO& Write(const CImg<T> &img, bbox chunk, bool BadValCheck=false) {
			typename std::vector< GeoRasterIO<T> >::iterator iBand;
			int i(0);
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				CImg<T> tmp = img.get_channel(i++);
				iBand->WriteChunk(tmp,chunk, BadValCheck);
			}
			return *this;
		}

        //! Get NDVI
		CImg<float> NDVI(bbox chunk) const {
            CImg<float> red( (*this)["Red"].Ref(chunk) );
            CImg<float> nir( (*this)["NIR"].Ref(chunk) );
            CImg<float> ndvi( (nir-red).div(nir+red) );

            // Check for NoData
            /*if (*this["Red"].NoData() || *this["NIR"].NoData()) {
                double val = *this["Red"].NoDataValue();
                cimg_forXY(ndvi,x,y)
                    if (red(x,y) == val || nir(x,y) == val) ndvi(x,y) = val;
            }*/
            return ndvi;
		}

        //! Get NDSI
		CImg<float> NDSI(bbox chunk) const {
            CImg<float> green( (*this)["Green"].Ref(chunk) );
            CImg<float> swir1( (*this)["SWIR1"].Ref(chunk) );
            CImg<float> ndsi( (green-swir1).div(green+swir1) );

            // Check for NoData
            /*if (*this["Green"].NoData() || *this["SWIR1"].NoData()) {
                double val = *this["Green"].NoDataValue();
                cimg_forXY(ndsi,x,y)
                    if (green(x,y) == val || swir1(x,y) == val) ndsi(x,y) = val;
            }*/
            return ndsi;
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
                    _RasterIOBands[j].GetGDALRasterBand()->RasterIO(GF_Read, col, row, 1, 1, &pix, 1, 1, GDALType(), 0, 0);
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

		// MASKS
		//! NoData mask (all bands)
		CImg<unsigned char> NoDataMask(bbox chunk) const {
		    unsigned int c;
		    CImg<T> cube = Read(chunk, true);
		    CImg<unsigned char> mask(cube.width(),cube.height(),1,1,0);
            cimg_forXY(cube,x,y) {
                for (c=0; c<NumBands(); c++) {
                    if ( (*this)[c].NoData() && (cube(x,y,c) == (*this)[c].NoDataValue())) mask(x,y) = 1;
                }
            }
            return mask;
		}

		//! Return a mask of snow
		CImg<bool> SnowMask(bbox chunk) const {
            CImg<float> nir( operator[]("NIR").Ref(chunk) );
            CImg<float> green( operator[]("Green").Ref(chunk) );
            CImg<float> temp( operator[]("LWIR").Ref(chunk) );

            float th_nir = 0.11;
            float th_green = 0.1;
            float th_ndsi = 0.15;
            float th_temp = 3.80;

            CImg<float> ndsi = NDSI(chunk);
            CImg<bool> mask = ndsi.threshold(th_ndsi) & nir.threshold(th_nir)
                & green.threshold(th_green) & temp.threshold(th_temp,false,true)^=1;
            return mask;
		}

        //! Return a mask of water (and possibley clear-sky pixels)
		CImg<bool> WaterMask(bbox chunk) const {
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
        CImg<bool> HazeMask(bbox chunk) const {
            CImg<float> red( operator[]("Red").Ref(chunk) );
            CImg<float> blue( operator[]("Blue").Ref(chunk) );
            CImg<bool> mask( (blue - 0.5*red - 0.08).threshold(0.0) );
            return mask;
        }

        CImg<float> Whiteness(bbox chunk) const {
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

	protected:
		//! Vector of raster bands
		std::vector< GeoRasterIO<T> > _RasterIOBands;

	private:
		// Private default constructor prevents direct creation
		GeoImageIO() {}

		//! Returns GDAL Type corresponding to template type T
		GDALDataType GDALType() const {
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
