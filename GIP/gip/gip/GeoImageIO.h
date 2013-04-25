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
	using cimg_library::CImgList;

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
		CImg<T> ReadCube(bbox chunk) const {
			CImgList<T> images;
			typename std::vector< GeoRasterIO<T> >::const_iterator iBand;
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				images.insert( iBand->Read(chunk) );
			}
			//return images.get_append('c','p');
			return images.get_append('c');
		}

		//! Read Cube as list
		CImgList<T> ReadCubeAsList(bbox chunk) const {
			CImgList<T> images;
			typename std::vector< GeoRasterIO<T> >::const_iterator iBand;
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				images.insert( iBand->Read(chunk) );
			}
			return images;
		}

		//! Write Cube of data (for all bands)
		GeoImageIO& WriteCube(const CImg<T> &img, bbox chunk, bool BadValCheck=false) {
			typename std::vector< GeoRasterIO<T> >::iterator iBand;
			int i(0);
			for (iBand=_RasterIOBands.begin();iBand!=_RasterIOBands.end();iBand++) {
				CImg<T> tmp = img.get_channel(i++);
				iBand->Write(tmp,chunk, BadValCheck);
			}
			return *this;
		}

        //! Get NDVI
		CImg<float> NDVI(bbox chunk) const {
            CImg<float> red( (*this)["Red"].Read(chunk, REFLECTIVITY) );
            CImg<float> nir( (*this)["NIR"].Read(chunk, REFLECTIVITY) );
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
            CImg<float> green( (*this)["Green"].Read(chunk, REFLECTIVITY) );
            CImg<float> swir1( (*this)["SWIR1"].Read(chunk, REFLECTIVITY) );
            CImg<float> ndsi( (green-swir1).div(green+swir1) );

            // Check for NoData
            /*if (*this["Green"].NoData() || *this["SWIR1"].NoData()) {
                double val = *this["Green"].NoDataValue();
                cimg_forXY(ndsi,x,y)
                    if (green(x,y) == val || swir1(x,y) == val) ndsi(x,y) = val;
            }*/
            return ndsi;
		}

		// MASKS

		//! NoData mask (all bands)
		CImg<bool> NoDataMask(bbox chunk) const {
		    unsigned int c;
		    CImg<T> cube = ReadCube(chunk);
		    CImg<bool> mask(cube.width(),cube.height(),1,1,0);
            cimg_forXY(cube,x,y) {
                for (c=0; c<NumBands(); c++) {
                    if ( (*this)[c].NoData() && (cube(x,y,c) == (*this)[c].NoDataValue())) mask(x,y) = 1;
                }
            }
            return mask;
		}

		//! Return a mask of snow
		CImg<bool> SnowMask(bbox chunk) const {
            CImg<float> nir( *this["NIR"].Read(chunk, REFLECTIVITY) );
            CImg<float> green( *this["Green"].Read(chunk, REFLECTIVITY) );
            CImg<float> temp( *this["LWIR"].Read(chunk, REFLECTIVITY) );

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
            CImg<float> red( (*this)["Red"].Read(chunk, REFLECTIVITY) );
            CImg<float> nir( (*this)["NIR"].Read(chunk, REFLECTIVITY) );
            CImg<float> ndvi = (nir-red).get_div(nir+red);
            CImg<bool> mask(red.width(),red.height(),1,1,false);
            cimg_forXY(mask,x,y) {
                if ( ((ndvi(x,y) < 0.01) && (nir(x,y) < 0.11)) || ((ndvi(x,y) < 0.1) && (nir(x,y) < 0.05)) ) mask(x,y) = true;
            }
            return mask;
		}

        //! Return haze mask
        CImg<bool> HazeMask(bbox chunk) const {
            CImg<float> red( (*this)["Red"].Read(chunk, REFLECTIVITY) );
            CImg<float> blue( (*this)["Blue"].Read(chunk, REFLECTIVITY) );
            CImg<bool> mask( (blue - 0.5*red - 0.08).threshold(0.0) );
            return mask;
        }

        CImg<float> Whiteness(bbox chunk) const {
            // RAW or RADIANCE ?
            CImg<float> red = (*this)["Red"].Read(chunk, RAW);
            CImg<float> green = (*this)["Green"].Read(chunk, RAW);
            CImg<float> blue = (*this)["Blue"].Read(chunk, RAW);
            CImg<float> white(red.width(),red.height());
            float mu;
            cimg_forXY(white,x,y) {
                mu = (red(x,y) + green(x,y) + blue(x,y))/3;
                white(x,y) = (abs(red(x,y)-mu) + abs(green(x,y)-mu) + abs(blue(x,y)-mu))/mu;
            }
            // Saturation?  If pixel saturated make Whiteness 0 ?
            return white;
        }

        //! Copy from input band using specified units
		/*GeoImage& Process(const GeoImage& image) {
		    if (image.NumBands() != NumBands()) throw std::exception();
            for (unsigned int b=0;b<NumBands();b++) {
                _RasterIOBands[b].Process(GeoRasterIO<float>(image[b]));
            }
            return *this;
		}*/

	protected:
		//! Vector of raster bands
		std::vector< GeoRasterIO<T> > _RasterIOBands;

	private:
		// Private default constructor prevents direct creation
		GeoImageIO() {}

	}; // class GeoImageIO
} // namespace gip



#endif /* GEOIMAGEIO_H_ */
