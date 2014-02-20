#ifndef GIP_UTILS_H
#define GIP_UTILS_H

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace bgeo = boost::geometry;

//using namespace cimg_library;

namespace gip {
	typedef bgeo::model::d2::point_xy<float> point_2d;

	//const string DataTypes[] = {"Unknown","Byte","UInt16","Int16","UInt32","Int32","Float32","Float64","CInt16","CInt32","CFloat32","CFloat64"};

    // Multiply together all permutations
    /*GeoImage Permutations(const GeoImage& image1, const GeoImage& image2, string filename) {
        GeoImageIO<float> imageout(GeoImage(filename, image1, GDT_Float32, image1.NumBands()*image2.NumBands()));
        CImg<float> cimgout;
        int bandout(0);
        CImg<unsigned char> mask;
        float nodataout = -32768;
        imageout.SetNoData(nodataout);

        std::vector<bbox>::const_iterator iChunk;
        for (unsigned int i1=0; i1<image1.NumBands(); i1++) {
            GeoRasterIO<float> img1(image1[i1]);
            std::vector<bbox> Chunks = img1.Chunk();
            for (unsigned int i2=0; i2<image2.NumBands(); i2++) {
                GeoRasterIO<float> img2(image2[i2]);
                imageout[bandout].SetDescription(img1.Description()+'-'+img2.Description());
                for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
                    cimgout = img1.Read(*iChunk).mul(img2.Read(*iChunk));
                    mask = img1.NoDataMask(*iChunk)|=(img2.NoDataMask(*iChunk));
                    cimg_forXY(mask,x,y) if (mask(x,y)) cimgout(x,y) = nodataout;
                    imageout[bandout].Write(cimgout, *iChunk);
                }
            bandout++;
            }
        }
        return imageout;
    }*/

	//! Return string of all supported formats
	string Formats() {
		std::stringstream str;
		for(int i=0; i<GDALGetDriverCount(); i++ )
		{
			GDALDriverH hDriver = GDALGetDriver(i);
			const char *pszRWFlag;
			if( GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL ) )
				pszRWFlag = "rw+";
			else if( GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATECOPY, NULL ) )
				pszRWFlag = "rw";
			else
				pszRWFlag = "ro";
			str << GDALGetDriverShortName( hDriver ) << " (" << pszRWFlag << ") - "
					<< GDALGetDriverLongName( hDriver ) << std::endl;
		}
		return str.str();
	}

	//! Pixels of interest
	class POI {
	public:
		POI() {}
		~POI() {}
		bgeo::model::box<point_2d> Extent() const {
			//cout << "POI " << _xmin << ", " << _ymin << " - " << _xmax << ", " << _ymax << endl;
			return bgeo::model::box<point_2d>(point_2d(_xmin,_ymin),point_2d(_xmax,_ymax));
		}

		//! Get point
		point_2d operator[](int index) { return _Locs[index]; }
		const point_2d& operator[](int index) const { return _Locs[index]; }

		int NumPixels() const { return _Locs.size(); }
		void AddPixel(int x, int y) { AddPixel( point_2d(x,y) ); }
		void AddPixel(point_2d point) {
			_Locs.push_back(point);
			int x = point.x();
			int y = point.y();
			if (_Locs.size() == 0) {
				_xmax = _xmin = x;
				_ymax = _ymin = y;
			} else {
				if (x > _xmax) _xmax = x;
				if (x < _xmin) _xmin = x;
				if (y > _ymax) _ymax = y;
				if (y < _ymin) _ymin = y;
			}
		}
	private:
                vector<point_2d> _Locs;
		int _xmax, _xmin, _ymax, _ymin;	
        };

	//! Statistics for Pixels of Interest
/*	template<typename T> class GeoStats {
	public:
		GeoStats()
			: _Stats(12,1,1,1,0), _NumPixels(0) {}
		GeoStats(const GeoRaster& Image) {
			_NumPixels = Image.Size();
			//cout << Image.GetGeoImage()->Basename() << ": Calculating Stats" << endl;
			vector< bgeo::model::box<point_2d> > Chunks = Image.GetGeoImage()->Chunk();
			// Read in first chunk and get stats
			CImg<T> img = Image.Read(Chunks[0]);
			CImg<T> tmpstats;
			_Stats = img.stats();
			//cout << "chunk1 mean = " << _Stats[2] << endl;
			// weight by relative chunk size
			_Stats[2] = _Stats[2] * bgeo::area(Chunks[0]) / Image.Size();
			_Stats[3] = img.get_mul(img).sum() / Image.Size();
	                vector< bgeo::model::box<point_2d> >::const_iterator iChunk;
                	for (iChunk=Chunks.begin()+1; iChunk!=Chunks.end(); iChunk++) {
				img = Image.Read(*iChunk);
				tmpstats = img.stats();
				_Stats[2] = _Stats[2] + ( tmpstats[2] * bgeo::area(*iChunk) / Image.Size() );
				_Stats[3] = _Stats[3] + ( img.get_mul(img).sum() / Image.Size() );
				if (tmpstats[0] < _Stats[0]) {
					_Stats[0] = tmpstats[0];
					_Stats[4] = tmpstats[4];
					_Stats[5] = tmpstats[5];
				}
				if (tmpstats[1] > _Stats[1]) {
					_Stats[1] = tmpstats[1];
					_Stats[8] = tmpstats[8];
					_Stats[9] = tmpstats[9];
				}
				_Stats[3] = _Stats[3] - (_Stats[2]*_Stats[2]);
				//cout << "chunk mean = " << _Stats[2] << endl;
			}
		}
		GeoStats(const GeoRaster& Image, const POI& pixels, bool meanonly=false, GeoRaster *Mask=NULL)
			: _Pixels(pixels.NumPixels()), _NumPixels(0), _NumTossedPixels(0), _Stats() {
			bgeo::model::box<point_2d> extent = pixels.Extent();
			//cout << "read extent" << endl;
			CImg<T> img = GeoRasterIO<T>(Image).Read(extent);
			CImg<T> mask;  if (Mask!=NULL) mask = Mask->Read(extent);
			//cout << "NumPixels" << pixels.NumPixels() << endl;
                        for (int i=0;i<pixels.NumPixels(); i++) {
				// Coordinates on img
				point_2d p( pixels[i].x() - extent.min_corner().x(), pixels[i].y() - extent.min_corner().y());
                                T val = img(p.x(),p.y());
                                //cout << "Pixel " << p.x() << ", " << p.y() << ", value = " << val << endl;
				if (mask == NULL)
					_Pixels(_NumPixels++) = val;
				else {
					if (mask(p.x(),p.y()) > 0) _NumTossedPixels++;
					else _Pixels(_NumPixels++) = val;
				}		 
                             //if (!Image.NoData() || (val != Image.NoDataValue())) {_Pixels(_NumPixels++) = val;
                               // } else if (Image.NoData() && (val == Image.NoDataValue())) {_NumNoDataPixels++;	}

                        // Crop remaining array
			if (_NumPixels != 0)
				_Pixels.crop(0,0,_NumPixels-1,0);
			else _Pixels.crop(0,0,0,0);
			//for (int i=0; i<_Pixels.size();i++) cout << " " << _Pixels[i];
			//cout << endl;
			switch (_NumPixels) {
				case 0: break;
				case 1: _Stats = CImg<double>(12,1,1,1,0);
					_Stats[0] = _Stats[1] = _Stats[2] = _Pixels[0];
					break;
				default: 
					if (meanonly) {
						_Stats = CImg<double>(12,1,1,1,0);
						_Stats[2] = _Pixels.sum()/_NumPixels;
					} else _Stats = _Pixels.stats();
					break;
			}
			//cout << "Mean " << Mean() << endl;
                }
		~GeoStats() {}

		T Min() const { return _Stats[0]; }
		T Max() const { return _Stats[1]; }
		double Mean() const { return _Stats[2]; }
		double Variance() const { return _Stats[3]; }	
		double StdDev() const { return sqrt(_Stats[3]); }	

		//! Number of pixels used for statistics
		int NumPixels() const { return _NumPixels; }
		//! Number of pixels not used due to No Data value
		int NumTossedPixels() const { return _NumTossedPixels; }

		const CImg<double>& Stats() const { return _Stats; }

	private:
		CImg<double> _Stats;
		CImg<T> _Pixels;
		int _NumPixels;
		int _NumTossedPixels;
	}; // class GeoStats*/
} //namespace ags
#endif
