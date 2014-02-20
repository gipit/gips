#ifndef GIP_GEOMASK_H
#define GIP_GEOMASK_H

#include <gip_GeoRaster.h>

namespace gip {

	//template<typename T> class GeoMask;
	class GeoImage;
	//template<typename T> class GeoRaster;

	//! The GeoMask class allows chaining together of multiple masks
	template<typename T> class GeoMask {
	public:
		//! Default constructor
		GeoMask() : _Masks(0), _th(0), _invert(0) {}
		//! Constructor takes in ptr to a GeoRaster, a threshold, and an invert flag
		GeoMask(GeoRaster* mask, double th, bool invert=false) {
			if (!mask == NULL)
				AddMask(mask,th,invert);
		}
		//! Destructor
		~GeoMask() { Release(); }

		void Release() {
			//for (int i=0;i<_Masks.size();i++) {
	//			_Masks[i]->Release();
			//}
		}

		GeoMask& AddMask(GeoRaster* mask, double th, bool invert=false) {
			if (mask!=NULL) {
				_Masks.push_back(mask);
				_th.push_back(th);
				_invert.push_back(invert);
				mask->GetGeoImage()->Reference();
			}
			return *this;
		}

		bool Valid() const {
			if (_Masks.size() == 0) 
				return false;
			else if (_Masks[0] == NULL)
				return false;
			return true;
		}

		CImg<unsigned char> Read(const bgeo::model::box<point_2d>& chunk) const {
			CImg<unsigned char> Mask(chunk.max_corner().x()-chunk.min_corner().x()+1, chunk.max_corner().y()-chunk.min_corner().y()+1,1,1,1);
			for (int i=0;i<_Masks.size();i++) {
				// Read as normal, don't apply gain/offset
				CImg<T> tmpMask = _Masks[i]->Read(chunk,false);
				tmpMask.threshold(_th[i]+1);	
				if (_invert[i]) tmpMask = (tmpMask-1).abs();
				Mask = Mask.mul(tmpMask);
			}	
			return Mask;	
		}

		GeoImage<T> Write(string filename) const {
			if (_Masks.empty()) return GeoImage<T>();
			GeoImage<T>* In = _Masks[0]->GetGeoImage();
			GeoImage<T> Out(In->XSize(),In->YSize(),1,filename,"GTiff",GDT_Byte);
			Out.CopyMetadata(*In);
			vector< bgeo::model::box<point_2d> > Chunks = Out.Chunk();
			vector< bgeo::model::box<point_2d> >::const_iterator iChunk;
			for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
				CImg<T> img = Read(*iChunk);
				Out[0].Write( img, *iChunk );
			}
			Out[0].GetGDALRasterBand()->SetDescription("Mask");
			return Out;
		}

	private:
		//! Vector of pointers
		vector<GeoRaster*> _Masks;
		//! Vector of thresholds
		vector<double> _th;
		//! Vector of invert flags
		vector<bool> _invert;

	}; // class GeoMask

} // namespace ags

#endif
