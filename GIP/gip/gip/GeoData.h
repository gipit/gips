#ifndef GIP_GEODATA_H
#define GIP_GEODATA_H

#include <vector>
#include <string>
//#include <map>
#include <gdal/gdal_priv.h>
#include <boost/shared_ptr.hpp>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/filesystem.hpp>

#include <gip/Utils.h>

namespace gip {
	//class GeoData : boost::enable_shared_from_this<GeoData> {
	typedef boost::geometry::model::d2::point_xy<float> point;
	typedef std::map<std::string,std::string> dictionary;

	class GeoData {
	public:

		//! \name Constructors/Destructor
		//! Default constructor
		GeoData() : _GDALDataset() {}
		//! Open existing file
		GeoData(std::string, bool=true);
		//! Create new file on disk
		GeoData(int, int, int, GDALDataType, std::string, dictionary = dictionary());
		//! Copy constructor
		GeoData(const GeoData&);
		//! Assignment copy
		GeoData& operator=(const GeoData&);
		//! Destructor should flush GDALDataset if last open pointer
		~GeoData();

		//! \name File Information
		//! Full filename of dataset
		//boost::filesystem::path Filename() const { return _Filename; }
		std::string Filename() const { return _Filename.string(); }
		//! Filename without path
		std::string Basename() const { return _Filename.stem().string(); }
		//! File format of dataset
		std::string Format() const { return _GDALDataset->GetDriver()->GetDescription(); }
		//! Get Product name (used for determining Colors, etc
		std::string Product() const {
			std::string product = GetMeta("SHORTNAME");
			// TODO - Other product names?
			return product;
		}
		//! Get product-specific options
		//Options ProductOptions();

		//! \name Spatial Information
		//! X Size of image/raster, in pixels
		unsigned int XSize() const { return _GDALDataset->GetRasterXSize(); }
		//! Y Size of image/raster, in pixels
		unsigned int YSize() const { return _GDALDataset->GetRasterYSize(); }
		//! Total number of pixels
		unsigned long Size() const { return XSize() * YSize(); }
		//! Geolocated coordinates of a pixel
		point GeoLoc(float xloc, float yloc) const;
		//! Coordinates of top left
		point TopLeft() const { return GeoLoc(0,0); }
		//! Coordinates of bottom right
		point LowerRight() const { return GeoLoc(XSize()-1,YSize()-1); }

		//! \name Metadata functions
		//! Get metadata item
		std::string GetMeta(std::string key) const {
			const char* item = _GDALDataset->GetMetadataItem(key.c_str());
			if (item == NULL) return ""; else return item;
		}
		//! Set metadata item
		GeoData& SetMeta(std::string key, std::string item) {
			_GDALDataset->SetMetadataItem(key.c_str(),item.c_str());
			return *this;
		}
		//! Get metadata group
		std::vector<std::string> GetMetaGroup(std::string group,std::string filter="") const;
		//! Copy Meta data from input file.  Currently no error checking
		GeoData& CopyMeta(const GeoData& img);
		//! Copy collection of meta data
		//GeoData& CopyMeta(const GeoData&, std::vector<std::string>);
		//! Copy coordinate system
		GeoData& CopyCoordinateSystem(const GeoData&);

		//! \name Processing functions
		//! Break up image into chunks
		std::vector< boost::geometry::model::box<point> > Chunk(int overlap=0, unsigned int bytes=2) const;

		//! Add overviews
		GeoData& AddOverviews() {
            int anOverviewList[3] = { 2, 4, 8 };
            _GDALDataset->BuildOverviews( "NEAREST", 3, anOverviewList, 0, NULL, GDALDummyProgress, NULL );
            return *this;
		}

		//! Flush cache
		void Flush() { _GDALDataset->FlushCache(); }

	protected:
		//! Filename to dataset
		boost::filesystem::path _Filename;
		//! Underlying GDALDataset of this file
		boost::shared_ptr<GDALDataset> _GDALDataset;
		//! Options
		//Options _Options;
	}; //class GeoData

} // namespace ags

#endif
