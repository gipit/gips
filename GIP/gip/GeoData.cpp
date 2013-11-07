/*
 * gip_GeoData.cpp
 *
 *  Created on: Aug 26, 2011
 *  Author: mhanson
 */

#include <gip/GeoData.h>
#include <boost/make_shared.hpp>

#include <iostream>

namespace gip {
    using std::string;
    using std::vector;
    using std::map;
    using boost::filesystem::path;
    typedef boost::geometry::model::d2::point_xy<float> point;
    typedef boost::geometry::model::box<point> bbox;

	boost::filesystem::path Options::_ConfigDir("/usr/share/gip/");
	string Options::_DefaultFormat("GTiff");
	float Options::_ChunkSize(128.0);
	int Options::_Verbose(1);
	std::string Options::_WorkDir("/tmp/");

	/*!
	 * Open an existing file and returns a shared pointer
	 */
	GeoData::GeoData(string filename, bool Update) : _Filename(filename) {
		GDALAccess access = Update ? GA_Update : GA_ReadOnly;
		GDALDataset* ds = (GDALDataset*)GDALOpenShared(_Filename.string().c_str(), access);
		// Check if Update access not supported
		if (ds == NULL && CPLGetLastErrorNo() == 6)
			ds = (GDALDataset*)GDALOpenShared(_Filename.string().c_str(), GA_ReadOnly);
		if (ds == NULL) {
			throw std::runtime_error(to_string(CPLGetLastErrorNo()) + ": " + string(CPLGetLastErrorMsg()));
		}
		_GDALDataset.reset(ds);
		//ProductOptions();
		//std::cout << "GeoData Open (use_count = " << _GDALDataset.use_count() << ")" << std::endl;
	}

	/*!
	 * Creates new file on disk and returns shared pointer
	 */
	GeoData::GeoData(int xsz, int ysz, int bsz, GDALDataType datatype, string filename, dictionary options)
		:_Filename(filename) {
		string format = Options::DefaultFormat();
		//if (format == "GTiff" && datatype == GDT_Byte) options["COMPRESS"] = "JPEG";
		GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(format.c_str());
		// TODO check for null driver and create method
		// Check extension
		string ext = driver->GetMetadataItem(GDAL_DMD_EXTENSION);
		if (ext != "" && _Filename.extension().string() != ('.'+ext)) _Filename = path(_Filename.string() + '.' + ext);
		char **papszOptions = NULL;
		if (options.size()) {
            for (dictionary::const_iterator imap=options.begin(); imap!=options.end(); imap++)
                papszOptions = CSLSetNameValue(papszOptions,imap->first.c_str(),imap->second.c_str());
		}
		_GDALDataset.reset( driver->Create(_Filename.string().c_str(), xsz,ysz,bsz,datatype, papszOptions) );
		if (_GDALDataset.get() == NULL)
			std::cout << "Error creating " << _Filename.string() << CPLGetLastErrorMsg() << std::endl;
	}

	/*!
	 * Copy constructor
	 */
	GeoData::GeoData(const GeoData& geodata)
		: _Filename(geodata._Filename), _GDALDataset(geodata._GDALDataset) {
	}

	/*!
	 * Assignment copy
	 */
	GeoData& GeoData::operator=(const GeoData& geodata) {
	//GeoData& GeoData::operator=(const GeoData& geodata) {
		// Check for self assignment
		if (this == &geodata) return *this;
		_Filename = geodata._Filename;
		_GDALDataset = geodata._GDALDataset;
		return *this;
	}

	/*!
	 * Destructor
	 */
	GeoData::~GeoData() {
		//std::cout << "GeoData::~GeoData " << Basename() << "(use_count " << _GDALDataset.use_count() << ")" << std::endl;
		if (_GDALDataset.unique()) { _GDALDataset->FlushCache(); }
	}

	/*!
	 * Retrieves and loads Options based on current product, if available
	 */
	/*Options GeoData::ProductOptions() {
		if (!Product().empty()) {
			_Options = Options(Options::ConfigDir() + Product() + ".gipcfg");
		} else _Options = Options();
		return _Options;
	}*/

	/*!
	 * Using GDALDatasets GeoTransform get Geo-located coordinates
	 */
	point GeoData::GeoLoc(float xloc, float yloc) const {
		double Affine[6];
		_GDALDataset->GetGeoTransform(Affine);
		point Coord(Affine[0] + xloc*Affine[1] + yloc*Affine[2], Affine[3] + xloc*Affine[4] + yloc*Affine[5]);
		return Coord;
	}

	/*!
	 * Get group of metadata items
	 */
	vector<string> GeoData::GetMetaGroup(string group,string filter) const {
		char** meta= _GDALDataset->GetMetadata(group.c_str());
		int num = CSLCount(meta);
		vector<string> items;
		for (int i=0;i<num; i++) {
				if (filter != "") {
						string md = string(meta[i]);
						string::size_type pos = md.find(filter);
						if (pos != string::npos) items.push_back(md.substr(pos+filter.length()));
				} else items.push_back( meta[i] );
		}
		return items;
	}

	/*!
	 * Copy all meta data, or what is in Options
	 */
	GeoData& GeoData::CopyMeta(const GeoData& img) {
		//_GDALDataset->SetDescription(img.GetGDALDataset()->GetDescription());
		//vector<string> meta = _Options.Meta();
		//if (meta.empty())
		// This now just copies all metadata
			_GDALDataset->SetMetadata(img._GDALDataset->GetMetadata());
		/*else {
			for (vector<string>::const_iterator iMeta=meta.begin();iMeta!=meta.end();iMeta++) {
				_GDALDataset->SetMetadataItem(iMeta->c_str(),img.GetMeta(*iMeta).c_str());
			}
		}*/
		return *this;
	}

	/*!
	 * Copies coordinate system from another GeoData image
	 */
	GeoData& GeoData::CopyCoordinateSystem(const GeoData& img) {
		GDALDataset* ds = const_cast<GeoData&>(img)._GDALDataset.get();
		_GDALDataset->SetProjection(ds->GetProjectionRef());
		double Affine[6];
		ds->GetGeoTransform(Affine);
		_GDALDataset->SetGeoTransform(Affine);
		return *this;
	}

	/*!
	 * Chunk() breaks up the image into smaller size pieces, with
	 * each piece being no bigger than GeoData::_ChunkSize
	 */
	vector<bbox> GeoData::Chunk(int overlap, unsigned int bytes) const {
        // TODO - overlap!!  and bytes is fixed at 2!!!
		unsigned int rows = floor( (Options::ChunkSize()*1024*1024) / bytes / XSize() );
		rows = rows > YSize() ? YSize() : rows;
		int numchunks = ceil( YSize()/(float)rows );
		vector<bbox> Chunks;
		for (int i=0; i<numchunks; i++) {
			point p1(0,rows*i);
			point p2(XSize()-1, std::min((rows*(i+1)-1),YSize()-1));
			bbox chunk(p1,p2);
			//std::cout << "Chunk " << i << ": " << boost::g+eometry::dsv(chunk) << std::endl;
			Chunks.push_back(chunk);
		}
		return Chunks;
	}

}
