#ifndef GIP_GEOIMAGE_H
#define GIP_GEOIMAGE_H

#include <gip/Colors.h>
#include <gip/GeoData.h>
#include <gip/GeoRaster.h>

namespace gip {
	// Forward declaration
	class GeoRaster;

	//! GeoImage class
	/*!
		The GeoImage is a collection of GeoRaster objects
	*/
	class GeoImage : public GeoData {
		// friend so GeoRaster can access GeoImage freely
		//friend class GeoRaster;
	public:
		//! \name Constructors/Destructor
		//! Default constructor
		explicit GeoImage() : GeoData() {}
		//! Open file constructor
		explicit GeoImage(std::string filename, bool Update=true)
			: GeoData(filename, Update) {
			LoadBands();
		}
		//! Constructor for creating new file
		explicit GeoImage(std::string filename, int xsz, int ysz, int bsz, GDALDataType datatype=GDT_Byte) :
			GeoData(xsz, ysz, bsz, datatype, filename) {
			LoadBands();
		}
		//! Constructor for creating new file with same properties (xsize, ysize, metadata) as existing file
		explicit GeoImage(std::string filename, const GeoImage& image, GDALDataType datatype, int bsz) :
			GeoData(image.XSize(), image.YSize(), bsz, datatype, filename) {
			//if (datatype == GDT_Unknown) datatype = image->GetDataType();
			CopyMeta(image);
			CopyCoordinateSystem(image);
			LoadBands();
			//_Colors = image.GetColors();
		}
		//! Constructor for creating new file with same properties (xsize, ysize, bsize) as existing file
		explicit GeoImage(std::string filename, const GeoImage& image, GDALDataType datatype) :
			GeoData(image.XSize(), image.YSize(), image.NumBands(), datatype, filename) {
			//if (datatype == GDT_Unknown) datatype = image->GetDataType();
			CopyMeta(image);
			CopyCoordinateSystem(image);
			LoadBands();
			_Colors = image.GetColors();
		}
		//! Constructor for creating new file with given properties (xsize, ysize, bsize,datatype) as existing file
		explicit GeoImage(std::string filename, const GeoImage& image) :
			GeoData(image.XSize(), image.YSize(), image.NumBands(), image.DataType(), filename) {
			//if (datatype == GDT_Unknown) datatype = image->GetDataType();
			CopyMeta(image);
			CopyCoordinateSystem(image);
			LoadBands();
			//_Colors = image.GetColors();
		}
		//! Constructor to create new file based on input vector extents
		/*explicit GeoImage(std::string filename, std::string vector, float xres, float yres, GDALDataType datatype=GDT_Byte, int bsz=1) {
		    OGRDataSource *poDS = OGRSFDriverRegistrar::Open(vector.c_str());
            OGRLayer *poLayer = poDS->GetLayer(0);
		    OGREnvelope extent;
		    poLayer->GetExtent(&extent, true);
            int xsize = (int)(0.5 + (extent.MaxX - extent.MinX) / xres);
            int ysize = (int)(0.5 + (extent.MaxY - extent.MinY) / yres);
            GeoData::CreateNew(xsize, ysize, bsz, datatype, filename);
            double affine[6];
            affine[0] = extent.MinX;
            affine[1] = xres;
            affine[2] = 0;
            affine[3] = extent.MaxY;
            affine[4] = 0;
            affine[5] = -yres;
            _GDALDataset->SetGeoTransform(affine);
            char* wkt = NULL;
            poLayer->GetSpatialRef()->exportToWkt(&wkt);
            _GDALDataset->SetProjection(wkt);
            OGRDataSource::DestroyDataSource( poDS );
		}*/
		//! Copy constructor - copies GeoData and all bands
		GeoImage(const GeoImage& image);
		//! Assignment Operator
		GeoImage& operator=(const GeoImage& image) ;
		//! Destructor
		~GeoImage() { _RasterBands.clear(); }

		//! \name File Information
		//! Number of bands
		unsigned int NumBands() const { return _RasterBands.size(); }
		//! Get datatype of image (check all raster bands, return 'largest')
		GDALDataType DataType() const { return _RasterBands[0].DataType(); }
		//! Return information on image as string
		std::string Info(bool=true, bool=false) const;

		//! \name Bands and colors
		//! Get vector of band names
		std::vector<std::string> BandNames() const;
		//! Get raster band (0-based index)
		GeoRaster& operator[](int band) { return _RasterBands[band]; }
		//! Get raster band, const version
		const GeoRaster& operator[](int band) const { return _RasterBands[band]; }
		//! Get raster band by color
		GeoRaster& operator[](std::string col) {
			// Call const version
			return const_cast<GeoRaster&>(static_cast<const GeoImage&>(*this)[col]);
		}
		//! Get raster band by color, const version
		const GeoRaster& operator[](std::string col) const;
		//! Adds a band, at position bandnum (0-based)
		GeoImage& AddBand(const GeoRaster& band); //, unsigned int bandnum=0);
		//! Remove band
		GeoImage& RemoveBand(unsigned int bandnum);
		//! Prune bands to only provided colors
		GeoImage& PruneBands(std::vector<std::string>);
		//! Prune bands to RGB
		GeoImage& PruneToRGB() {
		    // TODO - get initializer_list working (c++11 standard)
		    std::vector<std::string> cols; cols.push_back("Red"); cols.push_back("Green"); cols.push_back("Blue");
		    return PruneBands(cols);
        }
		//! Retrieve colors class
		Colors GetColors() const { return _Colors; }
		//! Set a color
		void SetColor(std::string col,int bandnum) {
		    _Colors.SetColor(col,bandnum);
		    // Set color on individual band
            _RasterBands[bandnum-1].SetColor(col);
            _RasterBands[bandnum-1].SetDescription(col);
        }
		// TODO - dictionary?
        /*void SetColors(int blue, int green, int red, int nir=0) {
            SetColor("Blue",blue);
            SetColor("Green",green);
            SetColor("Red",red);
            if (nir > 0) SetColor("NIR",nir);
        }*/
        //! Copy color table from another image
		void CopyColorTable(const GeoImage& raster) {
		    if (NumBands() == 1) {
                GDALColorTable* table( raster[0].GetGDALRasterBand()->GetColorTable() );
                if (table != NULL) _RasterBands[0].GetGDALRasterBand()->SetColorTable(table);
		    }
		}

        //! \name Multiple band convenience functions
        //! Set gain for all bands
        void SetGain(float gain) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetGain(gain); }
        //! Set gain for all bands
        void SetOffset(float offset) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetOffset(offset); }
        //! Set  for all bands
        void SetUnits(std::string units) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetUnits(units); }
        //! Set UnitsOut for all bands
        void SetUnitsOut(std::string units) const { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetUnitsOut(units); }
        //! Clear atmosphere from all bands
        void ClearAtmosphere() { for(unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].ClearAtmosphere(); }
		//! Set NoData for all bands
		void SetNoData(double val) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetNoData(val); }
		//! Unset NoData for all bands
		void ClearNoData() { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].ClearNoData(); }

		//! \name Processing functions
        //! Process band into new file (copy and apply processing functions)
		GeoImage Process(std::string, GDALDataType = GDT_Unknown);

		//! Adds a mask band (1 for valid) to every band in image
		GeoImage& AddMask(const GeoRaster& band) {
		    for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].AddMask(band);
		    return *this;
        }
		//! Clear all masks
		void ClearMasks() { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].ClearMasks(); }

		//! Replace all 'Inf' or 'NaN' results with the bands NoData value
		GeoImage& FixBadPixels();

        // hmm, what's this do?
		//const GeoImage& ComputeStats() const;

		//! Add overviews
		GeoData& AddOverviews() {
            int panOverviewList[3] = { 2, 4, 8 };
            _GDALDataset->BuildOverviews( "NEAREST", 3, panOverviewList, 0, NULL, GDALDummyProgress, NULL );
            return *this;
		}

		//! Break up image into chunks
		/*void Chunk(unsigned int pad=0) const {
			GeoData::Chunk(pad);
			for (unsigned int b=0;b<NumBands();b++) _RasterBands.
		}*/

	protected:
		//! Vector of raster bands
		std::vector< GeoRaster > _RasterBands;

		//! Map indicating which bands are what colors
		Colors _Colors;

		//! Loads Raster Bands of this GDALDataset into _RasterBands vector
		void LoadBands();
	};

} // namespace gip
#endif
