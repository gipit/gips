#ifndef GIP_GEOIMAGE_H
#define GIP_GEOIMAGE_H

#include <gip/Colors.h>
#include <gip/GeoData.h>
#include <gip/GeoRaster.h>
#include <gip/GeoRasterIO.h>

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
		//! Constructor for creating new file with same properties (xsize, ysize, bsize,datatype) as existing file
		explicit GeoImage(std::string filename, const GeoImage& image) :
			GeoData(image.XSize(), image.YSize(), image.NumBands(), image.DataType(), filename) {
			//if (datatype == GDT_Unknown) datatype = image->GetDataType();
			CopyMeta(image);
			CopyCoordinateSystem(image);
			LoadBands();
			//_Colors = image.GetColors();
		}

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
		GDALDataType DataType() const;
		//! Return information on image as string
		std::string Info(bool=true, bool=false) const;
		//! Get vector of band names
		std::vector<std::string> BandNames() const;

		//! Retrieve colors class
		Colors GetColors() const { return _Colors; }
		//! Set a color
		void SetColor(std::string col,int bandnum) {
		    _Colors.SetColor(col,bandnum);
		    // Set color on individual band
            _RasterBands[bandnum-1].SetColor(col);
            _RasterBands[bandnum-1].SetDescription(col);
        }
		//! Set colors (blue,green,red,nir,swir1,swir2,lwir)
        void SetColors(int blue, int green, int red, int nir=0) {
            SetColor("Blue",blue);
            SetColor("Green",green);
            SetColor("Red",red);
            if (nir > 0) SetColor("NIR",nir);
        }


		// \name Band Operations
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

		// Raster functions to run on all bands
		const GeoImage& ComputeStats() const;

        //! Set gain for all bands
        void SetGain(float gain) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetGain(gain); }
        //! Set gain for all bands
        void SetOffset(float offset) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetOffset(offset); }

		//! Set NoData for all bands
		void SetNoData(double val) { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].SetNoData(val); }
		//! Unset NoData for all bands
		void ClearNoData() { for (unsigned int i=0;i<_RasterBands.size();i++) _RasterBands[i].ClearNoData(); }

	protected:
		//! Vector of raster bands
		std::vector< GeoRaster > _RasterBands;

		//! Map indicating which bands are what colors
		Colors _Colors;

		//! Loads Raster Bands of this GDALDataset into _RasterBands vector
		void LoadBands();

		//! Update colors in file, or at least RGB
		//void SaveColors();
	private:
		//! Default constructor, private so cannot be called
		explicit GeoImage() : GeoData() {}
	};

} // namespace gip
#endif
