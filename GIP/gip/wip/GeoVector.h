#ifndef GIP_GEOVECTOR_H
#define GIP_GEOVECTOR_H

#include "ogrsf_frmts.h"
#include "ogr_api.h"

#include <iostream>

using namespace std;
namespace bfs = boost::filesystem;

// windows.h defines min and max macros which clobbers proper functions
//#undef min
//#undef max

                /*GeoImage& AddVectorLayer(string filename, int layer=0) {
                        bfs::path vfilename(filename);
                        OGRDataSource *poDS;
                        poDS = OGRSFDriverRegistrar::Open(filename.c_str(),FALSE);
                        if (poDS == NULL) cout << filename << ": Error opening." << endl;
                        OGRLayer *poLayer;
                        poLayer = poDS->GetLayer(layer);
                        poLayer->ResetReading();

                        // Create coordinate transform between vector and GDALDataset coordinate system
                        OGRSpatialReference *SRSout = new OGRSpatialReference(_GDALDataset->GetProjectionRef());
                        OGRCoordinateTransformation *transform = OGRCreateCoordinateTransformation(poLayer->GetSpatialRef(),SRSout);
                        if (transform == NULL) cout << vfilename.leaf() << ": could not transform coordinate system to match raster." << endl;
                        OGRSpatialReference::DestroySpatialReference(SRSout);
                        
                        // Create new file for rasterization
                        GeoImage vImage(*this,1,"GTiff",vfilename.stem());
                        
                        // Loop through features
                        OGRFeature *poFeature;
                        double fid=1;
                        while( (poFeature = poLayer->GetNextFeature()) != NULL ) {
                                OGRGeometry *poGeometry;
                                poGeometry = poFeature->GetGeometryRef();
                                // apply transformation
                                poGeometry->transform(transform);
                                // rasterize into new file 
                                int bands[] = {1};
                                GDALRasterizeGeometries(vImage.GetGDALDataset(),1,&bands[0],1,(void**)&poGeometry,NULL,NULL,&fid,NULL,NULL,NULL);
                                fid++;
                                OGRFeature::DestroyFeature(poFeature);  
                        }
                        OGRCoordinateTransformation::DestroyCT(transform);

                        return *this;
                }*/



namespace gip {

	class GeoVector {
	public:
		//! \name Constructors / Assignment / Destructor 
		// @{
		//! Default Constructor.
		/*!
			The default constructor should be used with caution.  A GVector created with the default constructor
			is not a valid vector file and should not be used until it is assigned to an object. 
			Create a new GVector.  Will create a new data file on disk, with supplied format and size
			\param Options - A GIP::Options structure containing user supplied options
		*/
		//GVector(Options opts=Options()) 
		//	: m_Filename( bfs::path(opts.workdir/UniqueName()) ), m_Dataset(NULL) {}
		//! Existing GVector Constructor.
		/*!
			This constructor provides a filename for a Geospatial vector file to open.  An exception is thrown if
			the file doesn't exist.  The dataset is opened as Read Only
		*/
		GeoVector(string filename ) 
			: m_Filename(filename) {
			Open();
		}
		//! Copy constructor
		/*!
			Copy vector file
			\param image - The GVector to be copied
			\return Newly copied GVector
		*/
		GeoVector(const GeoVector& vector)
			: m_Filename(vector.GetDataset()->GetName()), m_Dataset(vector.GetDataset()) {
			int cnt = m_Dataset->Reference();
			//cout << Basename() << ": Copy - " << cnt << " references" << endl;
		}
                static GeoVector CreateNew(string filename = "") {
                //GVector (Options opts = Options(), string filename = "") {
                        OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");

                        if (filename == "") filename = "tmp.shp"; //(opts.workdir/UniqueName()).string() + ".shp";
                        stringstream err;
                        if( poDriver == NULL ) {
                                cout << "Unable to find driver " << " for creating new vector file" << endl;
                                throw std::exception(); //err.str().c_str());
                        }
                        if( !poDriver->TestCapability( ODrCCreateDataSource ) ) {
                                cout << "Driver does not support data source creation." << endl;
                                throw std::exception(); //err.str().c_str());
                        }
                        // Current bug in ogr - driver not attached to returned instance
                        OGRDataSource *poODS = poDriver->CreateDataSource( filename.c_str() );
                        poODS->SetDriver(poDriver);
                        if( poODS == NULL ) {
                                cout << "dRiver failed to create " << filename.c_str() << endl;
                                throw std::exception(); //err.str().c_str());
                        }
                        poODS->Reference();
                        cout << (bfs::path(filename)) << ": Creating new " << poDriver->GetName() << endl;
                        //m_Filename = filename;
                        //m_Dataset = poODS;
			OGRDataSource::DestroyDataSource(poODS);
                        return GeoVector(filename); //poODS, opts);
                }

		//! Assignment operator
		GeoVector& operator=(const GeoVector &vectors) {
			if (this == &vectors) return *this;
			// Remove old dataset
			Release();
			// Assign new dataset
			m_Filename = vectors.Filename();
			m_Dataset = vectors.GetDataset();
			int cnt = m_Dataset->Reference();
			return *this;
		}
		//! Destructor.
		~GeoVector() { Release();	}
		// @}

		//! \name File Information
		//! The filename of the dataset this image came from
		string Filename() const { return m_Filename.string(); }
		//! The base filename of the dataset
		string Basename() const { return m_Filename.leaf(); }
		//! File format name of input file
		string Format() const { return string(m_Dataset->GetDriver()->GetName()); }
		//! Return pointer to OGR Dataset
		OGRDataSource* GetDataset() const { return m_Dataset; }
		// @}

		//! \name Processing
		//! Copy Layers between datasets, transform if supplied
		void CopyLayers(OGRDataSource *src, OGRDataSource *dst, OGRSpatialReference *poOutputSRS = NULL, OGRCoordinateTransformation *poCT = NULL) {
			stringstream err;
			// Must step through here and transform each feature in each layer
			for(int iLayer = 0; iLayer < src->GetLayerCount(); iLayer++ ) {
				// Source Layer
				OGRLayer *poSourceLayer = src->GetLayer(iLayer);
				if( poSourceLayer == NULL ) { cout << "Could not process layer in vector file" << endl; }
				OGRFeatureDefn *poFDefn = poSourceLayer->GetLayerDefn();
				// Output Layer
				if (poOutputSRS == NULL) poOutputSRS = poSourceLayer->GetSpatialRef();
		        OGRLayer *poOutputLayer = dst->CreateLayer( poFDefn->GetName(), poOutputSRS, poFDefn->GetGeomType() );				
				// Copy over field definitions from source for all features in layer
				for(int iField = 0; iField < poFDefn->GetFieldCount(); iField++ ) poOutputLayer->CreateField( poFDefn->GetFieldDefn(iField) );

				// Transfer features
				OGRFeature  *poFeature;
				int         nFeaturesInTransaction = 0;
			    
				poSourceLayer->ResetReading();
				poOutputLayer->StartTransaction();
				while( (poFeature = poSourceLayer->GetNextFeature()) != NULL ) {
					OGRFeature *poDstFeature = OGRFeature::CreateFeature( poOutputLayer->GetLayerDefn() );

					if( poDstFeature->SetFrom( poFeature, TRUE ) != OGRERR_NONE ) {
						poOutputLayer->CommitTransaction();
						err << "Unable to translate feature " << poFeature->GetFID() << " from layer " << poFDefn->GetName() << endl;
						OGRFeature::DestroyFeature( poFeature );
						OGRFeature::DestroyFeature( poDstFeature );
						throw std::exception(); //err.str().c_str());
					}
					// Preserve feature ID
					poDstFeature->SetFID( poFeature->GetFID() );
			        
					if( poCT && poDstFeature->GetGeometryRef() != NULL ) {
						OGRErr eErr = poDstFeature->GetGeometryRef()->transform( poCT );
						if( eErr != OGRERR_NONE ) {
							err << "Failed to transform feature " << poFeature->GetFID() << endl;
							OGRFeature::DestroyFeature( poFeature );
							OGRFeature::DestroyFeature( poDstFeature );
							throw std::exception(); //err.str().c_str());
						}
					}
					poOutputLayer->CreateFeature( poDstFeature );
					OGRFeature::DestroyFeature( poFeature );
					OGRFeature::DestroyFeature( poDstFeature );
					/*CPLErrorReset();
					if( poDstLayer->CreateFeature( poDstFeature ) != OGRERR_NONE ) {
						if( nGroupTransactions )
							poDstLayer->RollbackTransaction();
						err << "Unknown problem with vector features " << poFeature->GetFID() << endl;
						OGRFeature::DestroyFeature( poDstFeature );
						throw std::exception(err.str());
						return FALSE;
					}*/	
				}
				poOutputLayer->CommitTransaction();
			}
		}

		//! Transform vector file to new projection, overwrite existing
		GeoVector& Transform(string wkt) {
			*this = get_Transform(wkt);
			return *this;
		}

		//! Transform the vector file to new projection, don't overwrite existing, returns new
		GeoVector get_Transform(string wkt) {
			stringstream err;
			cout << Filename() << ": Transforming:" << endl << wkt << endl;

			// This creates a new vector file in the working directory
			//cout << "Creating new" << endl;
			GeoVector newfile = GeoVector::CreateNew();
			//cout << "end create new" << endl;

			// Get source and output references and calculate transformation - using single transform for all layers
			bool delSourceSRS = false;
			OGRSpatialReference *poSourceSRS = m_Dataset->GetLayer(0)->GetSpatialRef();
			if (poSourceSRS == NULL) {
				delSourceSRS = true;
				poSourceSRS = new OGRSpatialReference(); //(OGRSpatialReference *)OSRNewSpatialReference( NULL );
				poSourceSRS->SetWellKnownGeogCS( "WGS84" );
			}
			OGRSpatialReference *poOutputSRS = new OGRSpatialReference(wkt.c_str()); //(OGRSpatialReference *)OSRNewSpatialReference( wkt.c_str() );
			OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation( poSourceSRS, poOutputSRS );
			/*if (m_Options.verbose) {
				char *pszWKT = NULL;
				poSourceSRS->exportToPrettyWkt( &pszWKT, FALSE ); cout << "Source SRS: " << pszWKT << endl;
				poOutputSRS->exportToPrettyWkt( &pszWKT, FALSE ); cout << "Output SRS: " << pszWKT << endl;
			}*/
			if (poCT == NULL) {
				cout << Basename() << " - could not transform coordinate system of vector file." << endl;
				throw std::exception(); //err.str().c_str());
			}
			CopyLayers(m_Dataset, newfile.GetDataset(), poOutputSRS, poCT);
			//CopyLayers(m_Dataset, newfile.GetDataset());
			//if (delSourceSRS) //delete poSourceSRS;
			//	OSRDestroySpatialReference( poSourceSRS );
			//delete poOutputSRS; 
			//OSRDestroySpatialReference( poOutputSRS );
			delete poCT;
			return newfile;
		}
		// @}

	protected:
		//! \name Data
		// @{
		//! The filename of the file
		bfs::path m_Filename;
		//! The underlying OGR Dataset
		OGRDataSource *m_Dataset;
		// @}

		//! Private constructor using OGRDataSource input
		/*!
			This constructor creates a new GVector given an already existing and open OGRDatasource
		*/
		GeoVector(OGRDataSource* ds )
			: m_Filename(ds->GetName()) {
			m_Dataset = ds;
		}
		void Open() {
			// Blank filename is OK, just means this is an empty image
			if ( m_Filename == "" ) return;
			m_Dataset = OGRSFDriverRegistrar::Open( m_Filename.string().c_str(), FALSE );
			if (m_Dataset == NULL) {
				stringstream err;
				cout << "Error opening " << m_Filename << endl;
				throw std::exception(); //err.str().c_str());
			}
			cout << Basename() << ": Opening as " << Format() << " (" << m_Dataset->GetLayerCount() << " layers)" << endl;
		}
		//! Release underlying dataset if there are no more references to it
		void Release() {
			int cnt = m_Dataset->Dereference();
			if (cnt == 0) {
				OGRSFDriver *driver = m_Dataset->GetDriver();
				delete m_Dataset;
				cout << Basename() << ": Freeeing dataset" << endl;
				// Delete underlying file if option set
				//if (m_Options.verbose) cout << Basename() << ": Deleting datafile" << endl;
				//driver->DeleteDataSource(Filename().c_str());
			}
		}
	}; // class GVector
}

#endif
