/*!
 * \ingroup ags_utilities
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <agsOptions.h>
#include <agsGeoImage.h>
#include <agsGeoAlgorithm.h>
#include <agsGeoMask.h>

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;
using namespace ags;

int main (int ac, char* av[]) {
	typedef float T;
	// Setup program options and parse command line
	boost::program_options::options_description localopts("Raster statistics options");
	localopts.add_options()
		("vector",bpo::value<string>()->default_value(""), "Rasterized vector file of regions of interest")
		("cloudband,c",bpo::value<int>(), "Band number of cloud classification data")
		("id",bpo::value<int>(), "Vector ID to calculate stats for")
		("nodata",bpo::value<T>(), "Value of nodata regions");
	ags::agsOptions Opts(localopts);
	try {
		Opts.ParseCmdLine(ac,av);
	} catch (std::exception& err) {
		cout << "Error: " << err.what() << endl;
		return 1;
	}
	if (Opts.Exit()) return 0;
	// End parse command line
	GDALAllRegister();
	CPLSetConfigOption("GDAL_PAM_ENABLED",Opts.PAM() ? "TRUE" : "FALSE");

	// Options
	string vfilename = Opts.varmap["vector"].as<string>();
	bool nodata = Opts.varmap.count("nodata") ? true : false;
	T nodatavalue;
	if (nodata) nodatavalue = Opts.varmap["nodata"].as<T>();
	int cloudband = Opts.varmap.count("cloudband") ? Opts.varmap["cloudband"].as<int>() : 0;

	// Main algorithm
	// vector file
	GeoImage<T> vfile(vfilename);
	GeoStats<T> vStats(vfile[0]);
	T maxID = vStats.Max();
	GeoRaster<T> *cloud;

	int p_lo = 1;
	int p_hi = maxID;
	if (Opts.varmap.count("id")) {
		int id = Opts.varmap["id"].as<int>();
		p_lo = id;
		p_hi = id;
	}

	int wsm(8);
	int w(12);
	cout << setw(wsm) << "WELDDOY";
	for (int b=0;b<Opts.Bands().size();b++) {
		cout << setw(wsm) << ("b"+to_string(Opts.Bands(b))+"_NumP")
			<< setw(w) << ("b" + to_string(Opts.Bands(b)) + "_Mean")
			<< setw(w) << ("b" + to_string(Opts.Bands(b)) + "_StdDev");
	}
	cout << endl;

	// Loop through input files
	for (int i=0;i<Opts.InputFiles().size(); i++) {
		GeoImage<T> Image(Opts.InputFiles(i),Opts.Gain(),Opts.Offset());
		vector< GeoMask<T> > Masks;
		vector< map<T,POI> > POIs;
		// data specific 
		Image[5].SetGainOffset(0.01);
		Image[6].SetGainOffset(0.01);
		string year= Image.Basename().substr(13,4);
		string doy = Image.Basename().substr(33,3);

		// Pointer to cloud band
		cloud = (cloudband > 0) ? &(Image[cloudband]) : NULL;
		GeoMask<T> mask(cloud,1,true);

		// Map Pixels of Interest
		for (int b=0;b<Opts.Bands().size(); b++) {
			// Set NoData value, if there is one, and map POIs for remaining good pixels
			Masks.push_back( GeoMask<T>(cloud,1,true) );
			if (nodata) {
				Image[Opts.Bands(b)-1].SetNoData(nodatavalue);
				Masks[b].AddMask(&Image[Opts.Bands(b)-1],nodatavalue,false);
			}
			POIs.push_back( vfile[0].MapPOIs(Masks[b]) );
		}
		// Loop through POI regions
		for (int p=p_lo;p<=p_hi;p++) {
			vector< GeoStats<T> > Stats;
			int numpixels = 0;
			for (int b=0;b<Opts.Bands().size(); b++) {
				map<T,POI>::iterator iPOI = POIs[Opts.Bands(b)-1].find(p);
				if (iPOI != POIs[Opts.Bands(b)-1].end())
					Stats.push_back( GeoStats<T>(Image[Opts.Bands(b)-1], iPOI->second ));
				else Stats.push_back( GeoStats<T>() );
				numpixels += Stats[b].NumPixels();
			}
			if (numpixels > 0) {
				// ID, Year, doy
				cout << setw(wsm) << doy;
				for (int b=0;b<Opts.Bands().size();b++)
					cout << setw(wsm) << Stats[b].NumPixels() << setw(w) << Stats[b].Mean() << setw(w) << Stats[b].StdDev();
				cout << endl;
			}
		}
	}
	// Two input files required: raster, then vector
	//Image.AddVectorLayer(Opts.InputFiles(1));
	
	return 0;
}
