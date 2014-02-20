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
		("stddev",bpo::value<bool>()->default_value(false)->implicit_value(true),"Include Standard Deviation");
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
	int cloudband = Opts.varmap.count("cloudband") ? Opts.varmap["cloudband"].as<int>() : 0;
	bool stddev = Opts.varmap["stddev"].as<bool>();
	// Main algorithm
	// vector file
	GeoImage<T> vfile(vfilename);
	GeoStats<T> vStats(vfile[0]);
	T maxID = vStats.Max();
	GeoRaster<T> *cloud;

	int id_lo = 1;
	int id_hi = maxID;
	if (Opts.varmap.count("id")) {
		int id = Opts.varmap["id"].as<int>();
		id_lo = id;
		id_hi = id;
	}
	//cout << "Calculating Stats for IDs " << id_lo << " - " << id_hi << endl;

	int wsm(8);
	int w(12);
	cout << setw(w) << "Datafile" << setw(wsm) << "FID" << setw(wsm) << "NumPix"; // << setw(wsm) << "Year" << setw(wsm) << "DOY";
	for (int b=0;b<Opts.Bands().size();b++) {
		//cout << setw(wsm) << ("b"+to_string(Opts.Bands(b))+"_NumP")
		cout << setw(w) << ("b" + to_string(Opts.Bands(b)) + "_Mean");
		if (stddev) cout << setw(w) << ("b" + to_string(Opts.Bands(b)) + "_StdDev");
	}
	cout << endl;

	// Loop through input files
	for (int i=0;i<Opts.InputFiles().size(); i++) {
		GeoImage<T> Image(Opts.InputFiles(i),Opts.Gain(),Opts.Offset());
		vector< GeoMask<T> > Masks;
		vector< map<T,POI> > POIs;
		// data specific hardcodes
		Image[5].SetGainOffset(0.01);
		Image[6].SetGainOffset(0.01);
		//string year= Image.Basename().substr(13,4);
		//string doy = Image.Basename().substr(33,3);

		// Pointer to cloud band
		cloud = (cloudband > 0) ? &(Image[cloudband-1]) : NULL;

		// Map Pixels of Interest
		for (int b=0;b<Opts.Bands().size(); b++) {
			// Set NoData value, if there is one, and map POIs for remaining good pixels
			Masks.push_back( GeoMask<T>(cloud,0,true) );
			if (Opts.NoData()) {
				Image[Opts.Bands(b)-1].SetNoData(Opts.NoDataValue());
				Masks[b].AddMask(&Image[Opts.Bands(b)-1],Opts.NoDataValue(),false);
			}
			POIs.push_back( vfile[0].MapPOIs(Masks[b]) );

			//Masks[b].Write("mask"+to_string(b));
			/*map<T,POI>::const_iterator imap;
			for (imap=POIs[b].begin();imap!=POIs[b].end();imap++) {
				int num = imap->second.NumPixels();
				if (num) cout << imap->first << ": " << num << endl;
			}*/
		}

		// Loop through POI regions
		for (int id=id_lo;id<=id_hi;id++) {
			vector< GeoStats<T> > Stats;
			int numpixels = 0;
			for (int b=0;b<Opts.Bands().size(); b++) {
				map<T,POI>::iterator iPOI = POIs[b].find(id);
				if (iPOI != POIs[b].end())
					Stats.push_back( GeoStats<T>(Image[Opts.Bands(b)-1], iPOI->second, !stddev ));
				else Stats.push_back( GeoStats<T>() );
				numpixels += Stats[b].NumPixels();
				//cout << id << " NumPixels = " << numpixels << endl;
			}
			if (numpixels > 0) {
				// ID, Year, doy
				cout << setw(w) << Image.Basename() << setw(wsm) << id << setw(wsm) << Stats[0].NumPixels(); // << setw(wsm) << year << setw(wsm) << doy;
				for (int b=0;b<Opts.Bands().size();b++) {
					//cout << setw(wsm) << Stats[b].NumPixels() 
					cout << setw(w) << Stats[b].Mean();
					if (stddev) cout << setw(w) << Stats[b].StdDev();
				}
				cout << endl;
			}
		}
	}
	
	return 0;
}
