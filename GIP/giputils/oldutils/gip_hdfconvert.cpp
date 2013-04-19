/*!
 * \ingroup ags_utilities
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <gip_Options.h>
#include <gip_GeoImage.h>
//#include <boost/any.hpp>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main (int ac, char* av[]) {
	// Setup program options and parse command line
	bpo::options_description localopts("Test Options");
	localopts.add_options()
		("info", bpo::value<bool>()->default_value(false)->implicit_value(true), "Print info on input files and exit")
		("meta,m", bpo::value<vector<string> >()->multitoken(), "Metadata keys to copy from input file");
	gip::Options Opts(localopts);
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
	//typedef boost::any T;
	typedef float T;

	// Options
	string prefix = Opts.varmap["prefix"].as<string>();
	bool info = Opts.varmap["info"].as<bool>();
	vector<string> meta;
	if (Opts.varmap.count("meta")) meta = Opts.varmap["meta"].as<vector<string> >();

	// Main algorithm
	// Loop through input files
	for(unsigned int iFiles=0;iFiles<Opts.InputFiles().size();iFiles++) {
		// Input Image
		gip::GeoImage<T> Images(Opts.InputFiles(iFiles).c_str(),Opts.Gain(),Opts.Offset());

		// Subdatasets 
		vector<string> subnames = Images.GetMetaGroup("SUBDATASETS","_NAME=");
		if (info) {
			cout << "Subdatasets for " << Images.Basename() << ":" << endl;
			for (unsigned int i=0; i<subnames.size(); i++) 
				cout << "\t" << i << ": " << subnames[i].substr(subnames[i].rfind(':')+1) << endl;
			// Skip rest of code, move to next input image
			continue;
		}
		// Create Output
		string fnameout = prefix + bfs::basename(bfs::path(Opts.InputFiles(iFiles)).leaf());
		GDALDataType datatype = (Opts.Gain() == 1.0 && Opts.Offset() == 0.0) ? 
			Images[0].GetGDALRasterBand()->GetRasterDataType() : gip::GDALType(&typeid(T));
		gip::GeoImage<T> OutputImages(Images[0].XSize(),Images[0].YSize(),Opts.Bands().size(),fnameout,Opts.Format(),datatype);

		// Copy metadata
		OutputImages.CopyCoordinateSystem(Images[0]);
		OutputImages.CopyMetadata(Images[0], meta);
		
		// Copy bands
		for (unsigned int b=0;b<Opts.Bands().size(); b++) {
			cout << "Copying band " << b << endl;
			string bandname = subnames[Opts.Bands(b)-1];
			OutputImages[b].Copy(Images[b]);
			// Set Description
			OutputImages[b].SetDescription(bandname.substr(bandname.rfind(":")+1));
			//OutputImages[b].SetDescription("Band "+to_string(b), bandname.substr(bandname.rfind(":")+1));
		}
		cout << "Created " << fnameout << endl;
		
	}
	return 0;
}
