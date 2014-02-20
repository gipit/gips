/*!
 * \ingroup gip_test
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "Options.h"
#include "gip.h"
#include "gip_GeoProcess.h"

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Setup program options and parse command line
	bpo::options_description localopts("Test Options");
	localopts.add_options()
		("copy", bpo::value<bool>()->implicit_value(true)->default_value(false), "Copy and assignment test")
		("outtest", bpo::value<bool>()->implicit_value(true)->default_value(false), "Output test");
	Options Opts(localopts);
	try {
		Opts.ParseCmdLine(ac,av);
	} catch (std::exception& err) {
		cout << "Error: " << err.what() << endl;
		return 1;
	}
	if (Opts.Exit()) return 0;
	// End parse command line
	GDALAllRegister();
	typedef float T;

	// Main algorithm
	//cout << "Input Files: " << endl; for (int i=0;i<Opts.InputFiles().size();i++) cout << "\t" << Opts.InputFiles(i) << endl;
	//cout << "Output Files: " << endl; for (int i=0;i<Opts.OutputFiles().size();i++) cout << "\t" << Opts.OutputFiles(i) << endl;
	//cout << "Bands: "; for (int i=0;i<Opts.Bands().size();i++) cout << " " << Opts.Bands(i); cout << endl;

	for(unsigned int iFiles=0;iFiles<Opts.InputFiles().size();iFiles++) {
		GeoImage GeoImageIn(Opts.InputFiles(iFiles));
		GeoImageIn.RemoveBand(1);
		// Output test
		if (Opts.varmap["copy"].as<bool>()) {
			/*GeoImage GeoImageOut(Opts.OutputFiles(iFiles),*GeoImageIn[0].GetGeoImage(),GeoImageIn.NumBands());
			// Copy bands
			for (int band=0;band<GeoImageOut.NumBands(); band++) {
				cout << "band = " << band << endl;
				GeoProcess<T> img(GeoImageOut[band]);
				img.Copy(GeoImageIn[band]);
			}*/
			GeoImage GeoImageOut = CreateCopy(GeoImageIn, Opts.OutputFiles(iFiles));
			cout << "Image In" << endl << GeoImageIn.Info() << endl;
			cout << "Image Out" << endl << GeoImageOut.Info() << endl;
		}
	}
	return 0;
}
