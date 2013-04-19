/*!
 * \ingroup ags_utilities
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

// STL
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <agsOptions.h>
#include <agsGDAL.h>

using namespace std;
using namespace ags;

int main (int ac, char* av[]) {
	// Setup program options and parse command line
	boost::program_options::options_description localopts("Test Options");
	localopts.add_options()
		("min", "Crop to minimum footprint")
		("max", "Crop to maximum footprint")
		("test,t", bpo::value<int>(), "Print values only");
	ags::agsOptions Opts(localopts);
	try {
		Opts.ParseCmdLine(ac,av);
	} catch (std::exception& err) {
		cout << "Error: " << err.what() << endl;
		return 1;
	}
	// End parse command line
	cout << "Input Files: "; for (int i=0;i<Opts.InputFiles().size();i++) cout << Opts.InputFiles(i) << " "; cout << endl;

	GDALAllRegister();

	// Main algorithm
	agsGDALDataset* img;
	for(int iFiles=0;iFiles<Opts.InputFiles().size();iFiles++) {
		img = (agsGDALDataset*)GDALOpen(Opts.InputFiles(iFiles).c_str(), GA_ReadOnly );
		cout << *(img->GetFileList()) << endl;
		cout << "\tTop Left Corner: " << setprecision(7) << bgeo::dsv(img->TopLeft()) << endl;
		cout << "\tLower Right Corner: " << bgeo::dsv(img->LowerRight()) << endl;
	}
}
