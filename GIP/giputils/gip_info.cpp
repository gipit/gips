/*!
 * \ingroup gip_test
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <iostream>
#include <gip/GeoImage.h>
#include <gip/GeoAlgorithms.h>
#include <gip/GeoRasterIO.h>

#include <gip/Options.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	using namespace boost::program_options;
	options_description localopts("Info Options");
	localopts.add_options()
		("stats", value<bool>()->implicit_value(true)->default_value(false), "Calculate Stats")
		("updatetest", value<bool>()->implicit_value(true)->default_value(false), "Update test")
		;
	gip::Options Opts(ac,av,localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	bool stats = Opts.varmap["stats"].as<bool>();

	for (unsigned int i=0;i<Opts.InputFiles().size();i++) {
		try {
			GeoImage img(Opts.InputFile(i), false);
			if (stats) img.ComputeStats();
			cout << img.Info(true, stats) << endl;
		} catch (exception const& err) {
			cout << "Error " << err.what() << endl;
		}
	}

	return 0;
}
