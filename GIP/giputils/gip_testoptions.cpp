/*!
 * \ingroup gip_test
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <iostream>
#include <gip/Options.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	bpo::options_description localopts("Test Options");
	localopts.add_options()
		("test", bpo::value<bool>()->implicit_value(true)->default_value(false), "Update test");
	gip::Options Opts(ac, av, localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	cout << Opts << endl;

	cout << "Config file test #1" << endl;
	Options config1(Options::ConfigDir() + "MCD43A4.gipcfg");
	cout << "config1: " << endl << config1 << endl;

	cout << "Config file test #2" << endl;
	Options config2(Options::ConfigDir() + "MCD43A4.gipcfg");
	cout << "config2: " << endl << config2 << endl;

	return 0;
}
