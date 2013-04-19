/*
 * Options.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: mhanson
 */

#include <vector>
#include <fstream>
//#include <iostream>

#include <gdal/gdal_priv.h>

#include <gip/Options.h>
#include <gip/Utils.h>

namespace gip {
    using std::string;
    using std::vector;
    using boost::program_options::value;
    using boost::program_options::options_description;
    using boost::program_options::positional_options_description;

	boost::filesystem::path Options::_ConfigDir("/usr/share/gip/");
	string Options::_DefaultFormat("GTiff");
	float Options::_ChunkSize(32.0);
	int Options::_Verbose(1);

	Options::Options(string filename) {
        options_description opts;
        opts.add(SetupInput());
        opts.add(SetupOutput());
        // TODO check for file existence
        std::ifstream in(filename.c_str());
        if (in.is_open()) {
            if (Verbose()) std::cout << "Reading configuration file " << filename << std::endl;
            boost::program_options::store(boost::program_options::parse_config_file(in,opts), varmap);
            notify(varmap);
            Parse(opts);
        }
    }

    Options::Options(int ac, char** av, boost::program_options::options_description localopts)
                : _Prefix(""), _Suffix("") {

        options_description opts;
        opts.add(SetupInput());
        opts.add(SetupOutput());
        opts.add(localopts);

        positional_options_description opts_pos;
        opts_pos.add("input",-1);

        //Parse command line
        options_description opts_cmdline("Command Line Options");
        opts_cmdline.add_options()
            ("verbose,v", value<int>()->implicit_value(1), "Print additional information")
            ("help,h", "Print help message");
        opts.add(opts_cmdline);
        boost::program_options::store(
            boost::program_options::command_line_parser(ac,av).options(opts).positional(opts_pos).run(), varmap);
        notify(varmap);
        Parse(opts);
        // Streaming input support - should only be on command line!!
        /*if (_Inputs.empty() && !_Exit) {
            string input;
            while(std::cin) {
                std::getline(std::cin,input);
                if (!input.empty()) _Inputs.push_back(input); else break;
            }
        }*/
        //if (_Inputs.empty() && !_Exit) throw bpo::error("Input file(s) missing");
    }

	options_description Options::SetupInput() {
		options_description opts_input("Input Options");
		opts_input.add_options()
			("input,i", value<vector<string> >()->multitoken(), "Input filename(s)")
			//("bands,b", value<string>(), "Input band numbers interest")
			("colors", value<string>(), "Color Bandnumbers (Blue Green Red NIR SWIR1 SWIR2 LWIR)")
			("gain", value<string>(), "Gain to apply to input data on reads")
			("offset", value<string>(), "Offset to apply to input data on reads")
			("nodata", value<double>(), "No Data value")
			("configdir", value<string>(), "Directory of configuration files")
			("chunksz", value<float>()->default_value(32.0), "Size of chunks (in MB)")
			("gdaldebug", value<bool>()->default_value(false)->implicit_value(true), "GDAL Debug (CPL_DEBUG)")
			;
        return opts_input;
		//_opts->add(opts_input);
		// Input as positional argument (for command line)
		//_opts_pos->add("input",-1);
	}

	options_description Options::SetupOutput() {
		options_description opts_output("Output Options");
		opts_output.add_options()
			("output,o", value<vector<string> >()->multitoken(), "Output filename(s)")
			("prefix,p", value<string>()->default_value(""), "Prefix/Path to prepend to output files")
			("suffix,s", value<string>()->default_value(""), "Suffix to prepend to output files")
			("nopam", value<bool>()->default_value(false)->implicit_value(true),"Disable PAM in GDAL output (auxilary metadata)")
			("meta,m", value<string>(), "Meta data keys to preserve on output")
			("format,f", value<string>()->default_value("GTiff"), "File format of outputs");
        return opts_output;
		//_opts->add(opts_output);
	}

	void Options::Parse(options_description opts) {
		_Exit = false;
		//_Verbose = varmap.count("verbose") ? varmap["verbose"].as<int>() : 0;
		if (varmap.count("help")) { std::cout << opts << std::endl; _Exit = true; }

		// Static options
		if (varmap.count("configdir")) Options::SetConfigDir(varmap["configdir"].as<string>());
		if (varmap.count("chunksz")) Options::SetChunkSize(varmap["chunksz"].as<float>());
		if (varmap.count("verbose")) Options::SetVerbose(varmap["verbose"].as<int>());
		if (varmap.count("format")) Options::SetDefaultFormat(varmap["format"].as<string>());

		// GDAL Config Options
		bool pam = varmap["nopam"].as<bool>() ? false : true;
		CPLSetConfigOption("GDAL_PAM_ENABLED",pam ? "TRUE" : "FALSE");
		bool gdaldebug = varmap["gdaldebug"].as<bool>() ? true : false;
		CPLSetConfigOption("CPL_DEBUG", gdaldebug ? "ON" : "OFF");
		//CPLPushErrorHandler(CPLQuietErrorHandler);

		if (varmap.count("input")) _Inputs = varmap["input"].as<vector<string> >();
        if (_Inputs.empty()) { std::cout << opts << std::endl; _Exit = true; }

		if (varmap.count("gain")) {
			 vector<string> strGain = Split(varmap["gain"].as<string>(), " ,");
			 for (vector<string>::const_iterator it=strGain.begin(); it!=strGain.end(); it++)
				 _Gain.push_back( atof(it->c_str()) );
		}

		if (varmap.count("offset")) {
			vector<string> strOffset = Split(varmap["offset"].as<string>(), " ,");
			for (vector<string>::const_iterator it=strOffset.begin(); it!=strOffset.end(); it++)
				_Offset.push_back(atof(it->c_str()));
		}
		_NoData = varmap.count("nodata") ? true : false;
		if (_NoData) _NoDataValue = varmap["nodata"].as<double>();

		// Get bands, ranges allowed e.g. 1-7 9
		/*if (varmap.count("bands")) {
            std::vector<string> str = Split(varmap["bands"].as<string>(), " ,");
            std::vector<string>::const_iterator iv;
            size_t loc;
            for (iv=str.begin();iv!=str.end();iv++) {
                loc = iv->find("-");
                if (loc==string::npos)
                    _Bands.push_back( atoi(iv->c_str()) );
                else {
                    int b1 = atoi(iv->substr(0,loc).c_str());
                    int b2 = atoi(iv->substr(loc+1).c_str());
                    for (int i=b1;i<=b2;i++) _Bands.push_back(i);
                }
            }
		}*/

        // TODO - FIX THIS!
		/*if (varmap.count("colors")) {
			vector<string> colors = Split(varmap["colors"].as<string>(), " ,");
			vector<string>::const_iterator iv;
			int col(1);
			for (iv=colors.begin();iv!=colors.end();iv++) _Colors.SetColor(col++,atoi(iv->c_str()));
		}*/
        //std::cout << "Prefix_" << std::endl;
		_Prefix = varmap.count("prefix") ? varmap["prefix"].as<string>() : "";
        _Suffix = varmap.count("suffix") ? varmap["suffix"].as<string>() : "";
        if (_Suffix == "") _Suffix = ".out";

		if (varmap.count("output")) _Outputs = varmap["output"].as<vector<string> >();
		// If no output files give, create new ones based on input filenames

		if (_Outputs.empty()) {
			for (unsigned int i=0;i<_Inputs.size();i++) {
				boost::filesystem::path out(_Inputs[i]);
				string root = out.branch_path().string();
				if (root != "") root = root + '/';
				_Outputs.push_back( root + _Prefix + boost::filesystem::basename(out)+_Suffix );
			}
		}
        //std::cout << "Parse4" << std::endl;
		if (varmap.count("meta")) _Meta = Split(varmap["meta"].as<string>(), " ,");

	}

	std::ostream& operator<<(std::ostream& out, const Options& opts) {
		using std::endl;
		out << "GIP Options" << endl;

		unsigned int i;
		out << "Input/Output Files:" <<endl;
		for (i=0;i<opts._Inputs.size();i++) out << "\t" << opts._Inputs[i] << " -> " << opts._Outputs[i] << endl;
		out << endl;

		//out << "Bands:";
		//for (i=0;i<opts._Bands.size();i++) out << "\t" << opts._Bands[i]; out << endl;
		//out << "Color: ";
		//for (i=0;i<opts._Bands.size();i++) out << "\t" << opts._Colors[i+1]; out << endl;
		out << "Gain: ";
		for (i=0;i<opts._Gain.size();i++) out << "\t" << opts._Gain[i]; out << endl;
		out << "Offset: ";
		for (i=0;i<opts._Offset.size();i++) out << "\t" << opts._Offset[i]; out << endl;
		out << "NoData: " << opts._NoData << " " << opts._NoDataValue << endl;
		out << endl;

		out << "Prefix: " << opts._Prefix << endl;
		out << "Suffix: " << opts._Suffix << endl;
		//out << "PAM: " << opts._PAM << endl;
		out << "Meta: " << endl;
		for (i=0;i<opts._Meta.size();i++) out << "\t" << opts._Meta[i] << endl;

		// Static options
		//out << "Format: " << opts._Format << endl;

		//out << "Local Opts:" << endl;
		//out << "Test: " << opts.varmap["test"].as<bool>() << endl;
		return out;
	}
}
