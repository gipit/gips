#ifndef GIP_OPTIONS_H
#define GIP_OPTIONS_H

#include <string>

#include <boost/program_options.hpp>
//#include <boost/program_options/errors.hpp>
#include <boost/filesystem.hpp>

#include <gip/Colors.h>

namespace gip {
	class Options {
	public:
        /*using std::vector;
        using std::string;
		using boost::program_options::value;
        using boost::program_options::options_description;
        using boost::program_options::positional_options_description;*/
		//! Default constructor sets up default options
		Options() {
			//SetupInput();
			//SetupOutput();
		}
		//! Constructor to add parse command line and add local options
		Options(int ac, char** av, boost::program_options::options_description localopts);

		//! Constructor to read from config file
		Options(std::string filename);

		//! Copy constructor
		Options(const Options& opts)
			: _Inputs(opts._Inputs), _Colors(opts._Colors),
			  _Gain(opts._Gain), _Offset(opts._Offset),
			  _NoData(opts._NoData), _NoDataValue(opts._NoDataValue),
			  _Outputs(opts._Outputs), _Prefix(opts._Prefix), _Suffix(opts._Suffix),
			  _Meta(opts._Meta) {}
		~Options() {}

		//! Copy assignment operator
		Options& operator=(const Options& opts) {
			if (this == &opts) return *this;
			_Inputs = opts._Inputs;
			//_Bands = opts._Bands;
			_Colors = opts._Colors;
			_Gain = opts._Gain;
			_Offset = opts._Offset;
			_NoData = opts._NoData;
			_NoDataValue = opts._NoDataValue;
			_Outputs = opts._Outputs;
			_Prefix = opts._Prefix;
            _Suffix = opts._Suffix;
			return *this;
		}

        //! Non-const version
		//Options& operator=(Options& opts) {
        //    return const_cast<Options&>(static_cast<const Options&>(*this));
		//}

		// Get Commands

		// Inputs
		std::string InputFile() const { return _Inputs[0]; }
		std::string InputFile(const int i) const { return _Inputs[i]; }
		std::vector<std::string> InputFiles() const { return _Inputs; }
		//std::vector<unsigned int> Bands() const { return _Bands; }
		//int Bands(unsigned int i) const { return _Bands[i]; }
		Colors GetColors() const { return _Colors; }
		std::vector<float> Gain() const { return _Gain; }
		float Gain(unsigned int i) const { return _Gain[i]; }
		std::vector<float> Offset() const { return _Offset; }
		float Offset(unsigned int i) const { return _Offset[i]; }
		bool NoData() const { return _NoData; }
		double NoDataValue() const { return _NoDataValue; }

		// Outputs
		std::string OutputFile() const { return _Outputs.empty() ? "" : _Outputs[0]; }
		std::string OutputFile(const int i) const { return _Outputs[i]; }
		std::vector<std::string> OutputFiles() const { return _Outputs; }
		std::vector<std::string> Meta() const { return _Meta; }
		std::string Meta(unsigned int i) const { return _Meta[i]; }
		std::string Prefix() const { return _Prefix; }
        std::string Suffix() const { return _Suffix; }
		//bool PAM() const { return _PAM; }

		// \name Global Options (static properties)
		//! Get Config directory
		static std::string ConfigDir() { return _ConfigDir.string(); }
		//! Set Config directory
		static void SetConfigDir(std::string dir) { _ConfigDir = dir; }
		//! Default format when creating new files
		static std::string DefaultFormat() { return _DefaultFormat; }
		//! Set default format when creating new files
		static void SetDefaultFormat(std::string str) { _DefaultFormat = str; }
		//! Default chunk size when chunking an image
		static float ChunkSize() { return _ChunkSize; }
		//! Set chunk size, used when chunking an image
		static void SetChunkSize(float sz) { _ChunkSize = sz; }
		//! Get verbose level
		static int Verbose() { return _Verbose; }
		//! Set verbose level
		static void SetVerbose(int v) { _Verbose = v; }

		bool Exit() const { return _Exit; }

		friend std::ostream& operator<<(std::ostream&, const Options&);

		//! Variable map of options
		boost::program_options::variables_map varmap;
	private:
		// Static options
		//! Configuration file directory
		static boost::filesystem::path _ConfigDir;
		//! Default format
		static std::string _DefaultFormat;
		//! Chunk size used when chunking up an image
		static float _ChunkSize;
		//! Verbosity level
		static int _Verbose;

		//! Exit flag (to indicate program should exit after parsing cmdline)
		bool _Exit;

		// Setting up and Parsing options
		//! Setup input options
		static boost::program_options::options_description SetupInput();
		//! Setup output options
		static boost::program_options::options_description SetupOutput();

		//! Parse the variable map (varmap)
		void Parse(boost::program_options::options_description);

		// Input Options
		//! Input filenames
		std::vector<std::string> _Inputs;
		//! Vector of bands of interest
		//std::vector<unsigned int> _Bands;
		//! Band Colors - Referenced to _Bands array
		Colors _Colors;
		//! Gain to apply to input
		std::vector<float> _Gain;
		//! Offset to apply to input
		std::vector<float> _Offset;
		//! Nodata value
		bool _NoData;
		double _NoDataValue;

		// Output Options
		//! Output filenames
		std::vector<std::string> _Outputs;
		//! Location of output
		std::string _Prefix;
        std::string _Suffix;
		//! Meta data of interest (to include on output)
		std::vector<std::string> _Meta;
		//! PAM support (Persistent Auxilary Metadata)
		//bool _PAM;

		// Other
		//! Verbosity
		//int _Verbose;
		//! Parse Python dictionary
		/*int ParseDict(string str) { //boost::python::dict dict) {
			std::cout << "HELLO " << str << std::endl;
			return 1;
		}*/
		//int Verbose() const { return _Verbose; }
	};
}
#endif
