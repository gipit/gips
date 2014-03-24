#ifndef FUNCTION_H
#define FUNCTION_H

#include <boost/function.hpp>
#include <gip/gip_CImg.h>

namespace gip {
	using cimg_library::CImg;
	
	class Function {

	public:
		Function() {}
		~Function() {}

		//! Unary member function
		boost::function< CImg<double>& (CImg<double>&) > F;
	private:


	}; // class Function

} // namespace GIP

#endif