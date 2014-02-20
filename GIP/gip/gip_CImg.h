#ifndef GIPCIMG_H
#define GIPCIMG_H

//#define cimg_debug 0
#define cimg_verbosity 1
#define cimg_display 0

#include <CImg.h>

namespace gip {

    using cimg_library::CImg;

    template<typename T> inline void cimg_printclasses(cimg_library::CImg<T> img, std::string prefix="Class") {
        for (int i=0; i<img.height(); i++) {
            std::cout << "\t" << prefix << " " << i+1 << ": ";
            cimg_forX(img,x) std::cout << img(x,i) << "  ";
            std::cout << std::endl;
        }
        return;
    }

    template<typename T> inline void cimg_printstats(cimg_library::CImg<T> img, std::string note="") {
        if (note != "") std::cout << note << "  ";
        CImg<float> stats = img.get_stats();
        std::cout << "Min/Max = " << stats(0) << ", " << stats(1)
            << " Mean/StdDev = " << stats(2) << " +/- " << stats(3)
            << " Sum = " << img.sum()
            << " NumPixels = " << img.size() << std::endl;
    }

    template<typename T> inline void cimg_print(cimg_library::CImg<T> img, std::string note="") {
        std::cout << note << ": ";
        //" (" << img.width() << ", " << img.height() << ", " << img.depth() << ", " << img.spectrum() << "): ";
        cimg_for(img,ptr,T) std::cout << *ptr << "  ";
        std::cout << std::endl;
    }

}

#endif
