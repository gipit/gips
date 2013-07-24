#ifndef GIPCIMG_H
#define GIPCIMG_H

#define cimg_debug 0
#define cimg_display 0

#include <CImg.h>

namespace gip {

    using cimg_library::CImg;

    template<typename T> inline void cimg_printclasses(CImg<T> img) {
        for (int i=0; i<img.height(); i++) {
            std::cout << "\tClass " << i+1 << ": ";
            cimg_forX(img,x) std::cout << img(x,i) << "  ";
            std::cout << std::endl;
        }
        return;
    }

}

#endif
