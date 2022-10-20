#ifndef __XAVIER_MACROS_HPP__
#define __XAVIER_MACROS_HPP__

#include <iostream>

namespace spurt {
namespace map_debug {
extern int verbose_level;
}

#define WARNING_MACRO(a, txt) {             \
        if (map_debug::verbose_level > a)           \
        {                                           \
            std::cout << txt << std::flush;         \
        }                                           \
    }

}

#endif

