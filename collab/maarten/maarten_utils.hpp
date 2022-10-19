#ifndef __XAVIER_MAARTENS_UTILS_HPP__
#define __XAVIER_MAARTENS_UTILS_HPP__

// STL
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>
// boost
#include <boost/shared_ptr.hpp>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
// xavier
#include <math/RBF.hpp>
#include <math/RBFbasis.hpp>
#include <image/nrrd_wrapper.hpp>
#include <format/format.hpp>
#include <maarten/data_IO.hpp>
#include <maarten/typedefs.hpp>

namespace xavier { namespace maarten {

inline double deg2rad(double deg) {
    static const double k = M_PI/180.;
    return deg*k;
}

// compute km distance between two points expressed in longitude / latitude
// formula found (by Kim) here: 
// http://bluemm.blogspot.com/2007/01/excel-formula-to-calculate-distance.html
// and successfully tested on Hui Huang's data
inline double distance_in_km(const nvis::vec2& lola0, const nvis::vec2& lola1) {
    static const double MEAN_EARTH_RADIUS = 6371;
    double theta0 = deg2rad(lola0[0]);
    double theta1 = deg2rad(lola1[0]);
    double phi0 = deg2rad(90-lola0[1]);
    double phi1 = deg2rad(90-lola1[1]);
    return acos(cos(phi0)*cos(phi1) + sin(phi0)*sin(phi1)*cos(theta0-theta1))*
           MEAN_EARTH_RADIUS;
}

inline double km_per_angle_ratio(const nvis::vec2& lola, 
                                 const nvis::vec2& direction,
                                 double eps=0.0001) {
    // note: selecting a value of epsilon less than 1.0e-4 leads to 
    //       wrong results
    const nvis::vec2& from = lola;
    nvis::vec2 to = lola + eps*direction;
    double km = distance_in_km(from, to);
    return km/eps;
}

class progress_message {
    std::string        _what;
    size_t             _size;
    
    mutable std::ostringstream _os;
    
public:
    progress_message(size_t size, std::string what = "")
    : _size(size), _what(what) {}
    
    std::string operator()(size_t n, double elapsed=0) const {
        _os.clear();
        _os.str("");
        _os << "\rCompleted " << 100.*(float)n/float(_size)
        << "\% of " << _size;
        if (_what.size()) {
            _os << " " << _what;
        }
        if (elapsed > 0) {
            _os << " in " << elapsed << " seconds (" << (float)n/elapsed << " Hz)";
        }
        _os << "                       \r";
        return _os.str();
    }
    
    void reset(const std::string& what = "") {
        if (what.size()) {
            _what = what;
        }
    }
};

} // maarten

} // xavier
#endif