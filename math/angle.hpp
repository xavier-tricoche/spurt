#ifndef __ANGLE_HPP__
#define __ANGLE_HPP__

#include <math/fixed_vector.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace {
inline double sign(const double& x)
{
    return (x >= 0 ? 1 : -1);
}
}

namespace xavier {
const double TWO_PI = 2*M_PI;

inline double signed_angle(const nvis::vec2& x)
{
    double l = nvis::norm(x);
    if (!l) {
        return 0;
    }
    
    double theta = acos(x[0] / l);
    
    if (theta > M_PI) {
        throw std::runtime_error("invalid angle value\n");
    }
    
    if (x[1] < 0) {
        theta *= -1;
    }
    return theta;
}

inline double positive_angle(const nvis::vec2& x)
{
    double theta = signed_angle(x);
    return theta < 0 ? theta + 2.*M_PI : theta;
}

inline double positive_modulo(double angle, double period = 2*M_PI)
{
    return angle >= 0 ? fmod(angle, period) : period + fmod(angle, period);
}

inline double modulo(double angle, double period = 2*M_PI)
{
    return fmod(angle, period);
}

inline double signed_angle(const nvis::vec2& v0, const nvis::vec2& v1)
{
    if (!nvis::norm(v0) || !nvis::norm(v1)) {
        return 0;
    }
    nvis::vec2 w0 = v0 / nvis::norm(v0);
    nvis::vec2 w1 = v1 / nvis::norm(v1);
    double dot = nvis::inner(w0, w1);
    if (dot <= -1) {
        return M_PI;
    } else if (dot >= 1) {
        return 0;
    }
    double theta = acos(dot);
    if (std::isnan(theta)) {
        std::ostringstream os;
        os << "signed angle return nan. v0 = " << v0 << ", v1 = " << v1 << ", w0 = " << w0 << ", w1 = " << w1 << std::endl;
        std::cerr << os.str();
    }
    double det = w0[0] * w1[1] - w0[1] * w1[0];
    return sign(det)*theta;
}

template<typename Vector_>
inline double signed_angle(const Vector_& v0, const Vector_& v1) {
    return signed_angle(nvis::vec2(v0[0], v0[1]), nvis::vec2(v1[0], v1[1]));
}

template<int N>
inline double unsigned_angle(const nvis::fixed_vector<double, N>& v0, const nvis::fixed_vector<double, N>& v1)
{
    typedef typename nvis::fixed_vector<double, N>  vec_type;
    
    if (!nvis::norm(v0) || !nvis::norm(v1)) {
        return 0;
    }
    vec_type w0 = v0 / nvis::norm(v0);
    vec_type w1 = v1 / nvis::norm(v1);
    return acos(nvis::inner(w0, w1));
}

inline double to_degrees(double rad)
{
    return rad/M_PI*180.;
}


}

#endif




