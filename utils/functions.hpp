#pragma once

#include <math/types.hpp>
#include <boost/math/special_functions/sinc.hpp>

namespace spurt
{

constexpr double PI = 3.14159265358979323846264338;
constexpr double TWO_PI=6.28318530717958647688;


inline double gauss(double r) {
    return exp(-r*r/2)/sqrt(TWO_PI);
}

inline double sqr(double x) {
    return x*x;
}
inline double cube(double x) {
    return x*x*x;
}

namespace functions3d {
    struct ABC_field {
        ABC_field(double a=sqrt(3), double b=sqrt(2), double c=1) : A(a), B(b), C(c) {}
    
        spurt::vec3 evaluate(const spurt::vec3& x) const {
            spurt::vec3 dxdt;
            dxdt[0] = A*sin(x[2]) + C*cos(x[1]);
            dxdt[1] = B*sin(x[0]) + A*cos(x[2]);
            dxdt[2] = C*sin(x[1]) + B*cos(x[0]);
            return dxdt;
        }

        double A, B, C;
    };

    inline double sphere(const spurt::vec3& p) {
        return norm(p);
    }

    inline double ball(const spurt::vec3& p) {
        double r = norm(p);
        if (r>1) return 0;
        else return gauss(r-1);
    }

    inline double ellipsoid(const spurt::vec3& p) {
        return norm(vec3(sqrt(3)*p[0], sqrt(2)*p[1], sqrt(5)*p[2]));
    }

    inline double torus(const spurt::vec3& p) {
        double theta = atan2(p[1], p[0]);
        double px = cos(theta);
        double py = sin(theta);
        return norm(vec3(p[0]-px, p[1]-py, p[2]*p[2]));
    }

    inline double helix(const spurt::vec3& p) {
        double theta = atan2(p[1], p[0]);
        double px = cos(theta);
        double py = sin(theta);
        double pz = 0.5*theta;
        return norm(vec3(p[0]-px, p[1]-py, p[2]-pz));
    }

    struct marshner_lobb{
        marshner_lobb(double alpha=0.25, double f_M=6) 
            : _M_alpha(alpha), _M_f_M(f_M) {}
    
        double rho_r(double r) const {
            return cos(2*PI*_M_f_M*cos(PI*r/2.));
        }
    
        double operator()(double x, double y, double z) const {
            return ((1-sin(PI*z/2)) + _M_alpha*(1 + rho_r(sqrt(x*x + y*y))))/(2*(1+_M_alpha));
        }
    
        double operator()(const spurt::vec3& p) const {
            return (*this)(p[0], p[1], p[2]);
        }
    
        double _M_alpha;
        double _M_f_M;
    };

} // functions3d

namespace functions2d {

    inline double mathworks1(const spurt::vec2& p) {
        const double& x = p[0];
        const double& y = p[1];
        return 3.*exp(-sqr(y+1)-sqr(x))*sqr(x-1)-exp(-sqr(x-1)-sqr(y))/3 + exp(-sqr(x)-sqr(y))*(10*cube(x)-2*x+10*pow(y,5));
    }

    inline double sinc(const spurt::vec2& p) {
        double r = sqrt(p[0]*p[0] + p[1]*p[1]);
        if (r==0) return 1;
        else return sin(r)/r;
    }

    inline double sincos(const spurt::vec2& p) {
        return sin(TWO_PI*p[0])*cos(TWO_PI*p[1]);
    }

    inline std::complex<double> complex_sinc(const spurt::vec2& p) {
        std::complex<double> z(p[0], p[1]);
        return boost::math::sinc_pi(z);
    }

    inline double z_marschner_lobb(const spurt::vec2& p) {
        double r = norm(p);
        return 2/PI*asin(cos(12*PI*cos(PI*r/2))/4);
    }

    inline double circle(const spurt::vec2& p) {
        return fabs(norm(p) - 1);
    }

    inline double radial(const spurt::vec2& p) {
        return norm(p);
    }

    inline double disk(const spurt::vec2& p) {
        double r = norm(p);
        if (r < 1) return 0;
        else return gauss(r - 1);
    }


} // functions2d
    
    
    
}