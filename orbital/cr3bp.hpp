#ifndef __cr3bp_hpp
#define __cr3bp_hpp

#include <cmath>
#include <vector>
#include <complex>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/bezier_spline.hpp>

#include "rhs_base.hpp"

namespace
{
double p_inf    = std::numeric_limits<double>::max();
double m_inf    = std::numeric_limits<double>::min();
};

const double Jupiter_Europa_mu = 2.528017705e-5;
const double Earth_Moon_mu = 0.012150571430596;

class cr3bp: public rhs_base<6>
{
public:

    // state variables
    // 0: x, 1: y, 2: dxdt, 4: dydt

    typedef nvis::fixed_vector<double, 6>     vec6;
    typedef nvis::fixed_vector<double, 7>     vec7;
    typedef nvis::bounding_box<vec6>          bbox6;

    double yd(double x, double xd) const {
        double ysq = 2*(1-m_mu)/std::abs(x+m_mu) + 2*m_mu/std::abs(x-1+m_mu) + x*x - xd*xd - m_C;
        if (ysq < 0) {
            throw std::runtime_error("invalid position exception");
        }
        else return sqrt(ysq);
    }

    cr3bp(double C=4.5, double mu=0.5) : m_bbox(vec6(m_inf), vec6(p_inf)) {
        m_bbox.min()[0] = -0.6985;
        m_bbox.min()[3] = -1.3335;
        m_bbox.max()[0] = -0.5080;
        m_bbox.max()[3] = 1.3335;
        m_C = C;
        m_mu = mu;
    }

    cr3bp(const bbox6& bbox)
            : m_bbox(bbox) {}

    vec6 operator()(const double& t, const vec6& y) const {
        vec6 r;
        if (!interpolate(y, r))
            throw undefined_point();

        return r;
    }

    nvis::vec2 project(const vec6& v) const {
        return nvis::vec2(v[0], v[3]);
    }

    vec6 unproject(const nvis::vec2& v) const {
        vec6 v6(0);
        v6[0] = v[0];
        v6[3] = v[1];
        v6[4] = yd(v[0], v[1]);
        return v6;
    }

    vec7 plane() const {
        vec7 p(0);
        p[1] = 1;
        return p;
    }

    cr3bp* clone() const {
        return new cr3bp(m_bbox);
    }

    bbox6 bounds() const {
        return m_bbox;
    }

    // ---


private:

    bool interpolate(const vec6& p, vec6& v) const {
        const double& x =  p[0];
        const double& y =  p[1];
        const double& z =  p[2];
        const double& xd = p[3];
        const double& yd = p[4];
        const double& zd = p[5];

        double d  = sqrt((x+m_mu)*(x+m_mu) + y*y + z*z);
        double r  = sqrt((x-1+m_mu)*(x-1+m_mu) + y*y + z*z);
        double d3 = d*d*d; //pow(d, 3);
        double r3 = r*r*r; //pow(r, 3);

        v[0] =  xd;
        v[1] =  yd;
        v[2] =  zd;
        v[3] =  2 * yd  + x - (1 - m_mu) * (x + m_mu) / d3 - m_mu * (x - 1 + m_mu) / r3;
        v[4] = -2 * xd + y - (1 - m_mu) * y / d3 - m_mu * y / r3;
        v[5] = -(1 - m_mu) * z / d3 - m_mu * z / r3;

        return true;
    }

    nvis::bezier_spline<double, 6> rspline;
    bbox6 m_bbox;
    double m_C, m_mu;
};

#endif // __cr3bp_hpp
