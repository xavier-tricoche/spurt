#ifndef __xmt_poincare_map_hpp
#define __xmt_poincare_map_hpp

#include <vector>
#include <exception>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdexcept>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <math/dopri5.hpp>

#include <tokamak/map2d.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace xavier {
class standard_map {
    static constexpr double twopi = 2.*M_PI;

    inline void forward(double& x, double& y) const {
        y -= _k / twopi * sin(twopi * x);
        x += y;
    }

    inline void backward(double& x, double& y) const {
        x -= y;
        y += _k / twopi * sin(twopi * x);
    }

    inline void forward_orig(double& x, double& y) const {
        y -= _k * sin(x);
        x += y;
    }

    inline void backward_orig(double& x, double& y) const {
        x -= y;
        y += _k * sin(x);
    }

public:
    typedef nvis::vec2 state_type;

    inline nvis::mat2 fjacobian(double x, double y) const {
        nvis::mat2 J;
        J(0,0) = 1 - _k*cos(twopi * x);
        J(0,1) = 1;
        J(1,0) = -_k*cos(twopi*x);
        J(1,1) = 1;
        return J;
    }

    inline nvis::mat2 bjacobian(double x, double y) const {
        nvis::mat2 J;
        J(0,0) = 1;
        J(0,1) = -1;
        J(1,0) = _k*cos(twopi*(x-y));
        J(1,1) = 1 - _k*cos(twopi*(x-y));
        return J;
    }

private:
    inline state_type _map(const state_type& x, bool fwd = true) const {
        state_type y = x;
        if (fwd) {
            forward(y[0], y[1]);
        } else {
            backward(y[0], y[1]);
        }
        return y;
    }

    inline state_type _map_orig(const state_type& x, bool fwd = true) const {
        state_type y = x;
        if (fwd) {
            forward_orig(y[0], y[1]);
        } else {
            backward_orig(y[0], y[1]);
        }
        return y;
    }

    inline std::pair<state_type,nvis::mat2> _map_and_jacobian(const state_type& x, bool fwd = true) const {
        std::pair<state_type, nvis::mat2> r;
        r.first = x;
        if (fwd) {
            r.second = fjacobian(x[0], x[1]);
            forward(r.first[0], r.first[1]);
        } else {
            r.second = bjacobian(x[0], x[1]);
            backward(r.first[0], r.first[1]);
        }
        return r;
    }

public:
    struct map_undefined : public std::runtime_error {
        map_undefined() : std::runtime_error("undefined map at this location") {}
    };

    standard_map(double k) : _k(k) {}

    double K() const {
        return _k;
    }

    state_type map(const state_type& x, int n = 1) const {
        state_type y = x;
        for (unsigned int i = 0 ; i < abs(n) ; ++i) {
            y = _map(y, (n > 0));
        }
        return y;
    }

    state_type map_orig(const state_type& x, int n = 1) const {
        state_type y = x;
        for (unsigned int i = 0 ; i < abs(n) ; ++i) {
            y = _map_orig(y, (n > 0));
        }
        return y;
    }

    void map(const state_type& x, std::vector< state_type >& hits, int n = 1) const {
        hits.resize(abs(n));
        state_type y = x;
        for (unsigned int i = 0; i < abs(n) ; ++i) {
            y = _map(y, (n > 0));
            hits[i] = y;
        }
    }

    void map(const state_type& x, std::vector<std::pair<state_type, nvis::mat2> >& out, int niter, double eps=0) const {
        out.resize(abs(niter));
        state_type y = x;
        for (unsigned int i = 0; i < abs(niter) ; ++i) {
            out[i] = _map_and_jacobian(y, (niter > 0));
            if (i>0) {
                out[i].second = out[i].second * out[i-1].second;
            }
            y = out[i].first;
        }
    }

    std::pair<state_type, nvis::mat2> map_and_jacobian(const state_type& in, int niter, double eps=0) const {
        std::vector<std::pair<state_type, nvis::mat2> > out;
        map(in, out, niter);
        if (out.size() != std::abs(niter)) {
            throw map_undefined();
        }
        return out.back();
    }

    const nvis::bbox2 bounds() const {
        return nvis::bbox2(state_type(0,0), state_type(1,1));
    }

    standard_map* clone() const {
        return new standard_map(*this);
    }

    double precision() const {
        return 0;
    }
    void precision(double) {}

private:
    double _k;
    double _dummy;
};

class symplectic4D {
public:
    typedef Eigen::Matrix<double, 4, 1> state_type;
    typedef Eigen::Matrix<double, 4, 4> deriv_type;
    typedef nvis::bounding_box<state_type> bounds_type;

    static constexpr double onepi = 3.1415926535897932384626433;
    static constexpr double twopi = 6.2831853071795864769252868;
    double k1, k2, eps;

private:
    void forward(double& p1, double& p2, double& q1, double& q2) const {
        q1 += p1;
        q2 += p2;
        p1 += k1/twopi * sin(twopi * q1) + eps/twopi * sin(twopi * (q1 + q2));
        p2 += k2/twopi * sin(twopi * q2) + eps/twopi * sin(twopi * (q1 + q2));
    }

    void backward(double& p1, double& p2, double& q1, double& q2) const {
        p1 -= k1/twopi * sin(twopi * q1) + eps/twopi * sin(twopi *(q1 + q2));
        p2 -= k2/twopi * sin(twopi * q2) + eps/twopi * sin(twopi *(q1 + q2));
        q1 -= p1;
        q2 -= p2;
    }

    deriv_type forwardJ(double p1, double p2, double q1, double q2) const {
        q1 += p1;
        q2 += p2;
        deriv_type J = deriv_type::Zero();
        J(0,0) = 1; // dp1'/dp1
        J(0,2) = k1 * cos(twopi * q1) + eps * cos(twopi * (q1 + q2)); // dp1'/dq1
        J(0,3) = J(1,2) = eps*cos(twopi*(q1+q2)); // dp1'/dq2 = dp2'/dq1
        J(1,1) = 1; // dp2'/dp2
        J(1,3) = k2 * cos(twopi * q2) + eps * cos(twopi * (q1 + q2)); // dp2'/dq2
        J(2,0) = J(2,2) = 1; // dq1'/dp1 = dq1'/dq1
        J(3,1) = J(3,3) = 1; // dq2'/dp2 = dq2'/dq2
        return J;
    }

    deriv_type backwardJ(double p1, double p2, double q1, double q2) const {
        deriv_type J = deriv_type::Zero();
        J(0,0) = 1; // dp1/dp1'
        J(0,2) = -k1 * cos(twopi * q1) - eps * cos(twopi * (q1 + q2)); // dp1/dq1'
        J(0,3) = J(1,2) = -eps * cos(twopi*( q1 + q2)); // dp1/dq2' = dp2/dq1'
        J(1,1) = 1; // dp2/dq1'
        J(1,3) = -k2 * cos(twopi * q2) - eps * cos(twopi * (q1 + q2)); // dp2/dq2'
        J(2,0) = -1;
        J(2,2) = 1 + k1 * cos(twopi * q1) + eps * cos(twopi * (q1 + q2)); // dq1/dq1'
        J(2,3) = J(3,2) = eps * cos(twopi * (q1 + q2)); // dq1/dq2' = dq2/dq1'
        J(3,1) = -1;
        J(3,3) = 1 + k2 * cos(twopi * q2) + eps * cos(twopi * (q1 + q2));  // dq2/dq2'
        return J;
    }

public:

    struct map_undefined : public std::runtime_error {
        map_undefined() : std::runtime_error("undefined map at this location") {}
    };

    symplectic4D(double k1, double k2, double eps) : k1(k1), k2(k2), eps(eps) {}

    state_type map(const state_type& x, int n = 1) const {
        state_type y(x);

        if (n>0) {
            for (int i=0; i<n; ++i)
                forward(y[0], y[1], y[2], y[3]);
        }
        else if (n<0) {
            for (int i=n; i<0; ++i)
                backward(y[0], y[1], y[2], y[3]);
        }
        return y;
    }

    void map(const state_type& x, std::vector< state_type >& hits, int n = 1) const {
        hits.resize(abs(n));
        state_type y(x);
        if (n>0) {
            for (int i=0; i<n; ++i) {
                forward(y[0], y[1], y[2], y[3]);
                hits[i] = y;
            }
        }
        else if (n<0) {
            for (int i=n; i<0; ++i) {
                backward(y[0], y[1], y[2], y[3]);
                hits[i] = y;
            }
        }
    }

    void map(const state_type& x, std::vector<std::pair<state_type, deriv_type> >& out, int niter, double eps=0) const {
        out.resize(abs(niter));
        state_type y = x;
        if (niter > 0) {
            for (int i=0; i<niter; ++i) {
                deriv_type J = forwardJ(y[0], y[1], y[2], y[3]);
                if (i>0) out[i].second = out[i-1].second * J;
                else out[i].second = J;
                forward(y[0], y[1], y[2], y[3]);
                out[i].first = y;
            }
        }
        else if (niter < 0) {
            for (int i=niter; i<0; ++i) {
                deriv_type J = backwardJ(y[0], y[1], y[2], y[3]);
                if (i>niter) out[i].second = out[i-1].second * J;
                else out[i].second = J;
                backward(y[0], y[1], y[2], y[3]);
                out[i].first = y;
            }
        }
    }

    std::pair<state_type, deriv_type> map_and_jacobian(const state_type& in, int niter, double eps=0) const {
        std::vector<std::pair<state_type, deriv_type> > out;
        map(in, out, niter);
        if (out.size() != std::abs(niter)) {
            throw map_undefined();
        }
        return out.back();
    }

    const bounds_type bounds() const {
        return bounds_type(state_type(-0.5,-0.5,-0.5,-0.5), state_type(0.5,0.5,0.5,0.5));
    }

    symplectic4D* clone() const {
        return new symplectic4D(*this);
    }

    double precision() const {
        return 0;
    }
    void precision(double) {}
};

template<typename Field>
class xmt_poincare_map {
public:
    typedef Field           field_type;
    typedef nvis::vec2      state_type;

    struct map_undefined : public std::runtime_error {
        map_undefined() : std::runtime_error("undefined map at this location") {}
    };

    xmt_poincare_map(const field_type& field) :
        __field(field), __prec(1e-6) {
    }

    xmt_poincare_map(const xmt_poincare_map& other) :
        __field(other.__field), __prec(other.__prec) {
    }

    void precision(const double& prec) {
        __prec = prec;
    }

    double precision() const {
        return __prec;
    }

    void map(const state_type& in, std::vector<state_type>& out, int niter) const {
        __map_rk56(in, out, niter);
    }

    void map(const state_type& in, std::vector<std::pair<state_type, nvis::mat2> >& out, int niter, double eps=0) const {
        __map_and_jacobian_rk56(in, out, niter, eps);
    }

    state_type map(const state_type& in, int niter) const {
        std::vector<state_type> out;
        map(in, out, niter);
        if (out.size() != std::abs(niter)) {
            throw map_undefined();
        }
        return out.back();
    }

    std::pair<state_type, nvis::mat2> map_and_jacobian(const state_type& in, int niter, double eps=0) const {
        std::vector<std::pair<state_type, nvis::mat2> > out;
        map(in, out, niter);
        if (out.size() != std::abs(niter)) {
            throw map_undefined();
        }
        return out.back();
    }

    const nvis::bbox2& bounds() const {
        return __field.bounds();
    }

    nvis::bvec2 periodic() const {
        return __field.periodic();
    }

    xmt_poincare_map* clone() const {
        return new xmt_poincare_map(*this);
    }

    const field_type& field() const {
        return __field;
    }

private:

    void __map_rk56(const state_type& in, std::vector<state_type>& out, int niter) const;
    void __map_and_jacobian_rk56(const state_type& in,
                                 std::vector<std::pair<state_type, nvis::mat2> >& out,
                                 int niter, double eps) const;

    nvis::mat3 __integrate_jacobian(const nvis::mat3& w,
                                    const nvis::dopri5<nvis::vec3>::step& step,
                                    double tmax, double eps) const;

    const field_type&    __field;
    double               __prec;
};


template<typename Field>
void xmt_poincare_map<Field>::__map_rk56(const state_type& in,
        std::vector<state_type>& out, int niter) const
{
    out.clear();
    out.reserve(std::abs(niter));

    nvis::dopri5<nvis::vec3> intg;

    intg.t = 0;
    intg.y = __field.unproject(in);
    intg.t_max = niter < 0 ? -1e9 : 1e9;
    intg.h = 1e-1;
    intg.abstol = __prec;
    intg.reltol = __prec;

    nvis::dopri5<nvis::vec3>::result res;
    nvis::dopri5<nvis::vec3>::step   step;

    unsigned int nsteps = 0;
    while (out.size() < std::abs(niter)) {
        ++nsteps;

        try {
            res = intg.do_step(__field, step);

            if (res != nvis::dopri5<nvis::vec3>::OK) {
                break;
            }
            if (__field.intersect(step.y0(), step.y1())) {
                double t0 = step.t0();
                double t1 = step.t1();

                nvis::vec3 y0 = step.y0();
                nvis::vec3 y1 = step.y1();

                double tmid;
                nvis::vec3   ymid;

                for (unsigned int i = 0; i < 10; ++i) {
                    tmid = (t0 + t1) / 2;
                    ymid = step.y(tmid);

                    if (__field.intersect(y0, ymid)) {
                        t1 = tmid;
                        y1 = ymid;
                    } else if (__field.intersect(ymid, y1)) {
                        t0 = tmid;
                        y0 = ymid;
                    } else {
                        assert(false);
                    }
                }
                out.push_back(__field.project(ymid));
            }
        } catch (std::runtime_error& e) {
#ifdef VERBOSE
            std::cerr << "exception caught: " << e.what() << "\n";
#endif
            break;
        }
    }
}

template<typename Field>
nvis::mat3 xmt_poincare_map<Field>::__integrate_jacobian(const nvis::mat3& w,
        const nvis::dopri5<nvis::vec3>::step& step, double tmax, double eps) const
{
    nvis::mat3 _w = w;
    double t0 = step.t0();
    double t1 = tmax;
    int n = 10;
    double dt = (t1 - t0) / 10.;
    if (eps > 0) {
        int m = floor(nvis::norm(step.y0() - step.y(t1))/eps);
        if (m > n) {
            dt = (t1-t0)/(double)m;
            n = m;
        }
    }
    nvis::mat3 J;
    for (int i=0 ; i<n ; ++i) {
        double t = t0 + (double)i * dt;
        nvis::vec3 x = step.y(t);
        J = __field.derivative(0, x);
        _w += dt * (J * _w);
    }
    return _w;
}

template<typename Field>
void xmt_poincare_map<Field>::__map_and_jacobian_rk56(const state_type& in,
        std::vector<std::pair<state_type, nvis::mat2> >& out, int niter, double eps) const
{
    out.clear();
    out.reserve(std::abs(niter));

    typedef std::pair<state_type, nvis::mat2>   value_type;

    nvis::dopri5<nvis::vec3> intg;

    intg.t = 0;
    intg.y = __field.unproject(in);
    intg.t_max = niter < 0 ? -1e9 : 1e9;
    intg.h = 1e-1;
    intg.abstol = __prec;
    intg.reltol = __prec;

    nvis::dopri5<nvis::vec3>::result res;
    nvis::dopri5<nvis::vec3>::step   step;

    nvis::mat3 w, dvdx, aux;
    w = nvis::mat3::identity(); // flow map derivative wrt initial conditions at time 0
    unsigned int nsteps = 0;
    while (out.size() < std::abs(niter)) {
        ++nsteps;

        try {
            res = intg.do_step(__field, step);

            if (res != nvis::dopri5<nvis::vec3>::OK) {
                break;
            }
            if (__field.intersect(step.y0(), step.y1())) {
                double t0 = step.t0();
                double t1 = step.t1();

                nvis::vec3 y0 = step.y0();
                nvis::vec3 y1 = step.y1();

                double tmid;
                nvis::vec3   ymid;

                for (unsigned int i = 0; i < 10; ++i) {
                    tmid = (t0 + t1) / 2;
                    ymid = step.y(tmid);

                    if (__field.intersect(y0, ymid)) {
                        t1 = tmid;
                        y1 = ymid;
                    } else if (__field.intersect(ymid, y1)) {
                        t0 = tmid;
                        y0 = ymid;
                    } else {
                        assert(false);
                    }
                }

                value_type v;
                v.first = __field.project(ymid);
                aux = __integrate_jacobian(w, step, tmid, eps);
                v.second = __field.project(aux);
                out.push_back(v);
            }

            w = __integrate_jacobian(w, step, step.t1(), eps);
        } catch (std::runtime_error& e) {
#ifdef VERBOSE
            std::cerr << "exception caught: " << e.what() << "\n";
#endif
            break;
        }
    }
}

} // namespace xavier

#endif // __xmt_poincare_map_hpp
