#ifndef __MEASURE_WRAPPER_HPP__
#define __MEASURE_WRAPPER_HPP__

#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <vector>
#include <map>
#include <image/probe.hpp>
#include "crease.hpp"

//#define DEBUG

namespace spurt {

struct measured_values {
    vec3 g, Hg, e[3];
    double val, eval[3], conf;
};

struct eps_order : public eps_lexicographical_order {
    eps_order() : eps_lexicographical_order(1.0e-8) {}
};

class MeasureWrapper {
public:
    typedef gage_interface::kernel_idx_t kernel_idx_t;
    
    MeasureWrapper(const Nrrd* nrrd, int aniso = 0, kernel_idx_t = gage_interface::C4_HEXIC);

    double confidence(const vec3& p) const;
    double value(const vec3& p) const;
    vec3 eigenvector(const vec3& p, unsigned int idx) const;
    vec3 gradient(const vec3& p) const;
    vec3 Hgradient(const vec3& p) const;
    double eigenvalue(const vec3& p, unsigned int idx) const;
    vec6 hessian(const vec3& p) const;

    int what_measure() const;

    // central difference based computation
    vec3 my_gradient(const vec3& p, double eps = 0.001) const;
    vec6 my_hessian(const vec3& p, double eps = 0.001) const;
    vec3 my_eigenvector(const vec3& p, unsigned int idx, double eps = 0.001) const;

    void turn_buffer_on() const;
    void turn_buffer_off() const;

private:
    gage_interface::scalar_wrapper *_scal_wrap;
    gage_interface::tensor_wrapper *_tens_wrap;
    mutable std::vector< vec3 > _evecs;
    mutable std::vector< double > _evals;
    mutable vec3 _grad;
    int _measure;

    double _confidence(const vec3& p) const;
    double _value(const vec3& p) const;
    vec3 _eigenvector(const vec3& p, unsigned int idx) const;
    vec3 _gradient(const vec3& p) const;
    double _eigenvalue(const vec3& p, unsigned int idx) const;

    bool is_stored(measured_values& values, const vec3& p) const;
    void store(measured_values& mv, const vec3& p) const;

    mutable bool active_buffer;
    mutable std::map< vec3, measured_values, eps_order > buffer;

    double my_directional_derivative(const vec3& p, const unsigned int dim, double eps = 0.001) const;
};

void test_derivatives(const Nrrd* nrrd, int aniso = 0);
};

inline double spurt::MeasureWrapper::
my_directional_derivative(const vec3& p, const unsigned int dim, double eps) const
{
    vec3 q(p);
    q[dim] += eps;
    double vplus = this->value(q);
    q[dim] -= 2 * eps;
    double vminus = this->value(q);
    return (vplus -vminus) / (2*eps);
}

inline spurt::vec3 spurt::MeasureWrapper::my_gradient(const vec3& p, double eps) const
{
    vec3 g;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        g[i] = my_directional_derivative(p, i, eps);
    }
    return g;
}

inline spurt::vec6 spurt::MeasureWrapper::my_hessian(const vec3& p, double eps) const
{
    vec6 h;
    vec3 g[6];
    vec3 q(p);
    for (unsigned int dim = 0 ; dim < 3 ; ++dim) {
        q[dim] += eps;
        g[2*dim] = my_gradient(q, eps);
        q[dim] -= 2 * eps;
        g[2*dim+1] = my_gradient(q, eps);
        q[dim] += eps;
    }
    unsigned int n = 0;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        for (unsigned int j = i ; j < 3 ; ++j) {
            h[n++] = (g[2*j][i] - g[2*j+1][i]) / (2 * eps);
        }
    }
    return h;
}

inline spurt::vec3 spurt::MeasureWrapper::
my_eigenvector(const vec3& p, unsigned int idx, double eps) const
{
    vec6 h = my_hessian(p, eps);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    return vec3(evec[3*idx], evec[3*idx+1], evec[3*idx+2]);
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

inline int spurt::MeasureWrapper::what_measure() const
{
    return _measure;
}

// --------------------------------------------------------------------------

inline void spurt::MeasureWrapper::turn_buffer_on() const
{
    active_buffer = true;
}

inline void spurt::MeasureWrapper::turn_buffer_off() const
{
    buffer.clear();
    active_buffer = false;
}

inline bool spurt::MeasureWrapper::is_stored(measured_values& values, const vec3& p) const
{
    std::map< vec3, measured_values, eps_order >::const_iterator it =
        buffer.find(p);
    if (it != buffer.end()) {
        values = it->second;
        return true;
    }

    return false;
}

inline void spurt::MeasureWrapper::store(measured_values& mv, const vec3& p) const
{
    mv.conf = this->_confidence(p);
    mv.val = this->_value(p);
    mv.g = this->_gradient(p);
    mv.Hg = spurt::crease::prod(this->hessian(p), mv.g);
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        mv.e[i] = this->_eigenvector(p, i);
        mv.eval[i] = this->_eigenvalue(p, i);
    }

    buffer[p] = mv;
}

// --------------------------------------------------------------------------

inline double spurt::MeasureWrapper::confidence(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);
        return mv.conf;
    }

    return this->_confidence(p);
}

// --------------------------------------------------------------------------

inline double spurt::MeasureWrapper::_confidence(const vec3& p) const
{
    double v;
    switch (_measure) {
    case 0: {
        v = 1;
        break;
    }
    case 1: {
        _tens_wrap->confidence(p, v);
        break;
    }
    case 2: {
        _tens_wrap->confidence(p, v);
        break;
    }
    default: {
        assert(false);
    }
    }

    return v;
}

// --------------------------------------------------------------------------

inline double spurt::MeasureWrapper::value(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);
        return mv.val;
    }

    return this->_value(p);
}

// --------------------------------------------------------------------------

inline double spurt::MeasureWrapper::_value(const vec3& p) const
{
    double v;

    switch (_measure) {
    case 0: {
        _scal_wrap->value(p, v);
        break;
    }
    case 1: {
        _tens_wrap->fa(p, v);
        break;
    }
    case 2: {
        _tens_wrap->mode(p, v);
        break;
    }
    default: {
        assert(false);
    }
    }

    return v;
}

// --------------------------------------------------------------------------

inline spurt::vec3 spurt::MeasureWrapper::eigenvector(const vec3& p, unsigned int idx) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.e[idx];
    }

    return this->_eigenvector(p, idx);
}

// --------------------------------------------------------------------------

inline spurt::vec3
spurt::MeasureWrapper::_eigenvector(const vec3& p, unsigned int idx) const
{
#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->hess_evecs(p, _evecs);
        break;
    }
    case 1: {
        _tens_wrap->fa_hess_evecs(p, _evecs);
        break;
    }
    case 2: {
        _tens_wrap->mode_hess_evecs(p, _evecs);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec6 h = my_hessian(p);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    vec3 ev(evec[3*idx], evec[3*idx+1], evec[3*idx+2]);
    // std::cout << "measured eigenvector = " << _evecs[idx] << ", approx eigenvector = " << ev << std::endl;

    return ev;
#endif

    return _evecs[idx];
}

// --------------------------------------------------------------------------

inline double spurt::MeasureWrapper::eigenvalue(const vec3& p, unsigned int idx) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.eval[idx];
    }

    return this->_eigenvalue(p, idx);
}

// --------------------------------------------------------------------------

inline double
spurt::MeasureWrapper::_eigenvalue(const vec3& p, unsigned int idx) const
{
    vec3 _tmp;

#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->hess_evals(p, _tmp);
        break;
    }
    case 1: {
        _tens_wrap->fa_hess_evals(p, _tmp);
        break;
    }
    case 2: {
        _tens_wrap->mode_hess_evals(p, _tmp);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec6 h = my_hessian(p);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    // std::cout << "measured eigenvalue = " << _tmp[idx] << ", approx eigenvalue = " << eval[idx] << std::endl;
    return eval[idx];
#endif

    return _tmp[idx];
}

// --------------------------------------------------------------------------

inline spurt::vec3 spurt::MeasureWrapper::gradient(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.g;
    }

    return this->_gradient(p);
}

// --------------------------------------------------------------------------

inline spurt::vec3 spurt::MeasureWrapper::Hgradient(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.Hg;
    }

    return spurt::crease::prod(this->hessian(p), this->_gradient(p));
}

// --------------------------------------------------------------------------

inline spurt::vec3 spurt::MeasureWrapper::_gradient(const vec3& p) const
{
#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->gradient(p, _grad);
        break;
    }
    case 1: {
        _tens_wrap->fa_grad(p, _grad);
        break;
    }
    case 2: {
        _tens_wrap->mode_grad(p, _grad);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec3 g = my_gradient(p);
    // std::cout << "measured gradient = " << _grad << ", approx gradient = " << g << std::endl;
    return g;
#endif

    return _grad;
}

// --------------------------------------------------------------------------

inline spurt::vec6 spurt::MeasureWrapper::hessian(const vec3& p) const
{
    vec6 _tmp;

#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->hessian(p, _tmp);
        break;
    }
    case 1: {
        _tens_wrap->fa_hess(p, _tmp);
        break;
    }
    case 2: {
        _tens_wrap->mode_hess(p, _tmp);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec6 h = my_hessian(p);
    // std::cout << "measured hessian = " << h << ", approx hessian = " << h << std::endl;
    return h;
#endif

    return _tmp;
}

#endif






































