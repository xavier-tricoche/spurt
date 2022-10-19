#ifndef __MEASURE_WRAPPER_HPP__
#define __MEASURE_WRAPPER_HPP__

#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <vector>
#include <map>
#include <image/probe.hpp>
#include "crease.hpp"
#include <image/nrrd_wrapper.hpp>


//#define DEBUG

namespace xavier {

struct measured_values {
    nvis::vec3 g, Hg, e[3];
    double val, eval[3], conf;
};

struct eps_order : public nvis::eps_lexicographical_order {
    eps_order() : nvis::eps_lexicographical_order(1.0e-8) {}
};

class MeasureWrapper {
public:
    typedef gage_interface::kernel_idx_t kernel_idx_t;
    
    MeasureWrapper(const Nrrd* nrrd, int aniso = 0, kernel_idx_t = gage_interface::C4_HEXIC);

    double confidence(const nvis::vec3& p) const;
    double value(const nvis::vec3& p) const;
    nvis::vec3 eigenvector(const nvis::vec3& p, unsigned int idx) const;
    nvis::vec3 gradient(const nvis::vec3& p) const;
    nvis::vec3 Hgradient(const nvis::vec3& p) const;
    double eigenvalue(const nvis::vec3& p, unsigned int idx) const;
    nvis::vec6 hessian(const nvis::vec3& p) const;

    int what_measure() const;

    // central difference based computation
    nvis::vec3 my_gradient(const nvis::vec3& p, double eps = 0.001) const;
    nvis::vec6 my_hessian(const nvis::vec3& p, double eps = 0.001) const;
    nvis::vec3 my_eigenvector(const nvis::vec3& p, unsigned int idx, double eps = 0.001) const;

    void turn_buffer_on() const;
    void turn_buffer_off() const;

private:
    gage_interface::scalar_wrapper *_scal_wrap;
    gage_interface::tensor_wrapper *_tens_wrap;
    mutable std::vector< nvis::vec3 > _evecs;
    mutable std::vector< double > _evals;
    mutable nvis::vec3 _grad;
    int _measure;

    double _confidence(const nvis::vec3& p) const;
    double _value(const nvis::vec3& p) const;
    nvis::vec3 _eigenvector(const nvis::vec3& p, unsigned int idx) const;
    nvis::vec3 _gradient(const nvis::vec3& p) const;
    double _eigenvalue(const nvis::vec3& p, unsigned int idx) const;

    bool is_stored(measured_values& values, const nvis::vec3& p) const;
    void store(measured_values& mv, const nvis::vec3& p) const;

    mutable bool active_buffer;
    mutable std::map< nvis::vec3, measured_values, eps_order > buffer;

    double my_directional_derivative(const nvis::vec3& p, const unsigned int dim, double eps = 0.001) const;
};

void test_derivatives(const Nrrd* nrrd, int aniso = 0);
};

inline double xavier::MeasureWrapper::
my_directional_derivative(const nvis::vec3& p, const unsigned int dim, double eps) const
{
    nvis::vec3 q(p);
    q[dim] += eps;
    double vplus = this->value(q);
    q[dim] -= 2 * eps;
    double vminus = this->value(q);
    return (vplus -vminus) / (2*eps);
}

inline nvis::vec3 xavier::MeasureWrapper::my_gradient(const nvis::vec3& p, double eps) const
{
    nvis::vec3 g;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        g[i] = my_directional_derivative(p, i, eps);
    }
    return g;
}

inline nvis::vec6 xavier::MeasureWrapper::my_hessian(const nvis::vec3& p, double eps) const
{
    nvis::vec6 h;
    nvis::vec3 g[6];
    nvis::vec3 q(p);
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

inline nvis::vec3 xavier::MeasureWrapper::
my_eigenvector(const nvis::vec3& p, unsigned int idx, double eps) const
{
    nvis::vec6 h = my_hessian(p, eps);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    return nvis::vec3(evec[3*idx], evec[3*idx+1], evec[3*idx+2]);
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

inline int xavier::MeasureWrapper::what_measure() const
{
    return _measure;
}

// --------------------------------------------------------------------------

inline void xavier::MeasureWrapper::turn_buffer_on() const
{
    active_buffer = true;
}

inline void xavier::MeasureWrapper::turn_buffer_off() const
{
    buffer.clear();
    active_buffer = false;
}

inline bool xavier::MeasureWrapper::is_stored(measured_values& values, const nvis::vec3& p) const
{
    std::map< nvis::vec3, measured_values, eps_order >::const_iterator it =
        buffer.find(p);
    if (it != buffer.end()) {
        values = it->second;
        return true;
    }

    return false;
}

inline void xavier::MeasureWrapper::store(measured_values& mv, const nvis::vec3& p) const
{
    mv.conf = this->_confidence(p);
    mv.val = this->_value(p);
    mv.g = this->_gradient(p);
    mv.Hg = xavier::crease::prod(this->hessian(p), mv.g);
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        mv.e[i] = this->_eigenvector(p, i);
        mv.eval[i] = this->_eigenvalue(p, i);
    }

    buffer[p] = mv;
}

// --------------------------------------------------------------------------

inline double xavier::MeasureWrapper::confidence(const nvis::vec3& p) const
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

inline double xavier::MeasureWrapper::_confidence(const nvis::vec3& p) const
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

inline double xavier::MeasureWrapper::value(const nvis::vec3& p) const
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

inline double xavier::MeasureWrapper::_value(const nvis::vec3& p) const
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

inline nvis::vec3 xavier::MeasureWrapper::eigenvector(const nvis::vec3& p, unsigned int idx) const
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

inline nvis::vec3
xavier::MeasureWrapper::_eigenvector(const nvis::vec3& p, unsigned int idx) const
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
    nvis::vec6 h = my_hessian(p);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    nvis::vec3 ev(evec[3*idx], evec[3*idx+1], evec[3*idx+2]);
    // std::cout << "measured eigenvector = " << _evecs[idx] << ", approx eigenvector = " << ev << std::endl;

    return ev;
#endif

    return _evecs[idx];
}

// --------------------------------------------------------------------------

inline double xavier::MeasureWrapper::eigenvalue(const nvis::vec3& p, unsigned int idx) const
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
xavier::MeasureWrapper::_eigenvalue(const nvis::vec3& p, unsigned int idx) const
{
    nvis::vec3 _tmp;

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
    nvis::vec6 h = my_hessian(p);
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

inline nvis::vec3 xavier::MeasureWrapper::gradient(const nvis::vec3& p) const
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

inline nvis::vec3 xavier::MeasureWrapper::Hgradient(const nvis::vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.Hg;
    }

    return xavier::crease::prod(this->hessian(p), this->_gradient(p));
}

// --------------------------------------------------------------------------

inline nvis::vec3 xavier::MeasureWrapper::_gradient(const nvis::vec3& p) const
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
    nvis::vec3 g = my_gradient(p);
    // std::cout << "measured gradient = " << _grad << ", approx gradient = " << g << std::endl;
    return g;
#endif

    return _grad;
}

// --------------------------------------------------------------------------

inline nvis::vec6 xavier::MeasureWrapper::hessian(const nvis::vec3& p) const
{
    nvis::vec6 _tmp;

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
    nvis::vec6 h = my_hessian(p);
    // std::cout << "measured hessian = " << h << ", approx hessian = " << h << std::endl;
    return h;
#endif

    return _tmp;
}

#endif






































