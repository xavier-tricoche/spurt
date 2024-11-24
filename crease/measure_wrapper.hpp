#ifndef __MEASURE_WRAPPER_HPP__
#define __MEASURE_WRAPPER_HPP__

#include <teem/nrrd.h>
#include <math/types.hpp>
#include <math/bounding_box.hpp>
#include <vector>
#include <map>
#include <image/probe.hpp>
#include "crease.hpp"
#include <image/nrrd_wrapper.hpp>


//#define DEBUG

namespace spurt {

struct measured_values {
    vec3 g, Hg, e[3];
    double val, eval[3], conf;
};

struct eps_order : public eps_lexicographical_order {
    eps_order() : eps_lexicographical_order(1.0e-8) {}
};

inline mat3 to_mat3(const vec6& m) 
{
    return mat3({m[0], m[1], m[2], m[1], m[3], m[4], m[2], m[4], m[5]});
}

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

#endif