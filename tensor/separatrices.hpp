#ifndef __SEPARATRICES_HPP__
#define __SEPARATRICES_HPP__

#include <iostream>
#include <iomanip>

#include <vector>
#include <complex>
#include <sstream>
#include <math.h>

#include <boost/format.hpp>
#include <boost/limits.hpp>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>

#include <math/math.hpp>

#include <iostream>
#include <list>

#include <data/grid.hpp>
#include <data/raster_data.hpp>

#include <util/wall_timer.hpp>

#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <tensor/eigenvector_field.hpp>
#include <tensor/tensor_math.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spurt {
template<typename T>
nvis::vec3 eigenplane(const T& field, const nvis::vec3& x, bool linear_type)
{
    EigenvectorField<T> efield(field);
    
    if (linear_type) {
        return efield.evec(x, 0);
    } else {
        return efield.evec(x, 2);
    }
}

nvis::vec2 eigen2d(const nvis::mat2& m, bool major)
{
    double alpha = 0.5*(m(0,0)-m(1,1));
    double beta = m(1,0);
    
    if (major) {
        return nvis::vec2(beta, -alpha + sqrt(alpha*alpha + beta*beta));
    } else {
        return nvis::vec2(beta, -alpha - sqrt(alpha*alpha + beta*beta));
    }
}

void make_basis(nvis::vec3& b0, nvis::vec3& b1, const nvis::vec3& n)
{
    double maxnorm = 0;
    for (int i = 0 ; i < 3 ; ++i) {
        nvis::vec3 e(0);
        e[i] = 1;
        nvis::vec3 tmp = nvis::cross(e, n);
        if (nvis::norm(tmp) > maxnorm) {
            b0 = tmp;
            maxnorm = nvis::norm(tmp);
        }
    }
    
    b0 /= nvis::norm(b0);
    b1 = nvis::cross(n, b0);
}

nvis::mat2 project_matrix_on_plane(const nvis::mat3& m,
                                   const nvis::vec3& b0, const nvis::vec3& b1)
{
    nvis::mat2 m2;
    m2(0, 0) = nvis::inner(b0, m * b0);
    m2(1, 0) = nvis::inner(b1, m * b0);
    m2(0, 1) = nvis::inner(b0, m * b1);
    m2(1, 1) = nvis::inner(b1, m * b1);
    return m2;
}

template<typename T>
nvis::mat2 planar_interpolate(const T& field, const nvis::vec3& x,
                              const nvis::vec3& b0, const nvis::vec3& b1)
{
    nvis::fixed_vector<double, 7> t = field(x);
    nvis::mat3 m = to_matrix(t);
    return project_matrix_on_plane(m, b0, b1);
}

template<typename T>
nvis::mat2 compute_2d_derivative(const T& field, const nvis::vec2& x,
                                 const nvis::vec3& e0, const nvis::vec3& e1)
{
    typedef T                                       field_type;
    typedef typename field_type::data_type          tensor_type;
    typedef typename field_type::derivative_type    derivative_type;
    
    nvis::vec3 normal = nvis::cross(e0, e1);
    derivative_type dT = field.derivative(x);
}

template<typename T>
void eigendirections(std::vector<nvis::vec3>& dir, const T& field, const nvis::vec3& x, double h, bool linear_type)
{
    typedef EigenvectorField<T>     efield_type;
    
    efield_type efield(field);
    nvis::vec3 n = eigenplane(field, x, linear_type);
    nvis::vec3 b0, b1;
    make_basis(b0, b1, n);
    
    nvis::mat2 dTdx = 0.5 / h * (planar_interpolate(field, x + h * b0, b0, b1) -
                                 planar_interpolate(field, x - h * b0, b0, b1));
    nvis::mat2 dTdy = 0.5 / h * (planar_interpolate(field, x + h * b1, b0, b1) -
                                 planar_interpolate(field, x - h * b1, b0, b1));
                                 
    double alpha0 = 0.5*(dTdx(0,0) - dTdx(1,1));
    double alpha1 = 0.5*(dTdy(0,0) - dTdy(1,1));
    double beta0 = dTdx(0,1);
    double beta1 = dTdy(0,1);
    std::complex<double> r[3];
    cubic_equation(beta1, beta0+2*alpha1, 2*alpha0-beta1, -beta0, r);
    
    std::vector<double> u;
    for (int i=0 ; i<3 ; ++i) {
        if (r[i].imag() == 0) {
            u.push_back(r[i].real());
        }
    }
    
    dir.clear();
    for (int i=0 ; i<u.size() ; ++i) {
        double theta = atan(u[i]);
        nvis::vec3 step = cos(theta)*b0 + sin(theta)*b1;
        nvis::vec3 y = x + h*step;
        nvis::mat2 t = planar_interpolate(field, y, b0, b1);
        nvis::vec2 c = eigen2d(t, linear_type);
        nvis::vec3 e = c[0]*b0 + c[1]*b1;
        e /= nvis::norm(e);
        if (nvis::norm(nvis::cross(e, step)) < 0.1) {
            dir.push_back(e);
        } else {
            dir.push_back(-1.*e);
        }
    }
}

} // spurt



#endif












