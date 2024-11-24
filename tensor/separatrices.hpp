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

#include <math/types.hpp>
#include <misc/progress.hpp>

#include <math/basic_math.hpp>

#include <iostream>
#include <list>

#include <data/mesh.hpp>
#include <data/image.hpp>

#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <tensor/eigenvector_field.hpp>
#include <tensor/tensor_math.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spurt {
template<typename T>
spurt::vec3 eigenplane(const T& field, const spurt::vec3& x, bool linear_type)
{
    EigenvectorField<T> efield(field);
    
    if (linear_type) {
        return efield.evec(x, 0);
    } else {
        return efield.evec(x, 2);
    }
}

spurt::vec2 eigen2d(const spurt::mat2& m, bool major)
{
    double alpha = 0.5*(m(0,0)-m(1,1));
    double beta = m(1,0);
    
    if (major) {
        return spurt::vec2(beta, -alpha + sqrt(alpha*alpha + beta*beta));
    } else {
        return spurt::vec2(beta, -alpha - sqrt(alpha*alpha + beta*beta));
    }
}

void make_basis(spurt::vec3& b0, spurt::vec3& b1, const spurt::vec3& n)
{
    double maxnorm = 0;
    for (int i = 0 ; i < 3 ; ++i) {
        spurt::vec3 e = 0;
        e[i] = 1;
        spurt::vec3 tmp = spurt::cross(e, n);
        if (spurt::norm(tmp) > maxnorm) {
            b0 = tmp;
            maxnorm = spurt::norm(tmp);
        }
    }
    
    b0 /= spurt::norm(b0);
    b1 = spurt::cross(n, b0);
}

spurt::mat2 project_matrix_on_plane(const spurt::mat3& m,
                                    const spurt::vec3& b0, const spurt::vec3& b1)
{
    spurt::mat2 m2;
    m2(0, 0) = spurt::inner(b0, m * b0);
    m2(1, 0) = spurt::inner(b1, m * b0);
    m2(0, 1) = spurt::inner(b0, m * b1);
    m2(1, 1) = spurt::inner(b1, m * b1);
    return m2;
}

template<typename T>
spurt::mat2 planar_interpolate(const T& field, const spurt::vec3& x,
                              const spurt::vec3& b0, const spurt::vec3& b1)
{
    spurt::fixed_vector<double, 7> t = field(x);
    spurt::mat3 m = to_matrix(t);
    return project_matrix_on_plane(m, b0, b1);
}

template<typename T>
spurt::mat2 compute_2d_derivative(const T& field, const spurt::vec2& x,
                                  const spurt::vec3& e0, const spurt::vec3& e1)
{
    typedef T                                       field_type;
    typedef typename field_type::data_type          tensor_type;
    typedef typename field_type::derivative_type    derivative_type;
    
    spurt::vec3 normal = spurt::cross(e0, e1);
    derivative_type dT = field.derivative(x);
}

template<typename T>
void eigendirections(std::vector<spurt::vec3>& dir, const T& field, 
                     const spurt::vec3& x, double h, bool linear_type)
{
    typedef EigenvectorField<T>     efield_type;
    
    efield_type efield(field);
    spurt::vec3 n = eigenplane(field, x, linear_type);
    spurt::vec3 b0, b1;
    make_basis(b0, b1, n);
    
    spurt::mat2 dTdx = 0.5 / h * (planar_interpolate(field, x + h * b0, b0, b1) -
                                 planar_interpolate(field, x - h * b0, b0, b1));
    spurt::mat2 dTdy = 0.5 / h * (planar_interpolate(field, x + h * b1, b0, b1) -
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
        spurt::vec3 step = cos(theta)*b0 + sin(theta)*b1;
        spurt::vec3 y = x + h*step;
        spurt::mat2 t = planar_interpolate(field, y, b0, b1);
        spurt::vec2 c = eigen2d(t, linear_type);
        spurt::vec3 e = c[0]*b0 + c[1]*b1;
        e /= spurt::norm(e);
        if (spurt::norm(spurt::cross(e, step)) < 0.1) {
            dir.push_back(e);
        } else {
            dir.push_back(-1.*e);
        }
    }
}

} // spurt



#endif












