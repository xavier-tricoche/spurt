#ifndef ____BASIC_MATH_HPP
#define ____BASIC_MATH_HPP

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <complex>

namespace spurt {

inline double cross(const nvis::vec2& x, const nvis::vec2& y)
{
    return x[0]*y[1] - x[1]*y[0];
}

inline double degrees(const nvis::vec2& x)
{
    nvis::vec2 y = x / nvis::norm(x);
    double theta = acos(y[0]);
    if (y[1] < 0) {
        theta = 2 * M_PI - theta;
    }
    return theta * 180. / M_PI;
}

inline nvis::vec2 evec(double c00, double c01, double c11, double lambda)
{
    nvis::vec2 e;
    double a = c00 - lambda;
    double b = c11 - lambda;
    if (fabs(a) > fabs(b)) {
        e[0] = -c01;
        e[1] = a;
    } else {
        e[0] = b;
        e[1] = -c01;
    }
    return 1. / nvis::norm(e)*e;
}

inline nvis::vec2 nullspace(const nvis::mat2& A)
{
    nvis::vec2 row1(A[0]);
    nvis::vec2 row2(A[1]);
    nvis::vec2 e;
    if (nvis::norm(row1) > nvis::norm(row2)) {
        e = nvis::vec2(-row1[1], row1[0]);
    } else {
        e = nvis::vec2(-row2[1], row2[0]);
    }
    e /= nvis::norm(e);
    return e;
}

inline bool eigen(nvis::vec2 evecs[2], std::complex<double> evals[2], const nvis::mat2& J)
{
    typedef std::complex<double>    complex_type;
    
    double tr = nvis::trace(J);
    double det = nvis::det(J);
    double delta = tr * tr - 4 * det;
    
    if (delta >= 0) { // real eigenvalues: saddle type
        double lmin = 0.5 * (tr - sqrt(delta));
        double lmax = 0.5 * (tr + sqrt(delta));
        evals[0] = complex_type(lmin);
        evals[1] = complex_type(lmax);
        
        nvis::mat2 A(J);
        A[0][0] -= lmin;
        A[1][1] -= lmin;
        evecs[0] = nullspace(A);
        
        A = J;
        A[0][0] -= lmax;
        A[1][1] -= lmax;
        evecs[1] = nullspace(A);
        
        return true;
    } else {
        evals[0] = 0.5 * (tr - complex_type(0, sqrt(-delta)));
        evals[1] = 0.5 * (tr + complex_type(0, sqrt(-delta)));
    }
    return false;
}

inline bool eigen(nvis::vec2 evecs[2], const nvis::mat2& J)
{
    double tr = nvis::trace(J);
    double det = nvis::det(J);
    double delta = tr * tr - 4 * det;
    
    if (delta > 0) { // real eigenvalues: saddle type
        double lmin = 0.5 * (tr - sqrt(delta));
        double lmax = 0.5 * (tr + sqrt(delta));
        
        nvis::mat2 A(J);
        A[0][0] -= lmin;
        A[1][1] -= lmin;
        evecs[0] = nullspace(A);
        
        A = J;
        A[0][0] -= lmax;
        A[1][1] -= lmax;
        evecs[1] = nullspace(A);
        
        return true;
    }
    return false;
}

//  Compute the PCA of the displacement field in the vicinity of the
//  considered position.
template< typename Map_Step >
void PCA(const Map_Step& map_step, const nvis::vec2& x, double dx, double dy,
         nvis::vec2& emin, nvis::vec2& emax, double& lmin, double& lmax)
{
    nvis::vec2 e0(1, 0), e1(0, 1);
    emin = e0;
    emax = e1;
    nvis::vec2 f[8];
    lmin = lmax = 0.;
    try {
        int n = 0;
        for (int i = -1 ; i <= 1 ; ++i) {
            for (int j = -1 ; j <= 1 ; ++j) {
                // we do not include central position which is assumed to be
                // least reliable
                if (!i && !j) {
                    continue;
                }
                nvis::vec2 y = x + i * dx * e0 + j * dy * e1;
                f[n++] = map_step(y);
            }
        }
    } catch (...) {
        // something's wrong. done
        return;
    }
    
    // compute the PCA of these vectors
    double c00 = 0, c01 = 0, c11 = 0;
    for (unsigned int i = 0 ; i < 8 ; ++i) {
        c00 += f[i][0] * f[i][0];
        c11 += f[i][1] * f[i][1];
        c01 += f[i][0] * f[i][1];
    }
    
    double b = -(c00 + c11);
    double c = c00 * c11 - c01 * c01;
    double delta = b * b - 4.*c;
    lmin = 0.5 * (-b - sqrt(delta));
    lmax = 0.5 * (-b + sqrt(delta));
    
    emin = evec(c00, c01, c11, lmin);
    emax = evec(c00, c01, c11, lmax);
}



}


#endif







