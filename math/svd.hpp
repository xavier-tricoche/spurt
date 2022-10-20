#ifndef __NR_SVD_HPP__
#define __NR_SVD_HPP__

#include "math/fixed_matrix.hpp"
#include "math/fixed_vector.hpp"
#include <Eigen/Dense>
#include <float.h>

namespace spurt {
template<typename T, int M, int N>
nvis::fixed_matrix<T, N, M> pseudo_inv_diag(const nvis::fixed_matrix<T, M, N>& Sigma)
{
    const double eps = DBL_EPSILON;
    T max = 0;
    for (int i=0 ; i<std::min(M, N) ; ++i) {
        max = std::max(max, Sigma(i,i));
    }
    nvis::fixed_matrix<T, N, M> inv(0);
    for (int i=0 ; i<std::min(M, N); ++i) {
        inv(i,i) = (Sigma(i,i)/max > eps ? 1./Sigma(i,i) : 0.);
    }
    return inv;
}
}

namespace eigen_svd {

using namespace Eigen;

template<typename Type1, typename Type2, typename Type3>
void svdcmp(const MatrixBase<Type1>& A, MatrixBase<Type2>& U,
            MatrixBase<Type1>& w, MatrixBase<Type3>& V)
{
    JacobiSVD<Type1> _jsvd(A);
    _jsvd.compute(A, ComputeFullU | ComputeFullV);
    U = _jsvd.matrixU();
    V = _jsvd.matrixV();
    for (int i=0 ; i<std::min(A.rows(), A.cols()) ; ++i) {
        w(i,i) = _jsvd.singularValues()(i);
    }
}
}

namespace nr_svd {
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) (a == 0.0 ? 0.0 : a*a)

inline double pythag(double a, double b)
{
    double fabsa,fabsb;
    fabsa=fabs(a);
    fabsb=fabs(b);
    if (fabsa > fabsb) {
        return fabsa*sqrt(1.0+SQR(fabsb/fabsa));
    } else {
        return (fabsb == 0.0 ? 0.0 : fabsb*sqrt(1.0+SQR(fabsa/fabsb)));
    }
}

template<typename T, int M, int N>
void svdcmp(const nvis::fixed_matrix<T, M, N>& A, nvis::fixed_vector<T, N>& w,
            nvis::fixed_matrix<T, M, N>& U, nvis::fixed_matrix<T, N>& V, double eps=1.0e-12)
{
    // From Numerical Recipes:
    //
    // Given the matrix stored in A[0..m-1][0..n-1], this routine computes its
    // singular value decomposition, A = U.W.V^T and stores the results in the
    // matrices U and V, and the vector w.
    
    nvis::fixed_matrix<T, M, N> u(A);
    
    bool flag;
    int i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, s, scale, x, y, z;
    nvis::fixed_vector<T, N> rv1(0);
    g = scale = anorm = 0.0;
    int m = M;
    int n = N;
    
    for(i = 0; i < n; i++) {
        l = i + 2;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if(i < m) {
            // Householder reduction to bidiagonal form.
            for(k = i; k < m; k++) {
                scale += fabs(u[k][i]);
            }
            if(scale != 0.0) {
                for(k = i; k < m; k++) {
                    u[k][i] /= scale;
                    s += u[k][i] * u[k][i];
                }
                f = u[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                u[i][i] = f - g;
                for(j = l - 1; j < n; j++) {
                    for(s = 0.0, k = i; k < m; k++) {
                        s += u[k][i] * u[k][j];
                    }
                    f = s / h;
                    for(k = i; k < m; k++) {
                        u[k][j] += f * u[k][i];
                    }
                }
                for(k = i; k < m; k++) {
                    u[k][i] *= scale;
                }
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if(i + 1 <= m && i + 1 != n) {
            for(k = l - 1; k < n; k++) {
                scale += fabs(u[i][k]);
            }
            if(scale != 0.0) {
                for(k = l - 1; k < n; k++) {
                    u[i][k] /= scale;
                    s += u[i][k] * u[i][k];
                }
                f = u[i][l - 1];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                u[i][l - 1] = f - g;
                for(k = l - 1; k < n; k++) {
                    rv1[k] = u[i][k] / h;
                }
                for(j = l - 1; j < m; j++) {
                    for(s = 0.0, k = l - 1; k < n; k++) {
                        s += u[j][k] * u[i][k];
                    }
                    for(k = l - 1; k < n; k++) {
                        u[j][k] += s * rv1[k];
                    }
                }
                for(k = l - 1; k < n; k++) {
                    u[i][k] *= scale;
                }
            }
        }
        anorm = std::max(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    for(i = n - 1; i >= 0; i--) {   // Accumulation of right-hand transformations.
        if(i < n - 1) {
            if(g != 0.0) {
                for(j = l; j < n; j++) { // Double division to avoid possible underflow.
                    V[j][i]=(u[i][j]/u[i][l])/g;
                }
                for(j = l; j < n; j++) {
                    for(s = 0.0, k = l; k < n; k++) {
                        s += u[i][k] * V[k][j];
                    }
                    for(k = l; k < n; k++) {
                        V[k][j] += s * V[k][i];
                    }
                }
            }
            for(j = l; j < n; j++) {
                V[i][j] = V[j][i] = 0.0;
            }
        }
        V[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for(i = std::min(m, n) - 1; i >= 0; i--) {
        l = i + 1;
        g = w[i];
        for(j = l; j < n; j++) {
            u[i][j] = 0.0;
        }
        if(g != 0.0) {
            g = 1.0 / g;
            for(j = l; j < n; j++) {
                // Accumulation of left-hand transformations.
                for(s = 0.0, k = l; k < m; k++) {
                    s += u[k][i] * u[k][j];
                }
                f = (s / u[i][i]) * g;
                for(k = i; k < m; k++) {
                    u[k][j] += f * u[k][i];
                }
            }
            for(j = i; j < m; j++) {
                u[j][i] *= g;
            }
        } else for(j = i; j < m; j++) {
                u[j][i] = 0.0;
            }
        ++u[i][i];
    }
    for(k = n - 1; k >= 0; k--) {
        for(its = 0; its < 30; its++) {
            flag = true;
            for(l = k; l >= 0; l--) {
                nm = l - 1;
                // Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
                // Test for splitting.
                if(l == 0 || fabs(rv1[l]) <= eps * anorm) {
                    flag = false;
                    break;
                }
                if(fabs(w[nm]) <= eps * anorm) {
                    break;
                }
            }
            if(flag) {
                c = 0.0; // Cancellation of rv1[l], if l > 0.
                s = 1.0;
                for(i = l; i < k + 1; i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if(fabs(f) <= eps * anorm) {
                        break;
                    }
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for(j = 0; j < m; j++) {
                        y=u[j][nm];
                        z=u[j][i];
                        u[j][nm]=y*c+z*s;
                        u[j][i]=z*c-y*s;
                    }
                }
            }
            z = w[k];
            if(l == k) { // Convergence.
                if(z < 0.0) { // Singular value is made nonnegative.
                    w[k] = -z;
                    for(j = 0; j < n; j++) {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }
            if(its == 29) {
                throw("no convergence in 30 svdcmp iterations");
            }
            x = w[l];   // Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c = s = 1.0; // Next QR transformation:
            for(j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for(jj = 0; jj < n; jj++) {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x * c + z * s;
                    V[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z; // Rotation can be arbitrary if z = 0.
                if(z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for(jj = 0; jj < m; jj++) {
                    y = u[jj][j];
                    z = u[jj][i];
                    u[jj][j] = y * c + z * s;
                    u[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    
    U = u;
}

template<typename T, int N>
inline nvis::fixed_matrix<T, N> to_mat(const nvis::fixed_vector<T, N>& w)
{
    nvis::fixed_matrix<T, N> r = nvis::fixed_matrix<T, N>::identity();
    for (int i=0 ; i<N ; ++i) {
        r(i,i) = w[i];
    }
    return r;
}

template<typename T, int N>
inline nvis::fixed_matrix<T, N> to_inv_mat(const nvis::fixed_vector<T, N>& w, double eps=1.0e-12)
{
    nvis::fixed_matrix<T, N> r = nvis::fixed_matrix<T, N>::identity();
    double max = *std::max_element(w.begin(), w.end());
    for (int i=0 ; i<N ; ++i) {
        r(i,i) = (w[i]/max > eps ? 1./w[i] : 0);
    }
    return r;
}

template<typename T, int M, int N>
inline nvis::fixed_matrix<T, M, N> pseudoinv(const nvis::fixed_matrix<T, M, N>& A, double eps=1.0e-12)
{
    nvis::fixed_matrix<T, M, N> U;
    nvis::fixed_matrix<T, N> V;
    nvis::fixed_vector<T, N> w;
    svdcmp<T, M, N>(A, w, U, V, eps);
    return V*to_inv_mat<double, 3>(w, eps)*nvis::transpose(U);
}

}

#endif