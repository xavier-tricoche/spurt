#include <poincare/ls.hpp>

typedef xavier::matrix<double> mat;

inline double SIGN(double a, double b)
{
    return fabs(a) * ((b >= 0) ? 1.0 : -1.0);
}

inline double SQR(double a)
{
    return a*a;
}

inline double pythag(double a, double b)
{
    double absa = abs(a), absb = abs(b);
    return (absa > absb ? absa*sqrt(1.0 + SQR(absb / absa)) :
            (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa / absb))));
}

bool xavier::NR_svddcmp(const mat& A, std::vector<double> w, mat& U, mat& V)
{
// Given a matrix A[1..m][1..n], this routine computes its singular value
// decomposition, A = U.Diag(w).V^T

    unsigned int m = A._m, n = A._n;
    unsigned int flag, i, its, j, jj, k, l(0), nm(0);
    double anorm, c, f, g, h, s, scale, x, y, z;
    
    U = A;
    w.resize(n, n);
    V.resize(n, n);
    
    double rv1[n];
    
    g = scale = anorm = 0.0; // Householder reduction to bidiagonal form.
    for (i = 1; i <= n; i++) {
        l = i + 1;
        rv1[i-1] = scale * g;
        g = s = scale = 0.0;
        if (i <= m) {
            for (k = i; k <= m; k++) {
                scale += fabs(U(k - 1, i - 1));
            }
            if (scale) {
                for (k = i; k <= m; k++) {
                    U(k - 1, i - 1) /= scale;
                    s += U(k - 1, i - 1) * U(k - 1, i - 1);
                }
                f = U(i - 1, i - 1);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                U(i - 1, i - 1) = f - g;
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = i; k <= m; k++) {
                        s += U(k - 1, i - 1) * U(k - 1, j - 1);
                    }
                    f = s / h;
                    for (k = i; k <= m; k++) {
                        U(k - 1, j - 1) += f * U(k - 1, i - 1);
                    }
                }
                for (k = i; k <= m; k++) {
                    U(k - 1, i - 1) *= scale;
                }
            }
        }
        w[i-1] = scale * g;
        g = s = scale = 0.0;
        if (i <= m && i != n) {
            for (k = l; k <= n; k++) {
                scale += fabs(U(i - 1, k - 1));
            }
            if (scale) {
                for (k = l; k <= n; k++) {
                    U(i - 1, k - 1) /= scale;
                    s += U(i - 1, k - 1) * U(i - 1, k - 1);
                }
                f = U(i - 1, l - 1);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                U(i - 1, l - 1) = f - g;
                for (k = l; k <= n; k++) {
                    rv1[k-1] = U(i - 1, k - 1) / h;
                }
                for (j = l; j <= m; j++) {
                    for (s = 0.0, k = l; k <= n; k++) {
                        s += U(j - 1, k - 1) * U(i - 1, k - 1);
                    }
                    for (k = l; k <= n; k++) {
                        U(j - 1, k - 1) += s * rv1[k-1];
                    }
                }
                for (k = l; k <= n; k++) {
                    U(i - 1, k - 1) *= scale;
                }
            }
        }
        anorm = std::max(anorm, (fabs(w[i-1]) + fabs(rv1[i-1])));
    }
    for (i = n; i >= 1; i--) { //  Accumulation of right-hand transformations.
        if (i < n) {
            if (g) {
                for (j = l; j <= n; j++) { // Double division to avoid possible under ow.
                    V(j - 1, i - 1) = (U(i - 1, j - 1) / U(i - 1, l - 1)) / g;
                }
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = l; k <= n; k++) {
                        s += U(i - 1, k - 1) * V(k - 1, j - 1);
                    }
                    for (k = l; k <= n; k++) {
                        V(k - 1, j - 1) += s * V(k - 1, i - 1);
                    }
                }
            }
            for (j = l; j <= n; j++) {
                V(i - 1, j - 1) = V(j - 1, i - 1) = 0.0;
            }
        }
        V(i - 1, i - 1) = 1.0;
        g = rv1[i-1];
        l = i;
    }
    for (i = std::min(m, n); i >= 1; i--) { // Accumulation of left-hand transformations.
        l = i + 1;
        g = w[i-1];
        for (j = l; j <= n; j++) {
            U(i - 1, j - 1) = 0.0;
        }
        if (g) {
            g = 1.0 / g;
            for (j = l; j <= n; j++) {
                for (s = 0.0, k = l; k <= m; k++) {
                    s += U(k - 1, i - 1) * U(k - 1, j - 1);
                }
                f = (s / U(i - 1, i - 1)) * g;
                for (k = i; k <= m; k++) {
                    U(k - 1, j - 1) += f * U(k - 1, i - 1);
                }
            }
            for (j = i; j <= m; j++) {
                U(j - 1, i - 1) *= g;
            }
        } else for (j = i; j <= m; j++) {
                U(j - 1, i - 1) = 0.0;
            }
        ++U(i - 1, i - 1);
    }
    for (k = n; k >= 1; k--) { // Diagonalization of the bidiagonal form: Loop over
        // singular values, and over allowed iterations.
        for (its = 1; its <= 30; its++) {
            flag = 1;
            for (l = k; l >= 1; l--) { // Test for splitting.
                // Note that rv1[0] is always zero.
                nm = l - 1;
                if ((double)(fabs(rv1[l-1]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }
                if ((double)(fabs(w[nm-1]) + anorm) == anorm) {
                    break;
                }
            }
            if (flag) {
                c = 0.0; // Cancellation of rv1[0], if l > 1.
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i-1];
                    rv1[i-1] = c * rv1[i-1];
                    if ((double)(fabs(f) + anorm) == anorm) {
                        break;
                    }
                    g = w[i-1];
                    h = pythag(f, g);
                    w[i-1] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 1; j <= m; j++) {
                        y = U(j - 1, nm - 1);
                        z = U(j - 1, i - 1);
                        U(j - 1, nm - 1) = y * c + z * s;
                        U(j - 1, i - 1) = z * c - y * s;
                    }
                }
            }
            z = w[k-1];
            if (l == k) { // Convergence.
                if (z < 0.0) { // Singular value is made nonnegative.
                    w[k-1] = -z;
                    for (j = 1; j <= n; j++) {
                        V(j - 1, k - 1) = -V(j - 1, k - 1);
                    }
                }
                break;
            }
            if (its == 30) {
                return false;
            }
            
            x = w[l-1]; // Shift from bottom 2-by-2 minor.
            nm = k - 1;
            y = w[nm-1];
            g = rv1[nm-1];
            h = rv1[k-1];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0; // Next QR transformation:
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i-1];
                y = w[i-1];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j-1] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; jj++) {
                    x = V(jj - 1, j - 1);
                    z = V(jj - 1, i - 1);
                    V(jj - 1, j - 1) = x * c + z * s;
                    V(jj - 1, i - 1) = z * c - x * s;
                }
                z = pythag(f, h);
                w[j-1] = z;  // Rotation can be arbitrary if z = 0.
                if (z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                
                for (jj = 1; jj <= m; jj++) {
                    y = U(jj - 1, j - 1);
                    z = U(jj - 1, i - 1);
                    U(jj - 1, j - 1) = y * c + z * s;
                    U(jj - 1, i - 1) = z * c - y * s;
                }
            }
            rv1[l-1] = 0.0;
            rv1[k-1] = f;
            w[k-1] = x;
        }
    }
    
    return true;
}





