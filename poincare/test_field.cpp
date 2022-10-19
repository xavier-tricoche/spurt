#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <iostream>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <sstream>
#include <stdexcept>

nvis::mat3 convert(const nvis::vec9& J)
{
    nvis::mat3 r;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        for (unsigned int j = 0 ; j < 3 ; ++j) {
            r[i][j] = J[3*i+j];
        }
    }
    return r;
}

inline double val(const nvis::vec9& A, int i, int j)
{
    return A[3*i+j]; // row first
}

nvis::vec3 operator*(const nvis::vec9& A, const nvis::vec3& x)
{
    nvis::vec3 r(0, 0, 0);
    for (int i = 0 ; i < 3 ; ++i) {
        for (int j = 0 ; j < 3 ; ++j) {
            r[i] += val(A, i, j) * x[j];
        }
    }
}

template< typename T >
nvis::mat3 my_jacobian(const T& field, const nvis::vec3& x)
{
    nvis::mat3 J;
    double eps = 0.001;
    nvis::vec3 dx(eps, 0, 0);
    nvis::vec3 dy(0, eps, 0);
    nvis::vec3 dz(0, 0, eps);
    nvis::vec3 ddx = 0.5 / eps * (field(0, x + dx) - field(0, x - dx));
    nvis::vec3 ddy = 0.5 / eps * (field(0, x + dy) - field(0, x - dy));
    nvis::vec3 ddz = 0.5 / eps * (field(0, x + dz) - field(0, x - dz));
    
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        J[i][0] = ddx[i];
        J[i][1] = ddy[i];
        J[i][2] = ddz[i];
    }
    
    return J;
}

int main(int argc, char* argv[])
{
    tokamak_nimrod_parametric field("/scratch2/data/NIMROD/CDXU/hdf5/cdxub7n.h5",
                                    "2000");
                                    
    srand48(time(0));
    double error1 = 0., error2 = 0.;
    double maxerr1 = 0., maxerr2 = 0.;
    nvis::vec3 f, f1, f2, df, g;
    nvis::mat3 J1, J2;
    for (unsigned i = 0 ; i < 1000 ; ++i) {
        nvis::vec3 x(80.*drand48(), 120.*drand48(), 50.*drand48());
        try {
            f = field(0, x);
            J1 = field.derivative(0, x);
        } catch (...) {
            std::cerr << "exception caught at " << x << '\n';
            --i;
            continue;
        }
        nvis::vec3 dx(0.25*(2.*drand48() - 1.), 0.25*(2.*drand48() - 1), 0.25*(2.*drand48() - 1));
        f1 = f + J1 * dx;
        try {
            g = field(0, x + dx);
        } catch (...) {
            std::cerr << "exception caught at " << x + dx << '\n';
            --i;
            continue;
        }
        try {
            J2 = my_jacobian(field, x);
        } catch (...) {
            std::cerr << "exception caught at " << x << '\n';
            --i;
            continue;
        }
        f2 = f + J2 * dx;
        
        double err1 = nvis::norm(f1 - g) / nvis::norm(g);
        double err2 = nvis::norm(f2 - g) / nvis::norm(g);
        if (err1 > maxerr1) {
            maxerr1 = err1;
        }
        if (err2 > maxerr2) {
            maxerr2 = err2;
        }
        
        // std::cout << "f1 = " << f1 << ", f2 = " << f2 << ", freal = " << g << ",\nJ1 = " << J1
        //           << ",\n J2 = " << J2 << '\n';
        error1 += err1;
        error2 += err2;
        
        // std::cout << "|error1| = " << err1*100. << "%\n";
        // std::cout << "|error2| = " << err2*100. << "%\n";
    }
    
    error1 /= 1000.;
    error2 /= 1000.;
    
    std::cout << "average error 1 = " << error1 << '\n'
              << "max error 1 = " << maxerr1 << '\n';
    std::cout << "average error 2 = " << error2 << '\n'
              << "max error 2 = " << maxerr2 << '\n';
              
    return 0;
}




























