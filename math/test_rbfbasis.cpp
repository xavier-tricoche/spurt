#include "RBFbasis.hpp"
#include <iostream>
#include <limits>
#include <util/timer.hpp>

template<typename T>
inline T quadratic(T x) {
    return x*x;
}

template<typename T>
inline T dquadratic(T x) {
    return 2*x;
}

template<typename T>
inline T cubic(T x) {
    return x*x*x;
}

template<typename T>
inline T dcubic(T x) {
    return 3*x*x;
}

template<typename T>
inline T quartic(T x) {
    return x*x*x*x;
}

template<typename T>
inline T dquartic(T x) {
    return 4*x*x*x;
}

template<typename T>
inline T quintic(T x) {
    return x*x*x*x*x;
}

template<typename T>
inline T dquintic(T x) {
    return 5*x*x*x*x;
}

namespace xrbf = xavier::RBF;

template<typename T>
struct type_traits {};

template<>
struct type_traits<float> {
    static std::string name() { return "float"; }
};

template<>
struct type_traits<double> {
    static std::string name() { return "double"; }
};

template<typename T>
void compare(size_t N) {
    double time_mono = 0;
    double time_ref  = 0;
    
    const T LARGE = pow(std::numeric_limits<T>::max(), 1./6.);
    
    srand48(time(0));
    
    std::cout << "Performance comparison for type=" << type_traits<T>::name() << ":\n";
    
    nvis::timer _timer;
    T x, y, r, err=0;
    double tic, toc;
    for (size_t i=0 ; i<N ; ++i) {
        x = LARGE*(-0.5 + drand48());
        tic = _timer.elapsed();
        y = xrbf::dquadratic<T>(x);
        toc = _timer.elapsed();
        time_mono += toc - tic;
        tic = _timer.elapsed();
        r = dquadratic<T>(x);
        toc = _timer.elapsed();
        time_ref += toc - tic;
        T _err = fabs(r-y)/fabs(r);
        err += _err;
        if (_err > 1.0e-6) {
            std::cout << x << "^2 = " << y << " vs " << r << '\n';
        }
    }
    err /= (T)N;
    
    std::cout << "Quadratic:\n";
    std::cout << "\ttemplated monomial computation took " << time_mono << " s. (" << (double)N/time_mono << " Hz)\n";
    std::cout << "\thardcoded expression took " << time_ref << " s. (" << (double)N/time_ref << " Hz)\n";
    std::cout << "\taverage discrepancy = " << err << '\n';
    
    time_mono = 0;
    time_ref  = 0;
    err = 0;
    for (size_t i=0 ; i<N ; ++i) {
        x = LARGE*(-0.5 + drand48());
        tic = _timer.elapsed();
        y = xrbf::dcubic<T>(x);
        toc = _timer.elapsed();
        time_mono += toc - tic;
        tic = _timer.elapsed();
        r = dcubic<T>(x);
        toc = _timer.elapsed();
        time_ref += toc - tic;
        err += fabs(r-y)/fabs(r);
    }
    err /= (T)N;
    
    std::cout << "Cubic:\n";
    std::cout << "\ttemplated monomial computation took " << time_mono << " s. (" << (double)N/time_mono << " Hz)\n";
    std::cout << "\thardcoded expression took " << time_ref << " s. (" << (double)N/time_ref << " Hz)\n";
    std::cout << "\taverage discrepancy = " << err << '\n';
    
    time_mono = 0;
    time_ref  = 0;
    err = 0;
    for (size_t i=0 ; i<N ; ++i) {
        x = LARGE*(-0.5 + drand48());
        tic = _timer.elapsed();
        y = xrbf::dquartic<T>(x);
        toc = _timer.elapsed();
        time_mono += toc - tic;
        tic = _timer.elapsed();
        r = dquartic<T>(x);
        toc = _timer.elapsed();
        time_ref += toc - tic;
        err += fabs(r-y)/fabs(r);
    }
    err /= (T)N;
    
    std::cout << "Quartic:\n";
    std::cout << "\ttemplated monomial computation took " << time_mono << " s. (" << (double)N/time_mono << " Hz)\n";
    std::cout << "\thardcoded expression took " << time_ref << " s. (" << (double)N/time_ref << " Hz)\n";
    std::cout << "\taverage discrepancy = " << err << '\n';
    
    time_mono = 0;
    time_ref  = 0;
    err = 0;
    for (size_t i=0 ; i<N ; ++i) {
        x = LARGE*(-0.5 + drand48());
        tic = _timer.elapsed();
        y = xrbf::dquintic<T>(x);
        toc = _timer.elapsed();
        time_mono += toc - tic;
        tic = _timer.elapsed();
        r = dquintic<T>(x);
        toc = _timer.elapsed();
        time_ref += toc - tic;
        err += fabs(r-y)/fabs(r);
    }
    err /= (T)N;
    
    std::cout << "Quintic:\n";
    std::cout << "\ttemplated monomial computation took " << time_mono << " s. (" << (double)N/time_mono << " Hz)\n";
    std::cout << "\thardcoded expression took " << time_ref << " s. (" << (double)N/time_ref << " Hz)\n";
    std::cout << "\taverage discrepancy = " << err << '\n';
}

int main(int argc, char* argv[]) {
    size_t ntries = 1000000;
    if (argc == 2) ntries = atoi(argv[1]);
    
    compare<float>(ntries);
    compare<double>(ntries);
    return 0;
}

