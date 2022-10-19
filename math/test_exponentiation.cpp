#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <util/timer.hpp>
#include <boost/lexical_cast.hpp>

template<typename T>
inline T exponentiation_by_squaring(T x, unsigned int n)
{
    T r(1);
    while (n) {
        if (n % 2) { // n is odd
            r *= x;
            n--;
        }
        x *= x;
        n = n >> 1;
    }
    return r;
}

template<typename T>
inline T hardcoded_exponentiation(T x, unsigned int n) {
    switch (n) {
        case 0:  return 1;
        case 1:  return x;
        case 2:  return x*x;
        case 3:  return x*x*x;
        case 4:  return x*x*x*x;
        case 5:  return x*x*x*x*x;
        case 6:  return x*x*x*x*x*x;
        case 7:  return x*x*x*x*x*x*x;
        case 8:  return x*x*x*x*x*x*x*x;
        case 9:  return x*x*x*x*x*x*x*x*x;
        case 10: return x*x*x*x*x*x*x*x*x*x;
        default: {
            unsigned int m = n >> 1;
            T half_pow = hardcoded_exponentiation<T>(x, m);
            half_pow += half_pow;
            if (n % 2) half_pow *= x;
            return half_pow;
        }
    }
}

template<typename T>
struct type_traits {};

template<>
struct type_traits<int> {
    static std::string name() { return "int"; }
};

template<>
struct type_traits<float> {
    static std::string name() { return "float"; }
};

template<>
struct type_traits<double> {
    static std::string name() { return "double"; }
};

template<typename T>
void compare(unsigned int n, unsigned int maxp) {
    T max_value = std::pow((double)std::numeric_limits<T>::max(), 
                           (double)1/(double)maxp);
    
    nvis::timer _timer;
    double err;
    
    // loop over exponents
    for (unsigned int p=2 ; p<=maxp ; ++p) {
        double t_squaring=0, t_hardcoded=0, t_pow=0;    
        double err_squaring=0, err_hardcoded=0;
        
        // repeat n times
        for (unsigned int i=0 ; i<n ; ++i) {
            T x = (-0.5 + drand48())*max_value;
            T r_squaring, r_hardcoded, r_pow;
            
            _timer.restart();
            r_squaring = exponentiation_by_squaring<T>(x, p);
            t_squaring += _timer.elapsed();
            
            _timer.restart();
            r_hardcoded = hardcoded_exponentiation<T>(x, p);
            t_hardcoded += _timer.elapsed();
            
            _timer.restart();
            r_pow = std::pow((double)x, (double)p);
            t_pow += _timer.elapsed();
            
            err = (double)(r_squaring - r_pow)/(double)r_pow;
            err_squaring += err*err;
            
            err = (double)(r_hardcoded - r_pow)/(double)r_pow;
            err_hardcoded += err*err;
        }
        
        std::cout << "type " << type_traits<T>::name() 
                  << " \texponent " << p << '\n'
                  << " \t\tsquaring:  " << t_squaring << " s. (" 
                  << (double)n/t_squaring << " Hz) \t"
                  << "error: " << 1/(double)n*sqrt(err_squaring) << '\n'
                  << " \t\thardcoded: " << t_hardcoded << " s. ("
                  << (double)n/t_hardcoded << " Hz) \t"
                  << "error: " << 1/(double)n*sqrt(err_hardcoded) << '\n'
                  << " \t\tstd::pow:  " << t_pow << " s. (" 
                  << (double)n/t_pow << " Hz) \t"
                  << "error: 0 (by construction)\n\n";
    }
}

std::string me;
void usage() {
    std::cout 
        << me << ": Compare different methods of exponentiation.\n\n"
        << "USAGE: " << me << " [options]\n\n"
        << "Options:\n"
        << " -h         Print this message\n"
        << " -n <int>   Number of iterations for benchmarking (default: 1e+6)\n"
        << " -p <int>   Maximum exponent (default: 20)\n"
        << "\n";
    exit(0);
}

bool has_enough(int at, int size, int n=1) {
    return at + n < size;
}

int main(int argc, char* argv[]) {
    srand48(time(0));
    
    me = argv[0];
    unsigned int n = 1000000;
    unsigned int p = 20;
    
    for (int i=1; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-n") {
            if (!has_enough(i, argc)) usage();
            n = atoi(argv[++i]);
        }
        else if (arg == "-p") {
            if (!has_enough(i, argc)) usage();
            p = atoi(argv[++i]);
        }
        else {
            std::cerr << "Unrecognized argument: " << arg << '\n';
            usage();
        }
    }

    compare<float>(n, p);
    compare<double>(n, p);
    
    return 0;
}
