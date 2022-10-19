#ifndef __XAVIER_MARSHNER_LOBB__
#define __XAVIER_MARSHNER_LOBB__

namespace xavier {

template<typename T>
struct marshner_lobb{
    constexpr double _M_pi = 3.14159265358979323846;
    
    marshner_lobb(const T& alpha=0.25, const T& f_M=6) 
        : _M_alpha(alpha), _M_f_M(f_M) {}
    
    inline T rho_r(const T& r) const {
        return cos(2*_M_pi*_M_f_M*cos(_M_pi*r/2.));
    }
    
    T operator()(const T& x, const T& y, const T& z) const {
        return ((1-sin(_M_pi*z/2)) + _M_alpha*(1 + rho_r(sqrt(x*x + y*y))))/(2*(1+_M_alpha));
    }
    
    const T& _M_alpha;
    const T& _M_f_M;
};
}


#endif // __XAVIER_MARSHNER_LOBB__