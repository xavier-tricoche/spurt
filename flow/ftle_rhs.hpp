#ifndef __XAVIER_FTLE_RHS_HPP__
#define __XAVIER_FTLE_RHS_HPP__

namespace spurt {
    
template<typename Value_, typename State_>
struct Catseye {
    typedef Value_ value_t;
    typedef State_ state_t;
    
    Catseye(double c) : m_c(c), m_cc(sqrt(c*c-1)) {}

    // phi(x,y)=-log(c*cosh(y)+sqrt(c*c-1)*cos(x)):=-log(K(x,y))
    // dy phi(x,y)=-c*sinh(y)/K(x,y)
    // dx phi(x,y)=sqrt(c*c-1)*sin(x)/K(x,y)
    // W o phi(x,y)=1/K(x,y)
    
    void operator()(const state_t& x, state_t& dxdt, value_t) const {
        value_t tmp=1./(m_c*cosh(x[1]) + m_cc*cos(x[0]));
        dxdt[0] = -m_c*sinh(x[1])*tmp;
        dxdt[1] = m_cc*sin(x[0])*tmp;
        dxdt[2] = tmp;
    }

    value_t m_c, m_cc;
};

template<typename Value_, typename State_>
struct DoubleGyre {
    typedef Value_ value_t;
    typedef State_ state_t;
    
    static constexpr double PI=3.1415926535897932384626;
    
    DoubleGyre(value_t A=0.1, value_t eps=0.25, value_t omega=0.2*PI)
        : m_A(A), m_eps(eps), m_omega(omega) {}
    
    value_t a(value_t t) const {
        return m_eps*sin(m_omega*t);
    }
    
    value_t b(value_t t) const {
        return 1 - 2*m_eps*sin(m_omega*t);
    }
    
    value_t f(value_t x, value_t t) const {
        return a(t)*x*x + b(t)*x;
    }
    
    value_t dfdx(value_t x, value_t t) const {
        return 2*a(t)*x + b(t);
    }
    
    void operator()(const state_t& x, state_t& dxdt, value_t t) const {
        dxdt[0] = -PI*m_A*sin(PI*f(x[0],t))*cos(PI*x[1]);
        dxdt[1] = PI*m_A*cos(PI*f(x[0],t))*sin(PI*x[1])*dfdx(x[0], t);
    }
    
    value_t m_eps, m_omega, m_A;
};
    
    
} // spurt

#endif
