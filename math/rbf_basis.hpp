#ifndef __XAVIER_RBF_BASIS_HPP__
#define __XAVIER_RBF_BASIS_HPP__

#include <math.h>
#include <boost/static_assert.hpp>
#include <iostream>
#include <cmath>

namespace spurt { namespace RBF {

template<typename T>
inline T linear(T x) { return x; }

template<typename T>
inline T quadratic(T x) { return x*x; }

template<typename T>
inline T cubic(T x) { return x*x*x; }

template<typename T>
inline T quartic(T x) { return quadratic<T>(x*x); }

template<typename T>
inline T quintic(T x) { return cubic<T>(x)*quadratic<T>(x); }

template<typename T>
inline T wendland(T r, T radius) {
    if (r >= radius) return 0;
    T x = r/radius;
    return quartic<T>(1 - x)*(static_cast<T>(4)*x + 1);
}

template<typename T>
inline T gaussian(T r, T epsilon) {
    return exp(-quadratic<T>(epsilon*r));
}

template<typename T>
inline T multiquadric(T r, T epsilon) {
    return sqrt(1 + quadratic<T>(epsilon*r));
}

template<typename T>
inline T inverse_quadratic(T r, T epsilon) {
    return 1/(1 + quadratic<T>(epsilon*r));
}

template<typename T, unsigned int N>
inline T polyharmonic(T r) {
    BOOST_STATIC_ASSERT( N >= 2 && !(N % 2) );
    return std::pow(r, N)*log(r);
}

template<typename T>
inline T dlinear(T x) { return 1; }

template<typename T>
inline T dquadratic(T x) { return static_cast<T>(2)*x; }

template<typename T>
inline T dcubic(T x) { return static_cast<T>(3)*x*x; }

template<typename T>
inline T dquartic(T x) { return static_cast<T>(4)*cubic<T>(x); }

template<typename T>
inline T dquintic(T x) { return static_cast<T>(5)*quartic<T>(x); }

template<typename T>
inline T dwendland(T r, T radius) {
    if (r >= radius) return 0;
    T x = r/radius;
    return -20.*x*cubic<T>(1 - x)/radius;
}

template<typename T>
inline T dgaussian(T r, T epsilon) {
    return -2*quadratic<T>(epsilon)*r*gaussian<T>(r, epsilon);
}

template<typename T>
inline T dmultiquadric(T r, T epsilon) {
    return quadratic<T>(epsilon)*r/multiquadric<T>(r, epsilon);
}

template<typename T>
inline T dinverse_quadratic(T r, T epsilon) {
    return -multiquadric<T>(r, epsilon)*inverse_quadratic<T>(r, epsilon);
}

template<typename T, unsigned int N>
inline T dpolyharmonic(T r) {
    BOOST_STATIC_ASSERT( N >= 2 && !(N % 2) );
    return std::pow(r, N-1)*(N*log(r) + 1);
}


template<typename T>
struct linear_function {
    T operator()(T r) const {
        return linear<T>(r);
    }
    
    static T derivative(T r) {
        return dlinear<T>(r);
    }
};

template<typename T>
struct quadratic_function {
    T operator()(T r) const {
        return quadratic<T>(r);
    }
    
    static T derivative(T r) {
        return dquadratic<T>(r);
    }
};

template<typename T>
struct cubic_function {
    T operator()(T r) const {
        return cubic<T>(r);
    }

    static T derivative(T r) {
        return dcubic<T>(r);
    }
};

template<typename T>
struct quartic_function {
    T operator()(T r) const {
        return quartic<T>(r);
    }

    static T derivative(T r) {
        return dquartic<T>(r);
    }
};

template<typename T>
struct quintic_function {
    T operator()(T r) const {
        return quintic<T>(r);
    }

    static T derivative(T r) {
        return dquintic<T>(r);
    }
};

// smooth (C^2) in dim 2 and 3
template<typename T>
struct wendland_function {
    wendland_function(T radius=1) : _radius(radius) {}
    
    T radius() const {
        return _radius;
    }
    
    T operator()(T r) const {
        return wendland<T>(r, _radius);
    }
    
    T derivative(T r) const {
        return dwendland<T>(r, _radius);
    }
    
    T _radius;
};

template<typename T>
struct gaussian_function {
    gaussian_function(T eps=1.0) : _eps(eps) {}  
    
    T operator()(T r) const {
        return gaussian<T>(r, _eps);
    }
    
    T derivative(T r) const {
        return dgaussian<T>(r, _eps);
    }
    
    T _eps;
};

template<typename T>
struct truncated_gaussian_function {    
    truncated_gaussian_function(T eps, T radius) 
        : _radius(radius), _eps(eps) {}
        
    T radius() const {
        return _radius;
    }
        
    T operator()(T r) const {
        if (r >= _radius) return static_cast<T>(0);
        return gaussian<T>(r, _eps);
    }
    
    T derivative(T r) const {
        if (r >= _radius) return static_cast<T>(0);
        return dgaussian<T>(r, _eps);
    }
    
    T _eps;
    T _radius;
};

template<typename T>
struct multiquadric_function {
    multiquadric_function(T epsilon) : _eps(epsilon) {}
    
    T operator()(T r) const {
        return multiquadric<T>(r, _eps);
    }
    
    T derivative(T r) const {
        return dmultiquadric<T>(r, _eps);
    }
    
    T _eps;
};

template<typename T>
struct inverse_quadratic_function {
    inverse_quadratic_function(T epsilon) : _eps(epsilon) {}
    
    T operator()(T r) const {
        return inverse_quadratic<T>(r, _eps);
    }
    
    T derivative(T r) const {
        return dinverse_quadratic<T>(r, _eps);
    }
    
    T _eps;
};

template<typename T, unsigned int N>
struct polyharmonic_function {
    BOOST_STATIC_ASSERT( N >= 2 && !(N % 2) );

    T operator()(T r) const {
        return polyharmonic<T, N>(r);
    }
    
    static T derivative(T r) {
        return dpolyharmonic<T, N>(r);
    }
};

} // namespace RBF

} // namespace spurt


#endif
