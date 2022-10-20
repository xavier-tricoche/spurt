#ifndef __RECONSTRUCTION_FUNCTIONS_HPP__
#define __RECONSTRUCTION_FUNCTIONS_HPP__

#include <reconstruction/definitions.hpp>

namespace spurt { namespace reconstruction {

// template<typename _Kernel, size_t _Dim>
// struct reconstruction_kernel {
//     static const size_t dimension = _Dim;
//     static const size_t kernel_size = _Kernel::size;
//
//     typedef _Kernel                                     kernel_type;
//     typedef typename kernel_type::scalar_type           scalar_type;
//     typedef typename kernel_type::size_type             size_type;
//     typedef nvis::fixed_vector<size_t, dimension>       ivec_type;
//     typedef nvis::fixed_vector<scalar_type, dimension>  point_type;
//     typedef nvis::fixed_vector<scalar_type, dimension>  derivative_type;
//
//     template<typename _Raster>
//     static scalar_type operator()(const point_type& local_c,
//                                   const ivec_type& index_c,
//                                   const _Raster& raster) {
//         scalar_type r(0);
//         for (size_t d=0 ; d<dimension ; ++d) {
//
//         }
//     }
//
//     // 1D convolution
//     value_type pr(0);
//     for (long i=index_c[d]-kernel_size ; i<=index_c[d]+kernel_size ; ++i) {
//         r
//     }
// };

template<typename T>
inline T tent(T x) {
    return 1-x;
}

template<typename T>
inline T dtent(T x) {
    return -1;
}

template<typename T>
inline double tent(T x, T h) {
    return tent(x/h);
}

template<typename T>
inline double dtent(T x, T h) {
    return -1/h;
}

template<typename T>
inline T gauss(T x) {
    static constexpr double INV_SQRT_2PI = 0.39894228040143267794;
    return INV_SQRT_2PI*exp(-0.5*x*x);
}

template<typename T>
inline T gauss(T x, T sigma) {
    return gauss(x/sigma)/sigma;
}

template<typename T>
inline T gauss(T x, T sigma, T mu) {
    return gauss((x-mu)/sigma)/sigma;
}

template<typename T>
inline T dgauss(T x) {
    return -x*gauss(x);
}

template<typename T>
inline T dgauss(T x, T sigma) {
    return -x/(sigma*sigma)*gauss(x, sigma);
}

template<typename T>
inline T dgauss(T x, T sigma, T mu) {
    return dgauss(x-mu, sigma);
}

template<typename T>
inline T sinc(T x) {
    T y = M_PI*x;
    return sin(y)/y;
}

template<typename T>
struct tent_1d {
    typedef T value_type;

    value_type operator()(const value_type& x) {
        return tent(x);
    }
    value_type operator()(const value_type& x,
                          const value_type& h) {
        return tent(x/h);
    }
    value_type derivative(const value_type& x) {
        return dtent(x);
    }
    value_type derivative(const value_type& x,
                          const value_type& h) {
        return dtent(x, h);
    }
};

template<typename T>
struct gauss_1d {
    typedef T value_type;

    value_type operator()(const value_type& x) {
        return gauss(x);
    }
    value_type operator()(const value_type& x,
                          const value_type& sigma) {
        return gauss(x, sigma);
    }
    value_type derivative(const value_type& x) {
        return dgauss(x);
    }
    value_type derivative(const value_type& x,
                          const value_type& sigma) {
        return dgauss(x, sigma);
    }
};

template<typename K, size_t N>
struct kernel_base {
    static const size_t dimension = N;

    typedef K                                 kernel_type;
    typedef typename kernel_type::value_type  value_type;
    typedef size_t                            size_type;
    typedef nvis::fixed_vector<value_type, N> point_type;
    typedef nvis::fixed_vector<value_type, N> derivative_type;

    value_type operator()(const point_type&) = 0;
    value_type operator()(const point_type&, const point_type&) = 0;

    derivative_type derivative(const point_type&) = 0;
    derivative_type derivative(const point_type&, const point_type&) = 0;
};

template<typename K, size_t N>
struct separable_kernel : public kernel_base<K, N> {
    typedef kernel_base<K, N> base_type;
    typedef typename base_type::point_type point_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::derivative_type derivative_type;
    typedef typename base_type::kernel_type kernel_type;

    static const size_t dimension=base_type::dimension;

    value_type operator()(const point_type& x) {
        value_type r(1);
        for (size_t d=0 ; d<dimension ; ++d) r *= kernel_type(x[d]);
        return r;
    }

    value_type operator()(const point_type& x,
                          const point_type& scale) {
        value_type r(1);
        for (size_t d=0 ; d<dimension ; ++d) r *= kernel_type(x[d]/scale[d]);
        return r;
    }

    value_type derivative(const point_type& x) {
        derivative_type r;
        for (size_t i=0 ; i<dimension ; ++i) {
            r[i] = kernel_type::derivative(x);
            for (size_t j=0 ; j<dimension ; ++j) {
                if (j != i) r[i] *= kernel_type(x);
            }
        }
        return r;
    }

    derivative_type derivative(const point_type& x,
                               const point_type& scale) {
        derivative_type r;
        for (size_t i=0 ; i<dimension ; ++i) {
            r[i] = kernel_type::derivative(x, scale);
            for (size_t j=0 ; j<dimension ; ++j) {
                if (j != i) r[i] *= kernel_type(x, scale);
            }
        }
        return r;
    }
};

template<typename K, size_t N>
class radial_function_kernel : public kernel_base<K, N> {
    typedef kernel_base<K, N> base_type;
    typedef typename base_type::point_type point_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::derivative_type derivative_type;
    typedef typename base_type::kernel_type kernel_type;

    value_type operator()(const point_type& x) {
        return kernel_type(nvis::norm(x));
    }
    value_type operator()(const point_type& x,
                          const point_type& scale) {
        return kernel_type(nvis::norm(x/scale));
    }
    derivative_type derivative(const point_type& x) {
        value_type y = nvis::norm(x);
        return x*kernel_type::derivative(y)/y;
    }
};

template<typename T, size_t N> using tent_kernel =
    separable_kernel<tent_1d<T>, N>;

namespace franke_base {
    inline double f(double x, double a) {
        return (9.*x+a)*(9.*x+a);
    }

    inline double g(double x, double a, double k, double y, double b, double m) {
        return k*f(x,a) + m*f(y,b);
    }

    inline double f1(double x, double y) {
        return
            0.75*exp(g(x, -2, -0.25, y, -2, -0.25)) +
            0.75*exp(g(x, 1, -1./49., y, 1, -0.1)) +
            0.50*exp(g(x, -7, -0.25, y, -3, -0.25)) -
            0.20*exp(g(x, -4, -1, y, -7, -1));
    }

    inline double f7(double x, double y) {
        return
            2*cos(10*x)*sin(10*y) + sin(10*x*y);
    }

    inline double f10(double x, double y) {
        return
            exp(-.04*sqrt(sqr(80*x-40)+sqr(90*y-45))) *
            cos(.15*sqrt(sqr(80*x-40)+sqr(90*y-45)));
    }
}

    inline double franke(double x, double y, int fn=1)
    {
        using namespace franke_base;
        switch (fn) {
            case 1: return f1(x,y);
            case 7: return f7(x,y);
            case 10: return f10(x,y);
            default: {
                throw std::runtime_error("undefined Franke function");
            }
        }
    }

    namespace marschnerlobb_base {
        inline double rho_r(double r, double fM=6) {
            return cos(2*M_PI*fM*cos((M_PI*r)/2));
        }
    }

    inline double MarschnerLobb(double x, double y, double z, double alpha=0.25) {
        using namespace marschnerlobb_base;
        return
            (1-sin(M_PI*z/2)+alpha*(1+rho_r(sqrt(sqr(x)+sqr(y)))))/
            (2*(1+alpha));
    }

    inline double Sphere(double x, double y, double z, double radius=1) {
        return sqrt(x*x + y*y + z*z) - radius;
    }

    inline double Cylinder(double x, double y, double z, double radius=1) {
        return sqrt(x*x + y*y) - radius;
    }

    inline double Cone(double x, double y, double z, double radius=1, double height=1) {
        // at height h, radius is R(1-h/H)
        const double& h = z;
        double r = sqrt(x*x + y*y);
        if (h < 0) {
            if (r < radius) return -h;
            else return sqrt(h*h + (r-radius)*(r-radius));
        }
        else if (h < height) {
            double rr = radius*(1-h/height);
            return r-rr;
        }
        else return sqrt((h-height)*(h-height) + r*r);
    }

} // reconstruction
} // spurt

#endif
