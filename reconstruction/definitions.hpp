#ifndef __RECONSTRUCTION_DEFINTIONS_HPP__
#define __RECONSTRUCTION_DEFINTIONS_HPP__

#include <math/fixed_vector.hpp>

namespace xavier { namespace reconstruction {

template<typename T, int N>
class point {
public:
    typedef T                             value_type;
    typedef nvis::fixed_vector<T, N>      vector_type;
    typedef size_t                        index_type;

    point() : __i(), __v() {}    
    point(const index_type& i, const vector_type& v) : __i(i), __v(v) {}
    point(const point& p) : __i(p.__i), __v(p.__v) {}

    const vector_type& pos() const {
        return __v;
    }

    vector_type& pos() {
        return __v;
    }

    index_type idx() const {
        return __i;
    }

    index_type& idx() {
        return __i;
    }

    value_type distance_to(const point& p) const {
        return norm(__v - p.__v);
    }

    value_type operator[](size_t i) const {
        return __v[i];
    }

private:
    index_type    __i;
    vector_type   __v;
};

inline double sqr(double x){
    return x*x;
}

inline double cube(double x){
    return x*x*x;
}

} // reconstruction
} // xavier


#endif