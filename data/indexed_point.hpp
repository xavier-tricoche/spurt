#ifndef __INDEXED_POINT_HPP__
#define __INDEXED_POINT_HPP__

#include <iostream>
#include <math/fixed_vector.hpp>

/*
    A spatial coordinate with an index to be used in a tree data structure
*/


namespace spurt {

template<typename T, int N>
class indexed_point {
public:
    typedef T                           value_type;
    typedef fixed_vector<T, N>    vector_type;
    typedef size_t                      index_type;
    
    static size_t size() {
        return N;
    }
    
    indexed_point() : __v(), __idx(-1) {}
    indexed_point(const T& x) : __v(x), __idx(-1) {}
    indexed_point(const T& x, const T& y) : __v(x, y), __idx(-1) {}
    indexed_point(const T& x, const T& y, const T& z) : __v(x, y, z), __idx(-1) {}
    indexed_point(const T& x, const T& y, const T& z, const T& w) : __v(x, y, z, w), __idx(-1) {}
    indexed_point(const T& x, const T& y, const T& z, const T& w, const T& t) : __v(x, y, z, w, t), __idx(-1) {}
    indexed_point(const vector_type& v, index_type idx = -1) : __v(v), __idx(idx) {}
    indexed_point(const indexed_point& p) : __v(p.__v), __idx(p.__idx) {}
    
    const vector_type& pos() const {
        return __v;
    }
    
    vector_type& pos() {
        return __v;
    }
    
    index_type index() const {
        return __idx;
    }
    
    index_type& index() {
        return __idx;
    }
    
    value_type distance_to(const indexed_point& p) const {
        return norm(p.__v - __v);
    }
    
    value_type operator[](size_t n) const {
        return __v[n];
    }
    
private:
    vector_type     __v;
    index_type      __idx;
};

template<typename V, int N>
inline bool operator==(const indexed_point<V, N>& p0, const indexed_point<V, N>& p1)
{
    return (p0.pos() == p1.pos());
}

template<typename T, int N>
std::ostream& operator<<(std::ostream& os, const indexed_point<T, N>& pt)
{
    os << "(" << pt.pos() << ", " << pt.index() << ")";
    return os;
}

}

#endif
