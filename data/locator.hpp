#ifndef __XAVIER_KDTREE_HPP__
#define __XAVIER_KDTREE_HPP__

#include <kdtree++/kdtree.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <iostream>
#include <stdexcept>

namespace spurt {
// defines a data point associating a spatial coordinate with a value
template<typename T, typename V, int N>
class data_point {
public:
    typedef T                           value_type;
    typedef nvis::fixed_vector<T, N>    position_type;
    typedef V                           data_type;
    
    data_point() {}
    data_point(const position_type& p) : __p(p) {}
    data_point(const position_type& p, const data_type& d) : __p(p), __d(d) {}
    data_point(const data_point& dp) : __p(dp.__p), __d(dp.__d) {}
    
    const position_type& position() const {
        return __p;
    }
    
    position_type& position() {
        return __p;
    }
    
    const data_type& data() const {
        return __d;
    }
    
    data_type& data() {
        return __d;
    }
    
    value_type distance_to(const data_point& dp) const {
        return nvis::norm(dp.__p - __p);
    }
    
    value_type operator[](size_t n) const {
        return __p[n];
    }
    
private:
    position_type   __p;
    data_type       __d;
};

// kdtree-based point locator - can be used for query and insertion
template<typename T, typename V, int N>
class point_locator {
public:
    typedef V                                           data_type;
    typedef nvis::fixed_vector<T, N>                    vec_type;
    typedef data_point<T, V, N>                         data_point_type;
    typedef KDTree::KDTree<N, data_point_type>          kdtree_type;
    typedef typename kdtree_type::const_iterator        const_iterator;
    
    point_locator() {}
    
    point_locator(const point_locator& pl) : kdtree(pl.kdtree) {}
    
    void insert(const data_point_type& dp) const {
        kdtree.insert(dp);
    }
    
    data_type find_close_point(const vec_type& x) {
        if (kdtree.empty()) {
            throw std::runtime_error("invalid query on empty tree");
        }
        
        std::pair<const_iterator, T> nearest = kdtree.find_nearest(data_point_type(x));
        return (nearest.first)->data();
    }
    
    mutable kdtree_type kdtree;
};



}


#endif





