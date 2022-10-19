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
    typedef fixed_vector<T, N>    position_type;
    typedef V                           data_type;
    
    data_point() {}
    data_point(const position_type& p) : m_point(p) {}
    data_point(const position_type& p, const data_type& d) : m_point(p), m_data(d) {}
    data_point(const data_point& dp) : m_point(dp.m_point), m_data(dp.m_data) {}
    
    const position_type& position() const {
        return m_point;
    }
    
    position_type& position() {
        return m_point;
    }
    
    const data_type& data() const {
        return m_data;
    }
    
    data_type& data() {
        return m_data;
    }
    
    value_type distance_to(const data_point& dp) const {
        return norm(dp.m_point - m_point);
    }
    
    value_type operator[](size_t n) const {
        return m_point[n];
    }
    
private:
    position_type   m_point;
    data_type       m_data;
};

// kdtree-based point locator - can be used for query and insertion
template<typename T, typename V, int N>
class point_locator {
public:
    typedef V                                           data_type;
    typedef fixed_vector<T, N>                    vec_type;
    typedef data_point<T, V, N>                         data_point_type;
    typedef KDTree::KDTree<N, data_point_type>          kdtree_type;
    typedef typename kdtree_type::const_iterator        const_iterator;
    
    point_locator() {}
    
    point_locator(const point_locator& pl) : m_kdtree(pl.m_kdtree) {}
    
    void insert(const data_point_type& dp) const {
        m_kdtree.insert(dp);
    }
    
    data_type find_close_point(const vec_type& x) {
        if (m_kdtree.empty()) {
            throw std::runtime_error("invalid query on empty tree");
        }
        
        std::pair<const_iterator, T> nearest = m_kdtree.find_nearest(data_point_type(x));
        return (nearest.first)->data();
    }
    
    mutable kdtree_type m_kdtree;
};



}


#endif





