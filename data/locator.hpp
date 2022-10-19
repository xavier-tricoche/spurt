#ifndef __XAVIER_LOCATOR_HPP__
#define __XAVIER_LOCATOR_HPP__

#include <kdtree++/kdtree.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <stdexcept>
#include <limits>
#include <list>

namespace xavier {
    
// defines a data point associating a spatial coordinate with a value
template<typename T, typename V, int K>
class data_point {
public:
    typedef T                         value_type;
    typedef nvis::fixed_vector<T, K>  coord_type;
    typedef V                         data_type;

    data_point() {}
    data_point(const coord_type& c) : __c(c) {}
    data_point(const coord_type& c, const data_type& d) : __c(c), __d(d) {}
    data_point(const data_point& p) : __c(p.__c), __d(p.__d) {}
    const coord_type& coordinate() const {
        return __c;
    }
    coord_type& coordinate() {
        return __c;
    }
    const data_type& data() const {
        return __d;
    }
    data_type& data() {
        return __d;
    }
    value_type distance_to(const data_point& dp) const {
        return nvis::norm(dp.__c - __c);
    }
    value_type operator[](size_t n) const {
        return __c[n];
    }
    
private:
    coord_type __c;
    data_type  __d;
};

// kdtree-based point locator - can be used for query and insertion
template<typename T, typename V, int K>
class point_locator {
public:
    typedef data_point<T, V, K>                  point_type;
    typedef typename point_type::value_type      value_type;
    typedef typename point_type::coord_type      coord_type;
    typedef typename point_type::data_type       data_type;
    typedef KDTree::KDTree<K, point_type>        tree_type;
    typedef typename tree_type::const_iterator   const_iterator;
    typedef typename tree_type::_Region_         region_type;
    
    point_locator() {}
    
    point_locator(const point_locator& pl) : t(pl.t) {}
    
    template<typename Iterator_>
    point_locator(const Iterator_& begin, const Iterator_& end) 
        : t(begin, end) {}
        
    void insert(const point_type& p) const {
        t.insert(p);
    }

    size_t size() const {
        return t.size();
    }

    point_type find_nearest_point(const coord_type& c) {
        if (t.empty()) throw std::runtime_error("invalid query on empty tree");

        std::pair<const_iterator, value_type> found = t.find_nearest(point_type(c), std::numeric_limits<value_type>::max());
        return *found.first;
    }

    bool find_nearest_within_range(point_type& p, const coord_type& c, const value_type& dist) {
        if (t.empty()) return false;
        std::pair<const_iterator, value_type> found = t.find_nearest(point_type(c), dist);
        if (found.first == t.end()) return false;
        else {
            p = *found.first;
            return true;
        }
    }
    
    void find_within_range(std::list<point_type>& n, const coord_type& c, const value_type& dist) {
        t.find_within_range(point_type(c), dist, std::back_inserter(n));
    }

    void find_within_range(std::list<point_type>& n, const region_type& r) {
        t.find_within_range(r, std::back_inserter(n));
    }
    
    const_iterator begin() const { return t.begin(); }
    const_iterator end() const { return t.end(); }

    mutable tree_type t;
};
}


#endif