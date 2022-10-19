#ifndef __QUADTREE_HPP__
#define __QUADTREE_HPP__

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace xavier {

template<typename T1, typename T2>
class quadtree {
public:
    typedef T1                                  index_type;
    typedef T2                                  point_type;
    typedef nvis::bounding_box<point_type>      bounds_type;
    typedef std::pair<point_type, index_type>   data_type;
    
    quadtree(const bounds_type& bounds,
             size_t max_depth, size_t max_nb_pts, size_t depth = 0)
        : __is_leaf(true), __depth(depth), __bounds(bounds),
          __center(bounds.center()), __points(),
          __Max_nb_pts(max_nb_pts), __Max_depth(max_depth) {
        __child[0] = __child[1] = __child[2] = __child[3] = 0;
    }
    
    quadtree(const quadtree& to_copy)
        : __is_leaf(to_copy.__is_leaf), __depth(to_copy.__depth), __bounds(to_copy.__bounds),
          __center(to_copy.__center), __points(to_copy.__points),
          __Max_nb_pts(to_copy.__Max_nb_pts), __Max_depth(to_copy.__Max_depth) {
        if (!__is_leaf) {
            for (int i = 0 ; i < 4 ; ++i) {
                __child[i] = new quadtree(*(to_copy.__child[i]));
            }
        }
    }
    
    ~quadtree() {
        if (!__is_leaf) {
            for (int i = 0 ; i < 4 ; ++i) {
                delete __child[i];
                __child[i] = 0;
            }
        }
    }
    
    const quadtree& insert(const point_type& x, const index_type& id) {
        assert(__bounds.inside(x));
        if (!__is_leaf) {
            assert(which_child(x) < 4);
            return __child[which_child(x)]->insert(x, id);
        } else if ((__is_leaf && __points.size() < __Max_nb_pts) || __depth >= __Max_depth) {
            __points.push_back(data_type(x, id));
            return *this;
        } else {
            split();
            assert(which_child(x) < 4);
            return __child[which_child(x)]->insert(x, id);
        }
    }
    
    int which_child(const point_type& x) const {
        return 2*(x[1] > __center[1]) + (x[0] > __center[0]);
    }
    
    void split() {
        assert(__is_leaf && __depth < __Max_depth);
        const point_type& min = __bounds.min();
        const point_type& max = __bounds.max();
        point_type p0 = point_type(__center[0], min[1]);
        point_type p1 = point_type(max[0], __center[1]);
        point_type p2 = point_type(__center[0], max[1]);
        point_type p3 = point_type(min[0], __center[1]);
        __child[0] = new quadtree(bounds_type(min, __center), __Max_depth, __Max_nb_pts, __depth + 1);
        __child[1] = new quadtree(bounds_type(p0, p1), __Max_depth, __Max_nb_pts, __depth + 1);
        __child[2] = new quadtree(bounds_type(p3, p2), __Max_depth, __Max_nb_pts, __depth + 1);
        __child[3] = new quadtree(bounds_type(__center, max), __Max_depth, __Max_nb_pts, __depth + 1);
        
        for (int i = 0 ; i < __points.size() ; ++i) {
            const point_type& x = __points[i].first;
            __child[which_child(x)]->__points.push_back(__points[i]);
        }
        
        __points.clear();
        __is_leaf = false;
    }
    
    const quadtree& find(const point_type& x) const {
        if (__is_leaf) {
            return *this;
        } else {
            return __child[which_child(x)]->find(x);
        }
    }
    
    quadtree& find(const point_type& x) {
        if (__is_leaf) {
            return *this;
        } else {
            return __child[which_child(x)]->find(x);
        }
    }
    
    const std::vector<data_type>& data() const {
        return __points;
    }
    
    size_t depth() const {
        return __depth;
    }
    
    size_t max_depth() const {
        return __Max_depth;
    }
    
    size_t max_bucket_size() const {
        return __Max_nb_pts;
    }
    
    const bounds_type bounds() const {
        return __bounds;
    }
    
    bool is_leaf() const {
        return __is_leaf;
    }
    
    template<typename Func>
    void for_each_leaf(Func& f) const {
        if (__is_leaf) {
            f(*this);
        } else {
            __child[0]->for_each_leaf(f);
            __child[1]->for_each_leaf(f);
            __child[3]->for_each_leaf(f);
            __child[2]->for_each_leaf(f);
        }
    }
    
private:
    bool                        __is_leaf;
    size_t                      __depth;
    bounds_type                 __bounds;
    point_type                  __center;
    std::vector<data_type>      __points;
    const size_t                __Max_nb_pts;
    const size_t                __Max_depth;
    quadtree*                   __child[4];
};

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& _pair)
{
    os << "( " << _pair.first << " , " << _pair.second << " )";
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const quadtree<T1, T2>& qt)
{
    typedef typename quadtree<T1, T2>::data_type data_type;
    
    os << "bounds: min=" << qt.bounds().min() << "--max=" << qt.bounds().max() << "\n";
    os << "depth=" << qt.depth() << " (max depth=" << qt.max_depth() << "), " << (qt.is_leaf() ? "LEAF" : "NODE") << '\n';
    os << "points (max size= " << qt.max_bucket_size() << "): ";
    const std::vector<data_type>& data = qt.data();
    if (!data.size()) {
        os << "<empty>";
    } else {
        for (int i = 0 ; i < data.size() - 1 ; ++i) {
            os << data[i] << ", ";
        }
        os << data.back();
    }
    
    return os;
}

}

#endif














