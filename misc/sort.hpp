#ifndef __SORT_HPP__
#define __SORT_HPP__

#include <map>
#include <vector>
#include <list>
#include <math/fixed_vector.hpp>

namespace spurt {
template<typename T1, typename T2>
struct Lt_pair_lexicographic {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1) const {
        return (p0.first < p1.first) ||
               (p0.first == p1.first && p0.second < p1.second);
    }
};

template<typename T1, typename T2>
struct Lt_pair_first {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1) const {
        return (p0.first < p1.first);
    }
};

template<typename T1, typename T2>
struct Lt_pair_second {
    bool operator()(const std::pair<T1, T2>& p0, const std::pair<T1, T2>& p1) const {
        return (p0.second < p1.second);
    }
};

template<typename V, int N>
struct Lt_vector_N {
    bool operator()(const V& v0, const V& v1) const {
        return (v0[N] < v1[N]);
    }
};

struct Lt_fabs {
    bool operator()(double a, double b) const {
        return fabs(a) < fabs(b);
    }
};

// A simple fixed length container that keeps its elements in the order
// defined by std::less<T>
template<typename T, int N, class Comp = std::less<T> >
class fixed_sorted_vector : public nvis::fixed_vector<T, N> {
public:
    typedef nvis::fixed_vector<T, N>     base_type;
    typedef fixed_sorted_vector<T, N>    self_type;
    typedef Comp                        less_than;

    fixed_sorted_vector() : base_type(), original() {}

    fixed_sorted_vector(const base_type& v) : base_type(v), original(v) {
        std::sort(this->begin(), this->end(), less_than());
    }

    fixed_sorted_vector(const T& v) : base_type(v), original(v) {}

    fixed_sorted_vector(const T& v0, const T& v1)
        : base_type(v0, v1), original(v0, v1) {
        std::sort(this->begin(), this->end(), less_than());
    }

    fixed_sorted_vector(const T& v0, const T& v1, const T& v2)
        : base_type(v0, v1, v2), original(v0, v1, v2) {
        std::sort(this->begin(), this->end(), less_than());
    }

    fixed_sorted_vector(const T& v0, const T& v1, const T& v2, const T& v3)
        : base_type(v0, v1, v2, v3) , original(v0, v1, v2, v3){
        std::sort(this->begin(), this->end(), less_than());
    }

    fixed_sorted_vector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4)
        : base_type(v0, v1, v2, v3, v4), original(v0, v1, v2, v3, v4) {
        std::sort(this->begin(), this->end(), less_than());
    }

    template<typename V>
    fixed_sorted_vector(const V& v) : base_type(&v[0], &v[N]), original(&v[0], &v[N]) {
        std::sort(this->begin(), this->end(), less_than());
    }

    bool operator==(const self_type& sa) const {
        return nvis::all(static_cast<base_type>(*this) ==
                         static_cast<base_type>(sa));
    }

    bool operator<(const self_type& sa) const {
        nvis::lexicographical_order _order;
        return _order(static_cast<base_type>(*this),
                      static_cast<base_type>(sa));
    }

    const T& min() const {
        return (*this)[0];
    }

    const T& max() const {
        return (*this)[N-1];
    }

    base_type original;
};

template< typename T, typename Size_=size_t >
void sort(const std::vector< T >& values,
          std::vector< Size_ >& sorted);

template< typename T, typename Order, typename Size_=size_t >
void sort(const std::vector< T >& values,
          std::vector< Size_ >& sorted,
          const Order& order);

// template< typename T >
// void sort(const std::vector< T >::const_iterator& _begin,
//           const std::vector< T >::const_iterator& _end,
//           std::vector< unsigned int >& sorted);

template< typename T, typename Size_ >
bool test(const std::vector< T >& values,
          const std::vector< Size_ >& sorted);
};

// template< typename T >
// inline void spurt::sort(const T* _begin,
//                          const T* _end,
//                          std::vector< unsigned int >& sorted)
// {
//     sorted.clear();
//     std::map< T, unsigned int > _map;
//     std::vector< T >::const_iterator itv;
//     for (itv = _begin ; itv != _end ; itv++)
//         _map[*itv] = i;
//
//     typename std::map< T, unsigned int >::const_iterator itm;
//     for (itm = _map.begin() ; itm != _map.end() ; itm++)
//         sorted.push_back(itm->second);
// }

template< typename T, typename Size_ >
inline void spurt::sort(const std::vector< T >& values,
                         std::vector< Size_ >& sorted)
{
    sorted.clear();
    typedef std::list< Size_ >                        list_type;
    typename list_type::const_iterator                lit;
    std::map< T, list_type >                          _map;
    typename std::map< T, list_type >::iterator       mit;
    typename std::map< T, list_type >::const_iterator cmit;
    for (Size_ i = 0 ; i < values.size() ; i++) {
        mit = _map.find(values[i]);
        if (mit == _map.end())
            _map[values[i]].push_back(i);
        else {
            mit->second.push_back(i);
        }
    }

    for (cmit = _map.begin() ; cmit != _map.end() ; ++cmit) {
        for (lit = cmit->second.begin() ; lit != cmit->second.end(); ++lit)
            sorted.push_back(*lit);
    }
}

template< typename T, typename Order, typename Size_ >
inline void spurt::sort(const std::vector< T >& values,
                         std::vector< Size_ >& sorted,
                         const Order& order)
{
    sorted.clear();
    typedef std::list< Size_ >                               list_type;
    typename list_type::const_iterator                       lit;
    std::map< T, list_type, Order >                          _map;
    typename std::map< T, list_type, Order >::iterator       mit;
    typename std::map< T, list_type, Order >::const_iterator cmit;
    for (Size_ i = 0 ; i < values.size() ; i++) {
        mit = _map.find(values[i]);
        if (mit == _map.end())
            _map[values[i]].push_back(i);
        else {
            mit->second.push_back(i);
        }
    }

    for (cmit = _map.begin() ; cmit != _map.end() ; ++cmit) {
        for (lit = cmit->second.begin() ; lit != cmit->second.end(); ++lit)
            sorted.push_back(*lit);
    }
}

template< typename T, typename Size_ >
bool spurt::test(const std::vector< T >& values,
                  const std::vector< Size_ >& sorted)
{
    if (!sorted.size() && !values.size()) return true;

    T old, cur;
    old = values[sorted[0]];
    for (Size_ i = 1 ; i < sorted.size() ; i++) {
        cur = values[sorted[i]];
        if (old > cur) return false;

        old = cur;
    }

    return true;
}


#endif
