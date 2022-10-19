#ifndef __MISC_HPP__
#define __MISC_HPP__

#include <vector>
#include <list>
#include <boost/rational.hpp>
#include <math/fixed_vector.hpp>
#include <misc/misc_helper.hpp>

namespace xavier {
template<typename T>
inline void push_front(const T& val, std::vector<T>& _in)
{
    if (!_in.size()) {
        _in.push_back(val);
    } else if (nvis::all(val != _in[0])) {
        std::vector<T> _out;
        _out.reserve(_in.size() + 1);
        _out.push_back(val);
        std::copy(_in.begin(), _in.end(), std::back_inserter(_out));
        std::swap(_in, _out);
    }
}

template<typename T>
inline double value(const boost::rational<T>& r)
{
    return (double)r.numerator() / (double)r.denominator();
}

/*
template<typename T>
inline T sign(const T& t)
{
    return (t > 0 ? 1 : -1);
}
*/

template<typename Iterator>
inline double min_dist(double v, const Iterator& begin, const Iterator& end)
{
    std::list<double> dist;
    for (Iterator it = begin ; it != end ; ++it) {
        dist.push_back(fabs(*it - v));
    }
    return *std::min_element(dist.begin(), dist.end());
}

struct Edge {
    Edge() : i0(0), i1(0) {}
    Edge(size_t _i0, size_t _i1) : i0(_i0), i1(_i1) {
        if (i0 > i1) {
            size_t tmp = i0;
            i0 = i1;
            i1 = tmp;
        }
    }
    Edge(const Edge& e) : i0(e.i0), i1(e.i1) {}
    
    bool operator<(const Edge& e) const {
        if (i0 < e.i0) {
            return true;
        } else if (i0 > e.i0) {
            return false;
        }
        return i1 < e.i1;
    }
    
    size_t i0, i1;
};

} // xavier

#endif




