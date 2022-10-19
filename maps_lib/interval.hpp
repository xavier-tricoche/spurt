#ifndef __INTERVAL_HPP__
#define __INTERVAL_HPP__

#include <iostream>

namespace xavier {
template<typename T>
struct interval {
    interval() : __min(T(0)), __max(T(0)) {}
    interval(T min, T max) : __min(min), __max(max) {}
    interval(const interval<T>& _int) : __min(_int.__min), __max(_int.__max) {}
    
    bool inside(T val) const {
        return (val >= __min && val <= __max);
    }
    
    bool empty() const {
        return (__max <= __min);
    }
    
    T length() const {
        return std::max(0., __max - __min);
    }
    
    T __min, __max;
};
}

template<typename T>
inline xavier::interval<T> intersect(const xavier::interval<T>& i0, const xavier::interval<T>& i1)
{
    return xavier::interval<T>(std::max(i0.__min, i1.__min),
                               std::min(i0.__max, i1.__max));
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const xavier::interval<T>& i)
{
    os << "[ " << i.__min << ", " << i.__max << " ]";
    return os;
}


#endif


