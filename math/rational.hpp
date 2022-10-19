#ifndef __RATIONAL_HPP__
#define __RATIONAL_HPP__

#include <vector>
#include <list>
#include <boost/rational.hpp>

namespace xavier {

template<typename I, typename F>
inline F value(const boost::rational<I>& r)
{
    return (F)r.numerator() / (F)r.denominator();
}

template<typename I, typename F>
inline boost::rational<I> rational_approx(F v, I maxnum, I maxden =0)
{
    // naive method
    typedef I                       _int;
    typedef F                       _float;
    typedef boost::rational<_int>   _rational;
    if (!maxden) {
        maxden = maxnum;
    }
    std::map<_float, _rational> approx;
    for (_int num=1 ; num<=maxnum ; ++num) {
        for (_int den=1 ; den<=maxden ; ++den) {
            _rational r(num, den);
            _float err = fabs(v-value<_int, _float>(r));
            approx.insert(std::pair<_float, _rational>(err, r));
        }
    }
    return approx.begin()->second;
}

// template<typename I, typename F>
// inline boost::rational<I> approximate(F v, I num) {
//  typedef boost::rational<I> rational_type;
//  std::map<F, rational_type> approx;
// }

}



#endif