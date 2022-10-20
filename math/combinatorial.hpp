#ifndef __XAVIER_COMBINATORIAL_HPP__
#define __XAVIER_COMBINATORIAL_HPP__

#include <algorithm>
#include <cassert>

/*
    Combinatorial namespace: simple functions needed to compute the cardinal
    of subsets that arise in combinatorics problems.
 */

namespace spurt { namespace combinatorial {
    
template<typename _Int>
inline _Int factorial(_Int hi, _Int lo=1) {
    typedef _Int number_type;
    
    assert(hi>=0 && lo>=0);
    // Note: (hi < lo) is considered valid and returns 1 
    // (needed by binomial below)
    
    number_type c=lo, r=1;
    while (c <= hi) {
        r *= c++;
    }
    return r;
}

template<typename _Int>
inline _Int binomial(_Int n, _Int k) {
    typedef _Int                                number_type;
    typedef std::pair<number_type, number_type> pair_type;
    
    assert(n>=0 && k>=0 && n>=k);
    
    pair_type lohi = std::minmax(k, n-k);
    return factorial(n, lohi.second + 1)/factorial(lohi.first);
}

// compute multinomial coefficient for subsets of cardinal k
// chosen from a set of cardinal n
template<typename _Int>
inline _Int multinomial(_Int n, _Int k) {
    typedef _Int number_type;
    
    assert(n>0 && k>0);
    
    return binomial(n+k-1, k);
}


} // combinatorial

} // combinatorial


#endif