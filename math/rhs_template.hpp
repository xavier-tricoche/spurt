#ifndef __AUTOMATICALLY_GENERATED_HEADER_HPP__
#define __AUTOMATICALLY_GENERATED_HEADER_HPP__

#include <math/fixed_vector.hpp>

namespace custom_rhs {

template<typename T>
struct RHS {
    typedef T                        value_type;
    typedef nvis::fixed_vector<T, 3> vector_type;
    
    vector_type operator()(const vector_type& x, value_type t=0) const {
        vector_type rhs;
        /* insert expression here */
        <REPLACEME>
        /* end of expression */
        return rhs;
    }
};

}
