#ifndef __rhs_base_hpp
#define __rhs_base_hpp

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <stdexcept>

template< int N >
class rhs_base {
public:

    typedef nvis::fixed_vector<double, N>   vecN;
    typedef nvis::bounding_box<vecN>        bboxN;
    
    struct undefined_point : public std::runtime_error {
        undefined_point() : std::runtime_error( "undefined point" ) {
        }
    };
    
    virtual vecN operator()( const double& t, const vecN& y ) const = 0;
    
    virtual nvis::fixed_vector<double, N+1> plane() const = 0;
    
    virtual rhs_base* clone() const = 0;
    
    virtual bboxN bounds() const = 0;
};

#endif // __rhs_base_hpp
