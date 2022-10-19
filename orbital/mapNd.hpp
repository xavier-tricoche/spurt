#ifndef __mapNd_hpp
#define __mapNd_hpp

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <vector>

template< int N >
struct mapNd {

    typedef nvis::fixed_vector< double, N >     vecN;
    
    struct map_undefined {
    };
    
    virtual vecN map(const vecN& in, int niter) const = 0;
    virtual void map(const vecN& in, std::vector<vecN>& out, int niter) const = 0;
    virtual void map(const vecN& seed, std::vector<vecN>& steps, std::vector< double >& times, int niter) const = 0;
    
    virtual nvis::bounding_box<vecN> bounds() const = 0;
    
    virtual mapNd* clone() const = 0;
};

#endif // __map2d_hpp

