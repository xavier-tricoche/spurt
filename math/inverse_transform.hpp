#ifndef __inverse_transform_sampling__
#define __inverse_transform_sampling__

#include <iostream>
#include <vector>
#include <math/math.hpp>

namespace xavier {

template<typename T>
class inverse_transform_sampling {
public:

    typedef T   value_type;
    
    inverse_transform_sampling(const std::vector<value_type>& values) {
        size_t N = values.size();
        sort_ids(sorted_ids, values);
        cdf.resize(N);
        value_type sum(0);
        value_type min_val = values[sorted_ids.front()];
        value_type max_val = values[sorted_ids.back()];
        for (int i = 0 ; i < N ; ++i) {
            sum += values[sorted_ids[i]] - min_val;
            cdf[i] = sum;
        }
        for (int i = 0 ; i < N ; ++i) {
            cdf[i] /= sum;
        }
        srand48(time(0));
    }
    
    unsigned int sample() const {
        double u = drand48();
        size_t pos = std::distance(cdf.begin(), std::lower_bound(cdf.begin(), cdf.end(), u));
        return sorted_ids[pos];
    }
    
private:
    std::vector<unsigned int>   sorted_ids;
    std::vector<value_type>     cdf;
};

}



#endif



