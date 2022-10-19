#include <kdtree++/kdtree.hpp>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

class box_point {
public:
    typedef nvis::vec3      vector_type;
    typedef nvis::bbox3     box_type;
    typedef double          value_type;
    
    box_point() : __idx(-1) {}
    box_point(const bbox_type& box, size_t idx = -1) : _b(box), __idx(idx) {}
    box_point(const box_point& bp) : __b(bp.__b), __idx(bp.__idx) {}
    
    const box_type& box() const {
        return __b;
    }
    
    size_t index() const {
        return __idx;
    }
    
    value_type distance_to(const box_point& p) const {
    
    }
    
private:
    nvis::bbox3 __b;
    size_t      __idx;
};



