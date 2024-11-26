#pragma once

#include <assert.h>
// STL
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <memory>
#include <math/types.hpp>
#include <math/data_traits.hpp>

namespace spurt
{
template <typename Value>
struct bounding_box
{
    typedef data_traits<Value> value_traits; 

    typedef Value value_type;
    static const size_t dimension = value_traits::size();
    typedef typename value_traits::value_type scalar_type;

    bounding_box() 
        : m_min(value_traits::zero()), 
          m_max(value_traits::zero()) {}

    template<typename Vec>
    bounding_box(const Vec &min, const Vec &max) 
    {
        for (size_t i=0; i<dimension; ++i) {
            value_traits::value(m_min, i) = data_traits<Vec>::value(min, i);
            value_traits::value(m_max, i) = data_traits<Vec>::value(max, i);
        }
    }

    value_type size() const {
        return value_traits::minus(m_max, m_min); 
    }
    scalar_type diameter() const { return value_traits::norm(this->size()); }
    bool empty() const {
        return diameter() == 0;
    }

    template<typename Array>
    void add(const Array &p)
    {
        typedef data_traits<Array> array_traits;
        typedef typename array_traits::value_type other_scalar_type;
        for (auto i=0; i<dimension; ++i) {
            const other_scalar_type& v = array_traits::value(p, i);
            if (v < value_traits::value(m_min, i))
                value_traits::value(m_min, i) = v;
            if (v > value_traits::value(m_max, i))
                value_traits::value(m_max, i) = v;
        }
    }

    template< typename Array >
    bool inside(const Array &p) const
    {
        typedef data_traits<Array> array_traits;
        typedef typename array_traits::value_type other_scalar_type;
        for (auto i = 0; i < dimension; ++i)
        {
            const other_scalar_type& v = array_traits::value(p, i);
            if (v < value_traits::value(m_min, i)|| 
                v > value_traits::value(m_max, i))
                return false;
        }
        return true;
    }
    
    template <typename Array >
    scalar_type inf_distance(const Array& p) const {
        typedef data_traits<Array> array_traits;
        typedef typename array_traits::value_type other_scalar_type;
        scalar_type r = 0;
        for (auto i = 0; i < dimension; ++i)
        {
            const other_scalar_type& v = array_traits::value(p, i);
            const scalar_type& _min = value_traits::value(m_min, i);
            const scalar_type& _max = value_traits::value(m_max, i);
            if (v < _min) r = std::max(r, _min-v);
            else if (v > _max) r = std::max(r, v-_max);
        }
        return r;
    }

    const value_type& min() const { return m_min; }
    const value_type& max() const { return m_max; }
    value_type &min() { return m_min; }
    value_type &max() { return m_max; }

    value_type center() const { 
        return value_traits::mult(value_traits::plus(min(), max()), 0.5);
    }

    value_type m_min, m_max;
};

template<typename Value>
std::ostream& operator<<(std::ostream& os, const bounding_box<Value>& bbox) {
    os << "[ " << bbox.min() << " -> " << bbox.max() << " ]";
    return os;
}

typedef bounding_box<vec2> bbox2;
typedef bounding_box<vec3> bbox3;
typedef bounding_box<vec4> bbox4;
typedef bounding_box<vec5> bbox5;
typedef bounding_box<vec6> bbox6;

typedef bounding_box<fvec2> fbbox2;
typedef bounding_box<fvec3> fbbox3;
typedef bounding_box<fvec4> fbbox4;
typedef bounding_box<fvec5> fbbox5;
typedef bounding_box<fvec6> fbbox6;

typedef bounding_box<ivec2> ibbox2;
typedef bounding_box<ivec3> ibbox3;
typedef bounding_box<ivec4> ibbox4;
typedef bounding_box<ivec5> ibbox5;
typedef bounding_box<ivec6> ibbox6;

} // namespace spurt