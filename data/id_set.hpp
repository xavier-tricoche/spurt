#ifndef __ID_SET_HPP__
#define __ID_SET_HPP__

#include <math/fixed_vector.hpp>
#include <boost/static_assert.hpp>

namespace spurt
{
template<typename T, int N>
struct id_set {
    typedef T                           value_type;
    typedef fixed_vector<T, N>    vec_type;

    id_set() : __ids(-1) {}
    id_set(const value_type ids[N]) {
        for (int i = 0 ; i < N ; ++i) __ids[i] = ids[i];
        std::sort(__ids.begin(), __ids.end());
    }
    id_set(value_type i0, value_type i1) {
        BOOST_STATIC_ASSERT(N == 2);
        if (i0 < i1) {
            __ids[0] = i0;
            __ids[1] = i1;
        }
        else {
            __ids[0] = i1;
            __ids[1] = i0;
        }
    }
    id_set(value_type i0, value_type i1, value_type i2) {
        BOOST_STATIC_ASSERT(N == 3);
        __ids[0] = i0;
        __ids[1] = i1;
        __ids[2] = i2;
        std::sort(__ids.begin(), __ids.end());
    }
    id_set(value_type i0, value_type i1, value_type i2, value_type i3) {
        BOOST_STATIC_ASSERT(N == 4);
        __ids[0] = i0;
        __ids[1] = i1;
        __ids[2] = i2;
        __ids[3] = i3;
        std::sort(__ids.begin(), __ids.end());
    }
    id_set(const vec_type& ids) : __ids(ids) {
        std::sort(__ids.begin(), __ids.end());
    }

    vec_type __ids;
};

struct Lt_id_set {
    template<typename T, int N>
    bool operator()(const id_set<T, N>& f0, const id_set<T, N>& f1) {
        lexicographical_order Lt;
        return Lt(f0.__ids, f1.__ids);
    }
};

typedef id_set<int, 2>    id2;
typedef id_set<int, 3>    id3;
typedef id_set<int, 4>    id4;

} // xavier
