#ifndef __SPURT_MATH_VECTOR_MANIP_HPP__
#define __SPURT_MATH_VECTOR_MANIP_HPP__

#include <algorithm>
#include <cmath>
#include <cassert>

namespace spurt { namespace vector {
    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ normsq(const Vector_& v) {
        typedef Value_ value_type;
        size_t N=v.size();
        value_type n=static_cast<value_type>(0);
        for (size_t i=0; i<N; ++i) n+=v[i]*v[i];
        return n;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ norm(const Vector_& v) {
        return sqrt(normsq(v));
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ l1_norm(const Vector_& v) {
        typedef Value_ value_type;
        size_t N=v.size();
        value_type n=static_cast<value_type>(0);
        for (size_t i=0; i<N; ++i) n+=std::abs(v[i]);
        return n;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ distance(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        size_t N=v.size();
        assert(N == w.size());
        value_type n=static_cast<value_type>(0);
        for (size_t i=0; i<N; ++i) n+=(v[i]-w[i])*(v[i]-w[i]);
        return sqrt(n);
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ dot(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        size_t N=v.size();
        assert(N == w.size());
        value_type d=static_cast<value_type>(0);
        for (size_t i=0; i<N; ++i) d+=v[i]*w[i];
        return d;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Vector_ sum(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        size_t N=v.size();
        assert(N == w.size());
        Vector_ s;
        for (size_t i=0; i<N; ++i) {
            s[i]=v[i]+w[i];
        }
        return s;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Vector_ diff(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        size_t N=v.size();
        assert(N == w.size());
        Vector_ s;
        for (size_t i=0; i<N; ++i) {
            s[i]=v[i]-w[i];
        }
        return s;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ min(const Vector_& v) {
        return *std::min_element(&v[0], &v[v.size()]);
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ min_id(const Vector_& v) {
        return std::distance(&v[0], std::min_element(&v[0], &v[v.size()]));
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ max(const Vector_& v) {
        return *std::max_element(&v[0], &v[v.size()]);
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ max_id(const Vector_& v) {
        return std::distance(&v[0], std::max_element(&v[0], &v[v.size()]));
    }

    template<typename Vector1_, typename Vector2_,
             typename Value1_=typename Vector1_::value_type,
             typename Value2_=typename Vector2_::value_type>
    inline void copy(const Vector1_& src, Vector2_& dest, size_t shift1=0, size_t shift2=0) {
        const size_t N=src.size();
        assert(N-shift1 == dest.size()-shift2);
        for (size_t i=shift1; i<N; ++i) {
            dest[shift2+i]=src[i];
        }
    }

} // namespace vector
} // namespace spurt

#endif
