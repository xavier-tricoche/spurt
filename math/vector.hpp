#ifndef __SPURT_MATH_VECTOR_HPP__
#define __SPURT_MATH_VECTOR_HPP__

#include <algorithm>
#include <cmath>

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
        for (size_t i=0; i<N; ++i) n+=(v[i]-w[i])*(v[i]*w[i]);
        return sqrt(n);
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ dot(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        size_t N=v.size();
        assert(N == w.size());
        value_type d=static_cast<value_type>(0);
        for (size_t i=0; i<N; ++i) d+=v[i]*w[i];
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
    inline Value_ max(const Vector_& v) {
        return *std::max_element(&v[0], &v[v.size()]);
    }
} }

#endif
