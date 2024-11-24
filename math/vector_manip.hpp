#ifndef __SPURT_MATH_VECTOR_MANIP_HPP__
#define __SPURT_MATH_VECTOR_MANIP_HPP__

#include <algorithm>
#include <cmath>
#include <cassert>

namespace spurt { namespace vector {
    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ normsq(const Vector_& v) {
        typedef Value_ value_type;
        value_type n=static_cast<value_type>(0);
        std::for_each(v.begin(), v.end(), [&](const value_type& x){
           n += x*x; 
        });
        return n;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ norm(const Vector_& v) {
        return sqrt(normsq(v));
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ l1_norm(const Vector_& v) {
        typedef Value_ value_type;
        value_type n=static_cast<value_type>(0);
        std::for_each(v.begin(), v.end(), [&](const value_type& x){
           n += std::abs(x); 
        });
        return n;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ distance(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        value_type n=static_cast<value_type>(0);
        auto wit = w.begin();
        for (auto vit=v.begin(); vit!=v.end(); ++vit, ++wit) {
            n += ((*vit)-(*wit))*((*vit)-(*wit));
        }
        return sqrt(n);
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ dot(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        auto wit = w.begin();
        value_type d=static_cast<value_type>(0);
        for (auto vit=v.begin(); vit!=v.end(); ++vit, ++wit)
        {
            d += (*vit) * (*wit);
        }
        return d;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Vector_ sum(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        Vector_ s;
        auto wit = w.begin();
        auto sit = w.begin();
        for (auto vit=v.begin(); vit!=v.end() && wit!=w.end(); ++vit, ++wit, ++sit) {
            *sit=(*vit)+(*wit);
        }
        return s;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Vector_ diff(const Vector_& v, const Vector_& w) {
        typedef Value_ value_type;
        Vector_ s;
        auto wit = w.begin();
        auto sit = w.begin();
        for (auto vit=v.begin(); vit!=v.end(); ++vit, ++wit, ++sit) {
            *sit=(*vit)-(*wit);
        }
        return s;
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ min(const Vector_& v) {
        return *std::min_element(v.begin(), v.end());
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ min_id(const Vector_& v) {
        return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ max(const Vector_& v) {
        return *std::max_element(v.begin(), v.end());
    }

    template<typename Vector_, typename Value_=typename Vector_::value_type>
    inline Value_ max_id(const Vector_& v) {
        return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
    }

    template<typename Vector1_, typename Vector2_,
             typename Value1_=typename Vector1_::value_type,
             typename Value2_=typename Vector2_::value_type>
    inline void copy(const Vector1_& src, Vector2_& dest, size_t shift1=0, size_t shift2=0) {
        auto s1 = std::distance(src.begin(), src.end());
        auto s2 = std::distance(dest.begin(), dest.end())
        assert(s1-shift1 == s2-shift2);
        for (size_t i=shift1; i<s1; ++i) {
            dest[shift2+i]=src[i];
        }
    }

} // namespace vector
} // namespace spurt

#endif
