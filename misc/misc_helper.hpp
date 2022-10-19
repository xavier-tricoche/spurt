#ifndef __XAVIER_MISC_HELPER_HPP__
#define __XAVIER_MISC_HELPER_HPP__

namespace xavier {

template<typename T>
inline T sign(T a)
{
    return (a < 0) ? -1 : 1;
}

template<typename T>
inline T mindist(T a, T b)
{
    if (a > 0.5*b) {
        return a - b;
    } else if (a < -0.5*b) {
        return a + b;
    } else {
        return a;
    }
}

template<typename T>
inline T minimum(const T& t1, const T& t2) { return std::min<T>(t1, t2); }
template<typename T>
inline T minimum(const T& t1, const T& t2, const T& t3) {
    return std::min(std::min(t1, t2), t3);
}
template<typename T>
inline T minimum(const T& t1, const T& t2, const T& t3, const T& t4) {
    return std::min(std::min(t1, t2), std::min(t3, t4));
}
template<typename T> inline T 
minimum(const T& t1, const T& t2, const T& t3, const T& t4, const T& t5) {
    return std::min(minimum(t1, t2, t3), std::min(t4, t5));
}
template<typename T>
inline T minimum(const std::vector<T>& ts) {
    return *std::min_element(ts.begin(), ts.end());
}

template<typename T>
inline T maximum(const T& t1, const T& t2) { return std::max<T>(t1, t2); }
template<typename T>
inline T maximum(const T& t1, const T& t2, const T& t3) {
    return std::max(std::max(t1, t2), t3);
}
template<typename T>
inline T maximum(const T& t1, const T& t2, const T& t3, const T& t4) {
    return std::max(std::max(t1, t2), std::max(t3, t4));
}
template<typename T> inline T 
maximum(const T& t1, const T& t2, const T& t3, const T& t4, const T& t5) {
    return std::max(maximum(t1, t2, t3), std::max(t4, t5));
}
template<typename T>
inline T maximum(const std::vector<T>& ts) {
    return *std::max_element(ts.begin(), ts.end());
}

} // namespace xavier

#endif
