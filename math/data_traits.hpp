#pragma once

#include <type_traits>
#include <array>
#include <vector>
#include <Eigen/Core>
#include <math/types.hpp>
#include <numeric>

namespace spurt {

// uniform treatment of data attribute as (possibly unary) array
template <typename T, typename Enable = void>
struct data_traits {};

template <typename T>
struct data_traits<T, typename std::enable_if<std::is_scalar<T>::value>::type>
{
    typedef T value_type;
    typedef T data_type;

    constexpr static size_t size() { return 1; };
    constexpr static size_t ncols() { return 1; }
    constexpr static size_t nrows() { return 1; }
    constexpr static size_t mem_size() { return sizeof(value_type); }
    static const value_type &value(const data_type &v, size_t i = 0, size_t j = 0)
    {
        return v;
    }
    static value_type &value(data_type &v, size_t i = 0, size_t j = 0) { return v; }
    static data_type &assign(data_type &inout, value_type val)
    {
        inout = val;
        return inout;
    }
    static value_type norm(const data_type &v) { return std::labs(v); }
    static data_type zero() { return static_cast<value_type>(0); }
    static data_type one() { return static_cast<value_type>(1); }
    static data_type minus(const data_type &a, const data_type &b)
    {
        return a - b;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        return a + b;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        return a * b;
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        return a / b;
    }
    static data_type mult(const data_type &a, value_type b)
    {
        return a * b;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return a / b;
    }
};

// template<typename T, size_t N>
// struct data_traits< nvis::fixed_vector<T, N> > {
//     typedef T value_type;
//     typedef nvis::fixed_vector<value_type, N> data_type;

//     constexpr static size_t size() { return N; }
//     constexpr static size_t ncols() { return 1; }
//     constexpr static size_t nrows() { return N; }
//     constexpr static size_t mem_size() { return N*sizeof(value_type); }
//     static const value_type& value(const data_type& v, size_t i, size_t j=0) {
//         return v[i];
//     }
//     static value_type& value(data_type& v, size_t i, size_t j=0) {
//         return v[i]; }
//     static data_type& assign(data_type& inout, T val) {
//         std::fill(inout.begin(), inout.end(), val);
//         return inout;
//     }
//     static value_type norm(const data_type& v) { return nvis::norm(v); }
//     static data_type zero() { return data_type(0); }
//     static data_type one() { return data_type(1); }
//     static data_type minus(const data_type& a, const data_type& b) {
//         return a - b;
//     }
//     static data_type plus(const data_type &a, const data_type &b)
//     {
//         return a + b;
//     }
//     static data_type dotprod(const data_type& a, const data_type& b) {
//         return a*b;
//     }
//     static data_type dotdiv(const data_type &a, const data_type &b)
//     {
//         return a / b;
//     }
//     static data_type mult(const data_type& a, value_type b) {
//         return b*a;
//     }
//     static data_type div(const data_type& a, value_type b) {
//         return a/b;
//     }
// };

template <typename T, size_t N>
struct data_traits< std::array<T, N> >
{
    typedef T value_type;
    typedef std::array<value_type, N> data_type;

    constexpr static size_t size() { return N; }
    constexpr static size_t ncols() { return 1; }
    constexpr static size_t nrows() { return N; }
    constexpr static size_t mem_size() { return N * sizeof(value_type); }
    static const value_type &value(const data_type &v, size_t i, size_t j = 0)
    {
        return v[i];
    }
    static value_type &value(data_type &v, size_t i, size_t j = 0)
    {
        return v[i];
    }
    static data_type &assign(data_type &inout, value_type val)
    {
        std::fill(inout.begin(), inout.end(), val);
        return inout;
    }
    static value_type norm(const data_type &v)
    {
        return std::inner_product(v.begin(), v.end(), v.begin(), 0);
    }
    static data_type zero()
    {
        data_type a;
        std::fill(a.begin(), a.end(), static_cast<value_type>(0));
        return a;
    }
    static data_type one()
    {
        data_type a;
        std::fill(a.begin(), a.end(), static_cast<value_type>(1));
        return a;
    }
    static data_type minus(const data_type &a, const data_type &b)
    {
        data_type r;
        for (int i = 0; i < N; ++i)
            r[i] = a[i] - b[i];
        return r;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        data_type r;
        for (int i = 0; i < N; ++i)
            r[i] = a[i] + b[i];
        return r;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        data_type out;
        for (int i = 0; i < N; ++i)
            out[i] = a[i] * b[i];
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        data_type out;
        for (int i = 0; i < N; ++i)
            out[i] = a[i] / b[i];
    }
    static data_type mult(const data_type &a, value_type b)
    {
        data_type r;
        for (int i = 0; i < N; ++i)
            r[i] = a[i] * b;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return mult(a, 1. / b);
    }
};

// specialization for generic Eigen matrix type
template <typename Derived>
struct data_traits<Eigen::MatrixBase<Derived>>
{

    typedef Eigen::MatrixBase<Derived> data_type;
    typedef typename data_type::Scalar value_type;

    constexpr static size_t size() { return ncols() * nrows(); }
    constexpr static size_t ncols() { 
        if (data_type::ColsAtCompileTime == Eigen::Dynamic) {
            static_assert(false, "ERROR: trying to access static size of dynamic matrix.");
        }
        return data_type::ColumnsAtCompileTime; 
    }
    constexpr static size_t nrows() {
        if (data_type::RowsAtCompileTime == Eigen::Dynamic) {
            static_assert(false, "ERROR: trying to access static size of dynamic matrix.");
        }
        else return data_type::RowsAtCompileTime;
    }
    constexpr static size_t mem_size()
    {
        return size() * sizeof(value_type);
    }
    static const value_type &value(const data_type &v, size_t i, size_t j = 0)
    {
        return v(i, j);
    }
    static value_type &value(data_type &v, size_t i, size_t j = 0)
    {
        return v(i, j);
    }
    static data_type &assign(data_type &inout, value_type val)
    {
        inout.setConstant(val);
        return inout;
    }
    static value_type norm(const data_type &v)
    {
        return v.norm();
    }
    static data_type zero() { return data_type::Zero(); }
    static data_type one() { return data_type::One(); }
    static data_type minus(const data_type &a, const data_type &b)
    {
        return a - b;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        return a + b;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        return (a.array() * b.array()).matrix();
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        return (a.array() / b.array()).matrix();
    }
    static data_type mult(const data_type &a, value_type b)
    {
        return b * a;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return 1 / b * a;
    }
};

template <typename BinaryOp, typename LhsType, typename RhsType>
struct data_traits< Eigen::CwiseBinaryOp<BinaryOp, LhsType, RhsType> >
{
    typedef LhsType data_type;
    typedef typename data_type::Scalar value_type;
    typedef Eigen::CwiseBinaryOp<BinaryOp, LhsType, RhsType> op_type;

    constexpr static size_t size() { return op_type::size(); }
    constexpr static size_t ncols() { return op_type::cols(); }
    constexpr static size_t nrows() { return op_type::rows(); }
    constexpr static size_t mem_size()
    {
        return size() * sizeof(value_type);
    }
    static const value_type &value(const op_type &v, size_t i, size_t j = 0)
    {
        static_assert(false, "cannot access coefficients of binary expression");
        return value_type();
    }
    static value_type &value(data_type &v, size_t i, size_t j = 0)
    {
        static_assert(false, "cannot modify coefficients of binary expression");
        return value_type();
    }
    static data_type &assign(data_type &inout, value_type val)
    {
        static_assert(false, "cannot modify coefficients of binary expression");
        inout.setConstant(val);
        return inout;
    }
    static value_type norm(const data_type &v)
    {
        static_assert(false, "cannot compute norm of binary expression");
        return v.norm();
    }
    static data_type zero() { return data_type::Zero(); }
    static data_type one() { return data_type::One(); }
    static data_type minus(const data_type &a, const data_type &b)
    {
        return a - b;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        return a + b;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        return (a.array() * b.array()).matrix();
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        return (a.array() / b.array()).matrix();
    }
    static data_type mult(const data_type &a, value_type b)
    {
        return b * a;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return 1 / b * a;
    }
};

template <typename T, int M, int N>
struct data_traits<Eigen::Matrix<T, M, N>>
{
    typedef Eigen::Matrix<T, M, N> data_type;
    typedef T value_type;

    constexpr static size_t size() { return ncols() * nrows(); }
    constexpr static size_t ncols() { return M; }
    constexpr static size_t nrows() { return N; }
    constexpr static size_t mem_size()
    {
        return size() * sizeof(value_type);
    }
    static const value_type &value(const data_type &v, size_t i, size_t j = 0)
    {
        return v(i, j);
    }
    static value_type &value(data_type &v, size_t i, size_t j = 0)
    {
        return v(i, j);
    }
    static data_type &assign(data_type &inout, value_type val)
    {
        inout.setConstant(val);
        return inout;
    }
    static value_type norm(const data_type &v)
    {
        return v.norm();
    }
    static data_type zero() { return data_type::Zero(); }
    static data_type one() { return data_type::One(); }
    static data_type minus(const data_type &a, const data_type &b)
    {
        return a - b;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        return a + b;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        return (a.array() * b.array()).matrix();
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        return (a.array() / b.array()).matrix();
    }
    static data_type mult(const data_type &a, value_type b)
    {
        return b * a;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return 1 / b * a;
    }
};

template<typename T, size_t N>
struct data_traits<small_vector<T, N>>
{
    typedef small_vector<T, N> data_type;
    typedef T value_type;

    constexpr static size_t size() { return N; }
    constexpr static size_t ncols() { return 1; }
    constexpr static size_t nrows() { return N; }
    constexpr static size_t mem_size()
    {
        return size() * sizeof(value_type);
    }
    static const value_type &value(const data_type &v, size_t i, size_t j=0)
    {
        return v[i];
    }
    static value_type &value(data_type &v, size_t i, size_t j=0)
    {
        return v[i];
    }
    static data_type &assign(data_type &inout, value_type val)
    {
        inout = data_type(val);
    }
    static value_type norm(const data_type &v)
    {
        return spurt::norm(v);
    }
    static data_type zero() { return data_type(0); }
    static data_type one() { return data_type(1); }
    static data_type minus(const data_type &a, const data_type &b)
    {
        return a - b;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        return a + b;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        return spurt::inner(a, b);
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        return a / b;
    }
    static data_type mult(const data_type &a, value_type b)
    {
        return b * a;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return a / b;
    }
};

template<typename T, size_t M, size_t N>
struct data_traits<small_matrix<T, M, N>>
{
    typedef small_matrix<T, M, N> data_type;
    typedef T value_type;

    constexpr static size_t size() { return data_type::ncols * data_type::nrows; }
    constexpr static size_t ncols() { return M; }
    constexpr static size_t nrows() { return N; }
    constexpr static size_t mem_size()
    {
        return size() * sizeof(value_type);
    }
    static const value_type &value(const data_type &v, size_t i, size_t j = 0)
    {
        return v(i, j);
    }
    static value_type &value(data_type &v, size_t i, size_t j = 0)
    {
        return v(i, j);
    }
    static data_type &assign(data_type &inout, value_type val)
    {
        inout = data_type(val);
    }
    static value_type norm(const data_type &v)
    {
        return spurt::norm(v);
    }
    static data_type zero() { return data_type(0); }
    static data_type one() { return data_type(1); }
    static data_type minus(const data_type &a, const data_type &b)
    {
        return a - b;
    }
    static data_type plus(const data_type &a, const data_type &b)
    {
        return a + b;
    }
    static data_type dotprod(const data_type &a, const data_type &b)
    {
        return spurt::inner(a, b);
    }
    static data_type dotdiv(const data_type &a, const data_type &b)
    {
        return a / b;
    }
    static data_type mult(const data_type &a, value_type b)
    {
        return b * a;
    }
    static data_type div(const data_type &a, value_type b)
    {
        return a / b;
    }
};

template <typename T>
struct array_wrapper
{
    typedef data_traits<T> traits_type;
    typedef typename traits_type::value_type value_type;
    typedef T data_type;
    typedef value_type *iterator;
    typedef const value_type *const_iterator;
    typedef array_wrapper<T> self_type;

    array_wrapper(data_type &v) : m_v(v) {}
    int size() const { return traits_type::size(m_v); }
    constexpr static size_t mem_size() { return traits_type::mem_size(); }
    const value_type &operator[](size_t i) const
    {
        return traits_type::value(m_v, i);
    }
    const value_type &operator()(size_t i, size_t j = 0) const
    {
        return traits_type::value(m_v, i, j);
    }
    value_type &operator[](size_t i) { return traits_type::value(m_v, i); }
    value_type &operator()(size_t i, size_t j = 0)
    {
        return traits_type::value(m_v, i, j);
    }
    self_type &operator=(const value_type &v)
    {
        traits_type::assign(m_v, v);
        return *this;
    }
    self_type &operator=(const self_type &other)
    {
        m_v = other.m_v;
        return *this;
    }
    data_type &assign(value_type s)
    {
        traits_type::assign(m_v, s);
        return m_v;
    }
    data_type &set(value_type s) { return this->assign(s); }
    value_type norm() const { return traits_type::norm(m_v); }
    const_iterator begin() const { return &traits_type::value(m_v, 0); }
    const_iterator end() const
    {
        return &traits_type::value(m_v, 0) + this->size();
    }
    self_type operator+(const self_type &other) const
    {
        data_type r = traits_type::plus(m_v, other.m_v);
        return self_type(r);
    }
    self_type operator-(const self_type &other) const
    {
        data_type r = traits_type::minus(m_v, other.m_v);
        return self_type(r);
    }
    self_type dotprod(const self_type &other) const
    {
        return self_type(traits_type::dotprod(m_v, other.m_v));
    }

    data_type &m_v;
};

template <typename T>
struct const_array_wrapper
{
    typedef data_traits<T> traits_type;
    typedef typename traits_type::value_type value_type;
    typedef T data_type;
    typedef const value_type *const_iterator;
    typedef const_array_wrapper<T> self_type;

    const_array_wrapper(const data_type &v) : m_v(v) {}
    int size() const { return traits_type::size(m_v); }
    constexpr static size_t mem_size() { return traits_type::mem_size(); }
    const value_type &operator[](size_t i) const
    {
        return traits_type::value(m_v, i);
    }
    const value_type &operator()(size_t i, size_t j = 0) const
    {
        return traits_type::value(m_v, i, j);
    }
    value_type norm() const { return traits_type::norm(m_v); }
    const_iterator begin() const { return &traits_type::value(m_v, 0); }
    const_iterator end() const
    {
        return &traits_type::value(m_v, 0) + this->size();
    }
    self_type operator+(const self_type &other) const
    {
        data_type r = traits_type::plus(m_v, other.m_v);
        return self_type(r);
    }
    self_type operator-(const self_type &other) const
    {
        data_type r = traits_type::minus(m_v, other.m_v);
        return self_type(r);
    }
    self_type dotprod(const self_type &other) const
    {
        return self_type(traits_type::dotprod(m_v, other.m_v));
    }

    const data_type &m_v;
};

template <typename ArrayOut, typename ArrayIn>
ArrayOut as(const ArrayIn &in)
{
    ArrayOut out;
    typedef data_traits<ArrayOut> out_traits_t;
    typedef data_traits<ArrayIn> in_traits_t;

    for (size_t i = 0; i < out_traits_t::nrows(); ++i)
    {
        for (size_t j = 0; j < out_traits_t::ncols(); ++j)
            out_traits_t::value(out, i, j) = in_traits_t::value(in, i, j);
    }
    return out;
}

} // namespace spurt