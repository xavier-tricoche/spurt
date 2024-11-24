#pragma once

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <type_traits>

namespace spurt_internal {

// better_type:
// deciding which type should be
// used to perform a binary operation on a pair of types
template<typename T1, typename T2>
struct better_type {};

// always pick double
template<typename T, typename = 
typename std::enable_if<is_scalar<T>::value>::type >
struct better_type<T, double> {
    typedef double type;
};
template <typename T, typename =
typename std::enable_if<is_scalar<T>::value>::type >>
struct better_type<double, T>
{
    typedef double type;
};

// always pick floating point over integral
template <typename T1, typename T2, typename = 
typename std::enable_if < std::is_floating_point<T1>::value && 
                          std::is_integral<T2>::value
                        >::type >
struct better_type<T1, T2>
{
    typedef T1 type;
};
template <typename T1, typename T2, typename = typename std::enable_if<std::is_floating_point<T2>::value && std::is_integral<T1>::value>::type>
struct better_type<T1, T2>
{
    typedef T2 type;
};

// all integral types handled below:
// bool, 
// unsigned char, char, 
// unsigned short, short, 
// unsigned int, int, 
// unsigned long, long, 
// unsigned long long, long long

#define INTEGRAL_TYPES_MATCHES \
X(bool,               bool,               bool); \
X(bool,               unsigned char,      char); \
X(bool,               char,               char); \
X(bool,               unsigned short,     short); \
X(bool,               short,              short); \
X(bool,               unsigned int,       int); \
X(bool,               int,                int); \
X(bool,               unsigned long,      long); \
X(bool,               long,               long); \
X(bool,               unsigned long long, long long); \
X(bool,               long long,          long long); \
X(unsigned char,      bool,               char); \
X(unsigned char,      unsigned char,      char); \
X(unsigned char,      char,               char); \
X(unsigned char,      unsigned short,     short); \
X(unsigned char,      short,              short); \
X(unsigned char,      unsigned int,       int); \
X(unsigned char,      int,                int); \
X(unsigned char,      unsigned long,      long); \
X(unsigned char,      long,               long); \
X(unsigned char,      unsigned long long, long long); \
X(unsigned char,      long long,          long long); \
X(char,               bool,               char); \
X(char,               unsigned char,      char); \
X(char,               char,               char); \
X(char,               unsigned short,     short); \
X(char,               short,              short); \
X(char,               unsigned int,       int); \
X(char,               int,                int); \
X(char,               unsigned long,      long); \
X(char,               long,               long); \
X(char,               unsigned long long, long long); \
X(char,               long long,          long long); \
X(unsigned short,     bool,               short); \
X(unsigned short,     unsigned char,      char); \
X(unsigned short,     char,               char); \
X(unsigned short,     unsigned short,     short); \
X(unsigned short,     short,              short); \
X(unsigned short,     unsigned int,       int); \
X(unsigned short,     int,                int); \
X(unsigned short,     unsigned long,      long); \
X(unsigned short,     long,               long); \
X(unsigned short,     unsigned long long, long long); \
X(unsigned short,     long long,          long long); \
X(short,              bool,               short); \
X(short,              unsigned char,      short); \
X(short,              char,               short); \
X(short,              unsigned short,     short); \
X(short,              short,              short); \
X(short,              unsigned int,       int); \
X(short,              int,                int); \
X(short,              unsigned long,      long); \
X(short,              long,               long); \
X(short,              unsigned long long, long long); \
X(short,              long long,          long long); \
X(unsigned int,       bool,               int); \
X(unsigned int,       unsigned char,      int); \
X(unsigned int,       char,               int); \
X(unsigned int,       unsigned short,     int); \
X(unsigned int,       short,              unsigned int); \
X(unsigned int,       unsigned int,       int); \
X(unsigned int,       int,                int); \
X(unsigned int,       unsigned long,      long); \
X(unsigned int,       long,               long); \
X(unsigned int,       unsigned long long, long long); \
X(unsigned int,       long long,          long long); \
X(int,                bool,               int); \
X(int,                unsigned char,      int); \
X(int,                char,               int); \
X(int,                unsigned short,     int); \
X(int,                short,              int); \
X(int,                unsigned int,       int); \
X(int,                int,                int); \
X(int,                unsigned long,      long); \
X(int,                long,               long); \
X(int,                unsigned long long, long long); \
X(int,                long long,          long long); \
X(unsigned long,      bool,               long); \
X(unsigned long,      unsigned char,      long); \
X(unsigned long,      char,               long); \
X(unsigned long,      unsigned short,     long); \
X(unsigned long,      short,              long); \
X(unsigned long,      unsigned int,       long); \
X(unsigned long,      int,                long); \
X(unsigned long,      unsigned long,      long); \
X(unsigned long,      long,               long); \
X(unsigned long,      unsigned long long, long long); \
X(unsigned long,      long long,          long long); \
X(long,               bool,               long); \
X(long,               unsigned char,      long); \
X(long,               char,               long); \
X(long,               unsigned short,     long); \
X(long,               short,              long); \
X(long,               unsigned int,       long); \
X(long,               int,                long); \
X(long,               unsigned long,      long); \
X(long,               long,               long); \
X(long,               unsigned long long, long long); \
X(long,               long long,          long long); \
X(unsigned long long, bool,               long long); \
X(unsigned long long, unsigned char,      long long); \
X(unsigned long long, char,               long long); \
X(unsigned long long, unsigned short,     long long); \
X(unsigned long long, short,              long long); \
X(unsigned long long, unsigned int,       long long); \
X(unsigned long long, int,                long long); \
X(unsigned long long, unsigned long,      long long); \
X(unsigned long long, long,               long long); \
X(unsigned long long, unsigned long long, long long); \
X(unsigned long long, long long,          long long); \
X(long long,          bool,               long long); \
X(long long,          unsigned char,      long long); \
X(long long,          char,               long long); \
X(long long,          unsigned short,     long long); \
X(long long,          short,              long long); \
X(long long,          unsigned int,       long long); \
X(long long,          int,                long long); \
X(long long,          unsigned long,      long long); \
X(long long,          long,               long long); \
X(long long,          unsigned long long, long long); \
X(long long,          long long,          long long); \
#define X(type1,type2,type3) \
template <> \
struct better_type<type1,type2> \
{ \
    typedef type3 type; \
};    

}

namespace spurt {

    // fixed_array: 
    // a thin layer built on top of Eigen::Array
    template <typename T, size_t M, size_t N=1>
    class fixed_array: protected Eigen::Array<T, M, N>
    {
    public:
        typedef T scalar_type;
        typedef T value_type;
        typedef Eigen::Array<T, M, N> base_type;
        typedef Eigen::Matrix<T, M, N> matrix_type;
        typedef fixed_array<T, M, N> self_type;
        static constexpr size_t nrows = M;
        static constexpr size_t ncols = N;
        static constexpr bool is_vector = (N == 1);
        static constexpr bool is_matrix = (M > 1 && N > 1);
        template<typename T1>
        using other_array = Eigen::Array<T1, M, N>;
        template<typename T1>
        using other_matrix = Eigen::Matrix<T1, M, N>;
        template<typename T1>
        using other_farray = fixed_array<T1, M, N>;
        template<typename T1, size_t P>
        using rhs_farray = fixed_array<T1, N, P>;
        template<typename T1, size_t P>
        using lhs_farray = fixed_array<T1, P, M>;

        fixed_array(bool linalg=false) : m_linalg(linalg) { 
            base_type::setConstant(static_cast<scalar_type>(0));
        }

        template<typename T1>
        fixed_array(const other_farray<T1>& other) : base_type(other.as_eigen_array().template cast<scalar_type>()), m_linalg(other.m_linalg) {}

        template <typename T1>
        fixed_array(const other_array<T1> &a, bool linalg=false) : base_type(a.template cast<scalar_type>(), m_linalg(linalg)) {}
        template <typename T1>

        fixed_array(const other_matrix<T1> &m, bool linalg=false) : base_type(m.array().template cast<scalar_type>()), m_linalg(linalg) {}

        template <typename T1, typename =
                  typename std::enable_if<std::is_scalar<T1>::value>::type>
        fixed_array(const T1 &val) : m_linalg(false)
        {
            base_type::setConstant(static_cast<scalar_type>(val));
        }

        template<typename T1>
        fixed_array(const T1& val1, const T1& val2) : base_type(static_cast<scalar_type>(val1), static_cast<scalar_type>(val2)), m_linalg(false) {}

        template <typename T1>
        fixed_array(const T1 &val1, const T1 &val2, const T1 &val3) : base_type(val1, val2, val3), m_linalg(false) {}

        template <typename T1>
        fixed_array(const T1 &val1, const T1 &val2, const T1 &val3, const T1 &val4) : base_type(val1, val2, val3, val4), m_linalg(false) {}

        template <typename Op, typename Array>
        fixed_array(const Eigen::CwiseUnaryOp<Op, Array> &expr) : base_type(expr.eval().template cast<T>()), m_linalg(false) {}

        template<typename Op, typename RhsType, typename LhsType>
        fixed_array(const Eigen::CwiseBinaryOp<Op, RhsType, LhsType>& expr) : base_type(expr.eval().template cast<T>()), m_linalg(false) {}

        const scalar_type& operator[](size_t i) const {
            return base_type::operator[](i);
        }

        scalar_type& operator[](size_t i) {
            return base_type::operator[](i);
        }

        const scalar_type &operator()(size_t i, size_t j=0) const
        {
            return base_type::operator()(i,j);
        }

        scalar_type &operator ()(size_t i, size_t j=0)
        {
            return base_type::operator()(i,j);
        }

        size_t size() const
        {
            return base_type::size();
        }

        const base_type &as_eigen_array() const
        {
            return static_cast<const base_type &>(*this);
        }
        base_type &as_eigen_array()
        {
            return static_cast<base_type>(*this);
        }

        const matrix_type& linalg() const {
            return static_cast<const matrix_type&>(base_type::matrix());
        }

        matrix_type& linalg() {
            return static_cast<matrix_type&>(base_type::matrix());
        }

        // Comparison
        template<typename T1>
        other_farray<bool> operator<(const other_farray<T1>& other) const 
        {
            return other_farray<bool>(base_type::operator<(other.as_eigen_array().template cast<scalar_type>()));
        }
        template <typename T1>
        other_farray<bool> operator>(const other_farray<T1> &other) const
        {
            return other_farray<bool>(base_type::operator>(other.as_eigen_array().template cast<scalar_type>()));
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        other_farray<bool> operator<(const T1 &val) const
        {
            return other_farray<bool>(base_type::operator<(val));
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        other_farray<bool> operator>(const T1 &val) const
        {
            return other_farray<bool>(base_type::operator>(val));
        }
        template <typename T1>
        other_farray<bool> operator==(const other_farray<T1>& other) const 
        {
            return other_farray<bool>(base_type::operator==(other.as_eigen_array()));
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        other_farray<bool> operator==(const T1& val) const
        {
            return other_farray<bool>(base_type::operator==(val));
        }

        // Product
        template<typename T1>
        self_type& operator*=(const other_farray<T1>& other) {
            if (!m_linalg)
                base_type::operator*=(other.as_eigen_array().template cast<scalar_type>());
            else
                linalg() *= other.linalg();
            return (*this);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type &operator*=(const T1 &val)
        {
            base_type::operator*=(static_cast<scalar_type>(val));
            return (*this);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type operator*(const other_farray<T1> &other) const
        {
            if (!m_linalg)
                return self_type(base_type::operator*(other.as_eigen_array().template cast<T>()));
            else 
                return self_type(linalg() * other.linalg(), true);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type operator*(const T1 &val) const
        {
            return self_type(base_type::operator*(static_cast<scalar_type>(val)), m_linalg);
        }
        template <typename T1, size_t P>
        fixed_array<scalar_type, M, P> operator*(const rhs_farray<T1, P>& other) const
        {
            if (!m_linalg)
                std::cerr << "Warning: you are multiplying arrays of different shapes without <linalg> flag. Assuming implicit <true> value. \n";
            return fixed_array<T, M, P>(
                (base_type::matrix() * other.as_eigen_array().matrix()).eval(), m_linalg
            );
        }
        template <typename T1, size_t P>
        fixed_array<scalar_type, P, N> operator*(const lhs_farray<T1, P> &other) const
        {
            if (!m_linalg)
                std::cerr << "Warning: you are multiplying arrays of different shapes without <linalg> flag. Assuming implicit <true> value. \n";
            return fixed_array<T, P, M>(
                (other.as_eigen_array().matrix() * base_type::matrix()).eval(), m_linalg
            );
        }

        // Division
        template <typename T1>
        self_type &operator/=(const other_farray<T1> &other)
        {
            if (!m_linalg)
                base_type::operator/=(other.as_eigen_array().template cast<scalar_type>());
            else 
                linalg().operator*=(other.linalg());
            return (*this);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type &operator/=(const T1 &val)
        {
            base_type::operator/=(static_cast<scalar_type>(val));
            return (*this);
        }
        template<typename T1>
        self_type operator/(const other_farray<T1>& other) const
        {
            if (!m_linalg || other.m_linalg)
                return self_type(base_type::operator/(other.as_eigen_array().template cast<scalar_type>()), m_linalg);
            else 
                return self_type(linalg() * other.linalg().inverse(), true);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type operator/(const T1 &val) const
        {
            return self_type(base_type::operator/(static_cast<scalar_type>(val)), m_linalg);
        }

        // Addition
        template <typename T1>
        self_type &operator+=(const other_farray<T1> &other)
        {
            base_type::operator+=(other.as_eigen_array().template cast<scalar_type>());
            return (*this);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type &operator+=(const T1 &val)
        {
            base_type::operator+=(static_cast<scalar_type>(val));
            return (*this);
        }
        template<typename T1>
        self_type operator+(const other_farray<T1>& other) const
        {
            return self_type(base_type::operator+(other.as_eigen_array().template cast<scalar_type>()), m_linalg);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type operator+(const T1 &val) const
        {
            return self_type(base_type::operator+(static_cast<scalar_type>(val)));
        }

        // Subtraction
        template <typename T1>
        self_type &operator-=(const other_farray<T1> &other)
        {
            base_type::operator-=(other.as_eigen_array().template cast<scalar_type>());
            return (*this);
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type &operator-=(const T1 &val)
        {
            base_type::operator-=(static_cast<scalar_type>(val));
            return (*this);
        }
        template<typename T1>
        self_type operator-(const other_farray<T1>& other) const
        {
            return self_type(base_type::operator-(other.as_eigen_array().template cast<scalar_type>()));
        }
        template <typename T1, typename = typename std::enable_if<std::is_scalar<T1>::value>::type >
        self_type operator-(const T1 &val) const
        {
            return self_type(base_type::operator-(static_cast<scalar_type>(val)));
        }
    
        bool m_linalg;
    };

    template<typename Matrix>
    struct is_vector {
        static constexpr bool value = (Matrix::ColsAtCompileTime == 1);
    };

    template <typename T1, typename T2>
    fixed_array<T1, 3, 1> 
    cross(const fixed_array<T1, 3, 1> &a, const fixed_array<T2, 3, 1> &b)
    {
        return fixed_array<T1, 3, 1>(
                (a.as_eigen_array().matrix()).cross(b.as_eigen_array().template cast<T1>().matrix()));
    }

    template<typename T, size_t M, size_t N>
    T norm(const fixed_array<T, M, N>& a) {
        return a.as_eigen_array().matrix().norm();
    }

    template<typename T1, typename T2, size_t M, size_t N>
    T1
    inner(const fixed_array<T1, M, N>& a, const fixed_array<T2, M, N>& b) {
        return 
            a.as_eigen_array().matrix().dot(b.as_eigen_array().template cast<T1>().matrix());
    }

    template<typename T, size_t M>
    T det(const fixed_array<T, M, M>& a)
    {
        return a.as_eigen_array().matrix().determinant();
    }

    template <typename T1, typename T2, size_t M, size_t N,
              typename =
                  typename std::enable_if<std::is_scalar<T2>::value>::type>
    fixed_array<T1, M, N> 
    operator+(const T2 &val, const fixed_array<T1, M, N> &a)
    {
        return a + val;
    }

    template <typename T1, typename T2, size_t M, size_t N,
              typename =
                  typename std::enable_if<std::is_scalar<T2>::value>::type>
    fixed_array<T1, M, N> 
    operator-(const T2 &val, const fixed_array<T1, M, N> &a)
    {
        return a - val;
    }

    template <typename T1, typename T2, size_t M, size_t N, 
              typename = 
              typename std::enable_if<std::is_scalar<T2>::value>::type >
    fixed_array<T1, M, N> 
    operator*(const T2 &val, const fixed_array<T1, M, N> &a)
    {
        return a * val;
    }

    template <size_t M, size_t N>
    bool all(const fixed_array<bool, M, N> &a)
    {
        return a.as_eigen_array().all();
    }

    template<size_t M, size_t N>
    bool any(const fixed_array<bool, M, N>& a) {
        return a.as_eigen_array().any();
    }

    template<typename T, size_t M, size_t N>
    T min(const fixed_array<T, M, N>& a) {
        return a.as_eigen_array().minCoeff();
    }

    template <typename T, size_t M, size_t N>
    T max(const fixed_array<T, M, N> &a)
    {
        return a.as_eigen_array().maxCoeff();
    }

    template<typename T, size_t M, size_t N>
    fixed_array<T, M, N> abs(const fixed_array<T, M, N>& a) 
    {
        return fixed_array<T, M, N>(a.as_eigen_array().abs());
    }

    template <typename T, size_t M, size_t N>
    fixed_array<T, M, N> square(const fixed_array<T, M, N> &a)
    {
        return fixed_array<T, M, N>(a.as_eigen_array().square());
    }

    template<typename T, size_t N>
    using fixed_vector = fixed_array<T, N, 1>;

    template<typename T, size_t N>
    using fixed_matrix = fixed_array<T, N, N>;

    template<typename T, size_t M, size_t N>
    std::ostream& operator<<(std::ostream& os, const fixed_array<T, M, N>& a)
    {
        typedef fixed_array<T, M, N> a_type;
        typedef typename a_type::base_type b_type; 
        const b_type& arr = a.as_eigen_array();
        const size_t ncols = arr.cols();
        const size_t nrows = arr.rows();
        const size_t sz = arr.size();
        if (ncols == 1 || nrows==1)
        {
            os << "[";
            for (size_t i=0; i<sz-1; ++i) {
                os << arr.reshaped()[i] << ",";
            }
            os << arr.reshaped()[sz-1];
            os << "]";
        }
        else {
            os << "ncols = " << ncols << ", nrows = " << nrows << ' \n';
            os << " \n" << arr;
        }
        return os;
    }

    typedef fixed_vector<double, 1> vec1;
    typedef fixed_vector<double, 2> vec2;
    typedef fixed_vector<double, 3> vec3;
    typedef fixed_vector<double, 4> vec4;
    typedef fixed_vector<double, 5> vec5;
    typedef fixed_vector<double, 6> vec6;

    typedef fixed_vector<double, 1> dvec1;
    typedef fixed_vector<double, 2> dvec2;
    typedef fixed_vector<double, 3> dvec3;
    typedef fixed_vector<double, 4> dvec4;
    typedef fixed_vector<double, 5> dvec5;
    typedef fixed_vector<double, 6> dvec6;

    typedef fixed_vector<float, 1> fvec1;
    typedef fixed_vector<float, 2> fvec2;
    typedef fixed_vector<float, 3> fvec3;
    typedef fixed_vector<float, 4> fvec4;
    typedef fixed_vector<float, 5> fvec5;
    typedef fixed_vector<float, 6> fvec6;

    typedef fixed_vector<int, 1> ivec1;
    typedef fixed_vector<int, 2> ivec2;
    typedef fixed_vector<int, 3> ivec3;
    typedef fixed_vector<int, 4> ivec4;
    typedef fixed_vector<int, 5> ivec5;
    typedef fixed_vector<int, 6> ivec6;

    typedef fixed_vector<long int, 1> lvec1;
    typedef fixed_vector<long int, 2> lvec2;
    typedef fixed_vector<long int, 3> lvec3;
    typedef fixed_vector<long int, 4> lvec4;
    typedef fixed_vector<long int, 5> lvec5;
    typedef fixed_vector<long int, 6> lvec6;

    typedef fixed_vector<unsigned int, 1> uvec1;
    typedef fixed_vector<unsigned int, 2> uvec2;
    typedef fixed_vector<unsigned int, 3> uvec3;
    typedef fixed_vector<unsigned int, 4> uvec4;
    typedef fixed_vector<unsigned int, 5> uvec5;
    typedef fixed_vector<unsigned int, 6> uvec6;

    typedef fixed_vector<size_t, 1> svec1;
    typedef fixed_vector<size_t, 2> svec2;
    typedef fixed_vector<size_t, 3> svec3;
    typedef fixed_vector<size_t, 4> svec4;
    typedef fixed_vector<size_t, 5> svec5;
    typedef fixed_vector<size_t, 6> svec6;

    typedef fixed_matrix<double, 2> mat2;
    typedef fixed_matrix<double, 3> mat3;
    typedef fixed_matrix<double, 4> mat4;
    typedef fixed_matrix<double, 5> mat5;
    typedef fixed_matrix<double, 6> mat6;

    typedef fixed_matrix<double, 2> dmat2;
    typedef fixed_matrix<double, 3> dmat3;
    typedef fixed_matrix<double, 4> dmat4;
    typedef fixed_matrix<double, 5> dmat5;
    typedef fixed_matrix<double, 6> dmat6;

    typedef fixed_matrix<float, 2> fmat2;
    typedef fixed_matrix<float, 3> fmat3;
    typedef fixed_matrix<float, 4> fmat4;
    typedef fixed_matrix<float, 5> fmat5;
    typedef fixed_matrix<float, 6> fmat6;

    typedef fixed_matrix<int, 2> imat2;
    typedef fixed_matrix<int, 3> imat3;
    typedef fixed_matrix<int, 4> imat4;
    typedef fixed_matrix<int, 5> imat5;
    typedef fixed_matrix<int, 6> imat6;

    struct lexicographical_order
    {
        template <typename T, size_t M, size_t N>
        bool operator()(const fixed_array<T, M, N> &v0,
                        const fixed_array<T, M, N> &v1) const
        {
            auto a = v0.as_eigen_array().reshaped();
            auto b = v1.as_eigen_array().reshaped();
            for (auto i=0; i<a.size(); ++i) 
            {
                if (a[i] < b[i])
                    return true;
                else if (a[i] > b[i])
                    return false;
            }
            return false;
        }

        template <typename T, size_t M, size_t N>
        static bool equal(const fixed_array<T, M, N> &v0, 
                          const fixed_array<T, M, N> &v1)
        {
            auto a = v0.as_eigen_array().reshaped();
            auto b = v1.as_eigen_array().reshaped();

            for (auto i = 0; i < a.size(); ++i)
            {
                if (a[i] != b[i])
                    return false;
            }
            return true;
        }

        template <typename T, size_t M, size_t N>
        static bool equal_zero(const fixed_array<T, M, N> &v0)
        {
            auto a = v0.as_eigen_array().reshaped();
            for (auto i = 0; i < a.size(); ++i)
            {
                if (a[i] != 0)
                    return false;
            }
            return true;
        }
    };

    struct eps_lexicographical_order
    {
        eps_lexicographical_order(double eps) : m_eps(eps) {}

        template <typename T, size_t M, size_t N>
        bool operator()(const fixed_array<T, M, N> &v0,
                        const fixed_array<T, M, N> &v1) const
        {
            auto a = v0.as_eigen_array().reshaped();
            auto b = v1.as_eigen_array().reshaped();
            for (auto i = 0; i < a.size(); ++i)
            {
                if (a[i] < b[i] - m_eps)
                    return true;
                else if (a[i] > b[i] + m_eps)
                    return false;
            }
            return false;
        }

        template <typename T, size_t M, size_t N>
        bool equal(const fixed_array<T, M, N> &v0, 
                   const fixed_array<T, M, N> &v1) const
        {
            auto a = v0.as_eigen_array().reshaped();
            auto b = v1.as_eigen_array().reshaped();
            for (auto i = 0; i < a.size(); ++i)
            {
                if (fabs(a[i] - b[i]) > m_eps)
                    return false;
            }
            return true;
        }

        template <typename T, size_t M, size_t N>
        bool equal_zero(const fixed_array<T, M, N> &v0) const
        {
            auto a = v0.as_eigen_array().reshaped();
            for (auto i = 0; i < a.size(); ++i)
            {
                if (fabs(a[i]) > m_eps)
                    return false;
            }
            return true;
        }

        double m_eps;
    };

} // namespace spurt