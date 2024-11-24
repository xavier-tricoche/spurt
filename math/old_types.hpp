#pragma once

#include <Eigen/Eigen>
#include <misc/meta_utils.hpp>
#include <math/bounding_box.hpp>

namespace spurt {


    template <typename T, int N>
    using fixed_vector = Eigen::Vector<T, N>;

    template <typename T, int N>
    using fixed_matrix = Eigen::Matrix<T, N, N>;

    typedef Eigen::Vector<double, 1> vec1;
    typedef Eigen::Vector<double, 2> vec2;
    typedef Eigen::Vector<double, 3> vec3;
    typedef Eigen::Vector<double, 4> vec4;
    typedef Eigen::Vector<double, 5> vec5;
    typedef Eigen::Vector<double, 6> vec6;

    typedef Eigen::Vector<double, 1> dvec1;
    typedef Eigen::Vector<double, 2> dvec2;
    typedef Eigen::Vector<double, 3> dvec3;
    typedef Eigen::Vector<double, 4> dvec4;
    typedef Eigen::Vector<double, 5> dvec5;
    typedef Eigen::Vector<double, 6> dvec6;

    typedef Eigen::Vector<float, 1> fvec1;
    typedef Eigen::Vector<float, 2> fvec2;
    typedef Eigen::Vector<float, 3> fvec3;
    typedef Eigen::Vector<float, 4> fvec4;
    typedef Eigen::Vector<float, 5> fvec5;
    typedef Eigen::Vector<float, 6> fvec6;

    typedef Eigen::Vector<int, 1> ivec1;
    typedef Eigen::Vector<int, 2> ivec2;
    typedef Eigen::Vector<int, 3> ivec3;
    typedef Eigen::Vector<int, 4> ivec4;
    typedef Eigen::Vector<int, 5> ivec5;
    typedef Eigen::Vector<int, 6> ivec6;

    typedef Eigen::Vector<long int, 1> lvec1;
    typedef Eigen::Vector<long int, 2> lvec2;
    typedef Eigen::Vector<long int, 3> lvec3;
    typedef Eigen::Vector<long int, 4> lvec4;
    typedef Eigen::Vector<long int, 5> lvec5;
    typedef Eigen::Vector<long int, 6> lvec6;

    typedef Eigen::Vector<unsigned int, 1> uvec1;
    typedef Eigen::Vector<unsigned int, 2> uvec2;
    typedef Eigen::Vector<unsigned int, 3> uvec3;
    typedef Eigen::Vector<unsigned int, 4> uvec4;
    typedef Eigen::Vector<unsigned int, 5> uvec5;
    typedef Eigen::Vector<unsigned int, 6> uvec6;

    typedef Eigen::Vector<size_t, 1> svec1;
    typedef Eigen::Vector<size_t, 2> svec2;
    typedef Eigen::Vector<size_t, 3> svec3;
    typedef Eigen::Vector<size_t, 4> svec4;
    typedef Eigen::Vector<size_t, 5> svec5;
    typedef Eigen::Vector<size_t, 6> svec6;

    typedef Eigen::Matrix<double, 2, 2> mat2;
    typedef Eigen::Matrix<double, 3, 3> mat3;
    typedef Eigen::Matrix<double, 4, 4> mat4;
    typedef Eigen::Matrix<double, 5, 5> mat5;
    typedef Eigen::Matrix<double, 6, 6> mat6;

    typedef Eigen::Matrix<double, 2, 2> dmat2;
    typedef Eigen::Matrix<double, 3, 3> dmat3;
    typedef Eigen::Matrix<double, 4, 4> dmat4;
    typedef Eigen::Matrix<double, 5, 5> dmat5;
    typedef Eigen::Matrix<double, 6, 6> dmat6;

    typedef Eigen::Matrix<float, 2, 2> fmat2;
    typedef Eigen::Matrix<float, 3, 3> fmat3;
    typedef Eigen::Matrix<float, 4, 4> fmat4;
    typedef Eigen::Matrix<float, 5, 5> fmat5;
    typedef Eigen::Matrix<float, 6, 6> fmat6;

    typedef Eigen::Matrix<int, 2, 2> imat2;
    typedef Eigen::Matrix<int, 3, 3> imat3;
    typedef Eigen::Matrix<int, 4, 4> imat4;
    typedef Eigen::Matrix<int, 5, 5> imat5;
    typedef Eigen::Matrix<int, 6, 6> imat6;
    
    template<typename Derived>
    typename Eigen::MatrixBase<Derived>::Scalar 
    norm(const Eigen::MatrixBase<Derived>& m) {
        return m.norm();
    }

    template<typename Derived, typename OtherDerived>
    typename Eigen::MatrixBase<Derived>::Scalar
    inner(const Eigen::MatrixBase<Derived> &m0, 
          const Eigen::MatrixBase<OtherDerived> &m1)
    {
        return m0.dot(m1);
    }

    template <typename T>
    Eigen::Matrix<T, 3, 1>
    cross(const Eigen::Matrix<T, 3, 1> &m0,
          const Eigen::Matrix<T, 3, 1> &m1)
    {
        return m0.cross(m1);
    }

    template <typename Matrix>
    struct flatmat {
        flatmat(const Matrix& mat) : m_mat(mat) {} 
        const Matrix& m_mat;
    };

    template<typename Vector>
    std::string print_vector(const Vector& vec) {
        auto s = vec.size();
        std::ostringstream oss;
        oss << "[";
        for (auto i = 0; i < s; ++i)
        {
            oss << vec[i];
            if (i < s - 1)
            {
                oss << ", ";
            }
        }
        oss << "]";
        return oss.str();
    }

    template<typename Matrix>
    std::ostream& operator<<(std::ostream& os, const flatmat<Matrix>& fmat) {
        auto nrows = fmat.m_mat.rows();
        auto ncols = fmat.m_mat.cols();
        if (ncols == 1 || nrows == 1) {
            os << print_vector(fmat.m_mat);
        }
        else {
            os << "[ ";
            for (auto row=0 ; row<nrows; ++row) {
                os << print_vector(fmat.m_mat.row(row));
                if (row<nrows-1) os << ", ";
            }
            os << " ]";
        }
        return os;
    }

    template<typename Matrix>
    std::string to_str(const flatmat<Matrix>& m) {
        std::ostringstream oss;
        oss << m;
        return oss.str();
    }

    struct lexicographical_order
    {
        template<typename Vector>
        bool operator()(const Vector& v0, 
                        const Vector& v1) const {
            typedef data_traits<Vector> traits_t;
            for (size_t i=0; i<traits_t::size(); ++i) {
                if (traits_t::value(v0, i) < traits_t::value(v1, i)) 
                    return true; 
                else if (traits_t::value(v0, i) > traits_t::value(v1, i)) 
                    return false;
            }
            return false;
        }

        template<typename Vector>
        static bool equal(const Vector& v0, const Vector& v1) {
            typedef data_traits<Vector> traits_t;
            for (size_t i = 0; i < traits_t::size(); ++i) {
                if (traits_t::value(v0, i) != traits_t::value(v1, i)) return false;
            }
            return true;
        }

        template<typename Vector>
        static bool equal_zero(const Vector& v0) {
            typedef data_traits<Vector> traits_t;
            for (size_t i = 0; i < traits_t::size(); ++i)
            {
                if (traits_t::value(v0, i) != 0) return false;
            }
            return true;
        }
    };

    struct eps_lexicographical_order
    {
        eps_lexicographical_order(double eps) : m_eps(eps) {}
    
        template<typename Vector>
        bool operator()(const Vector& v0, 
                        const Vector& v1) const {
            typedef data_traits<Vector> traits_t;
            for (size_t i=0; i<traits_t::size(); ++i) {
                if (traits_t::value(v0,i) < traits_t::value(v1, i) - m_eps) 
                    return true;
                else if (traits_t::value(v0, i) > traits_t::value(v1,i) + m_eps)
                    return false;
            }
            return false;
        }

        template<typename Vector>
        bool equal(const Vector& v0, const Vector& v1) const {
            typedef data_traits<Vector> traits_t;
            for (size_t i=0; i<traits_t::size(); ++i) {
                if (fabs(traits_t::value(v0, i)-traits_t::value(v1, i)) > m_eps) 
                    return false;
            }
            return true;
        }

        template<typename Vector>
        bool equal_zero(const Vector& v0) const {
            typedef data_traits<Vector> traits_t;
            for (size_t i = 0; i < traits_t::size(); ++i)
            {
                if (fabs(traits_t::value(v0, i)) > m_eps)
                    return false;
            }
            return true;
        }

        double m_eps;
    };

    template<typename Value>
    std::ostream& operator<<(std::ostream& os, const bounding_box<Value>& bbox) {
        os << "[" << bbox.min() << " - " << bbox.max() << "]";
        return os;
    }

    typedef bounding_box<vec2> bbox2;
    typedef bounding_box<vec3> bbox3;
    typedef bounding_box<vec4> bbox4;
    typedef bounding_box<vec5> bbox5;
    typedef bounding_box<vec6> bbox6;

    template<typename Matrix>
    struct matrix : public Matrix {
        typedef Matrix base_type;
        typedef data_traits<Matrix> traits_type;
        typedef typename traits_type::value_type value_type;
        static const size_t dimension = traits_type::size(); 
        typedef matrix<Matrix> self_type;

        matrix() : base_type() {}
        matrix(const Matrix& m) : base_type(m) {}

        const value_type& operator()(size_t i=0, size_t j=0) const {
            return traits_type::value((*this), i, j);
        }

        value_type& operator()(size_t i=0, size_t j=0) {
            return traits_type::value((*this), i, j);
        }

        self_type dotprod(const self_type& other) {
            return traits_type::dotprod((*this), other);
        }
    };

    template<typename Matrix>
    bool any(const Matrix& m) {
        typedef data_traits<Matrix> traits;
        for (int i=0; i<traits::nrows(); ++i) {
            for (int j=0; j<traits::ncols(); ++j)
                if (traits::value(m, i, j)) return true;
        }
        return false;
    }

    template<typename Matrix>
    bool all(const Matrix& m) {
        typedef data_traits<Matrix> traits;
        for (int i=0; i<traits::nrows(); ++i) {
            for (int j=0; j<traits::ncols(); ++j)
                if (!traits::value(m, i, j)) return false;
        }
        return true;
    }

} // namespace spurt