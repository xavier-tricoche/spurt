#include <bitset>
#include <utility>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <misc/strings.hpp>
#include <misc/meta_utils.hpp>

namespace {
template< typename Float_, typename = typename std::enable_if< std::is_floating_point<Float_>::value, void >::type >
struct eps_compare {
    typedef Float_ float_t;
    static constexpr float_t eps = static_cast<float_t>(10)*std::numeric_limits<float_t>::epsilon();

    static bool equal(float_t u, float_t v) {
        // std::cout << std::setprecision(16) << eps << '\n';
        return ( std::abs(u-v)/std::abs(u) <= eps ) ||
               ( std::abs(u-v)/std::abs(v) <= eps );
    }

    static bool equal_zero(float_t u) {
        return std::abs(u) <= eps;
    }

    template< typename Compare_ = std::less< float_t > >
    static bool strict_compare(float_t u, float_t v, Compare_ comp ) {
        if ( equal(u, v) ) return false;
        return comp(u, v);
    }
};
}

namespace spurt {

template<typename Size_>
inline int fix_index(Size_ a, Size_ s) {
    if (a < 0) return 0;
    else if (a >= s) return s-1;
    else return a;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, typename Kernel_>
inline typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
convolve1D(const vec_type &u, const coord_typ& c, const coord_type& dorder) const
{
    typedef traits = spurt::data_traits<vec_type>;
    value_type r = traits::zero();
    scalar_type allw = 0;
    auto res = this->grid().m_res;
    const size_t width = kernel_type::size;
    for (int i = 1-width; i <= width; ++i) {
        int ii = fix_index(c[0] + i, res[0]);
        scalar_type w = m_kernel(u[0] - i, dorder[0]);
        allw += w;
        r +=  w * this->m_data[ii];
    }
    return r/allw;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, typename Kernel_>
inline typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
convolve2D(const vec_type &u, const coord_type& c, const coord_type& dorder) const
{
    typedef traits = spurt::data_traits<vec_type>;
    value_type r = traits::zero();
    scalar_type allw = 0;
    auto grid = this->grid();
    auto res = grid.m_res;
    const size_t width = kernel_type::size;
    for (int j = 1-width; j <= width; ++j)
    {
        int jj = fix_index(j+c[1], res[1]);
        scalar_type w = m_kernel(u[1] - j, dorder[1]);
        for (int i = 1-width; i <= width; ++i)
        {
            int ii = fix_index(i+c[0], res[0]);
            w *= m_kernel(u[0] - i, dorder[0]);
            allw += w;
            r += w * this->m_data[grid.index(ii, jj)];
        }
    }
    return r/allw;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, typename Kernel_>
inline typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
convolve3D(const vec_type &u, const coord_type& c, const coord_type& dorder) const
{
    typedef traits = spurt::data_traits<vec_type>;
    value_type r = traits::zero();
    scalar_type allw = 0;
    auto grid = this->grid();
    auto res = grid.m_res;
    const size_t width = kernel_type::size;
    for (int k = 1-width; k <= width; ++k)
    {
        int kk = fix_index(k+c[2], res[2]);
        scalar_type w = m_kernel(u[2] - k, dorder[2]);
        for (int j = 1-width; j <= width; ++j)
        {
            int jj = fix_index(j+c[1], res[1]);
            w *= m_kernel(u[1] - j, dorder[1]);
            for (int i = 1-width; i <= width; ++i)
            {
                int ii = fix_index(i+c[0], res[0]);
                w *= m_kernel(u[0] - i, dorder[0]);
                allw += w;
                r += w * this->m_data[grid.index(ii, jj, kk)];
            }
        }
    }
    return r / allw;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, typename Kernel_>
inline typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
convolve4D(const vec_type &u, const coord_type& c, const coord_type& dorder) const
{
    typedef traits = spurt::data_traits<vec_type>;
    value_type r = traits::zero();
    scalar_type allw = 0;
    auto grid = this->grid();
    auto res = grid.m_res;
    const size_t width = kernel_type::size;
    for (int l = 1-width; l <= width; ++l)
    {
        int ll = fix_index(l+c[3], res[3]);
        scalar_type w = m_kernel(u[3] - l, dorder[3]);
        for (int k = 1-width; k <= width; ++k)
        {
            int kk = fix_index(k+c[2], res[2]);
            w *= m_kernel(u[2] - k, dorder[2]);
            for (int j = 1-width; j <= width; ++j)
            {
                int jj = fix_index(j+c[1], res[1]);
                w *= m_kernel(u[1] - j, dorder[1]);
                for (int i = 1-width; i <= width; ++i)
                {
                    int ii = fix_index(i+c[0], res[0]);
                    w *= m_kernel(u[0] - i, dorder[0]);
                    allw += w;
                    r += w * this->m_data[grid.index(ii, jj, kk)];
                }
            }
        }
    }
    return r / allw;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
convolve(const vec_type &u, coord_type c, const coord_type &dorder, 
         size_t n) const
{
    switch (n)
    {
    case 0:
        throw std::runtime_error("invalid dimension in differential_image::convolve()");
    case 1:
        return convolve1D(u, c, dorder);
    case 2:
        return convolve2D(u, c, dorder);
    case 3:
        return convolve3D(u, c, dorder);
    case 4:
        return convolve4D(u, c, dorder);
    default:
        throw std::runtime_error("invalid dimension in differential_image::convolve()");
    }
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
interpolate(const vec_type &x, const coord_type &c) const
{
    return convolve(x, c);
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_,
          typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
value(const point_type &p) const
{
    std::pair<coord_type, vec_type> r = this->m_grid.locate(p);
    return interpolate(r.second, r.first);
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, 
          typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::value_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
value_in_voxel(const coord_type &vid, const point_type &p) const
{
    return interpolate(p, vid);
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_,
          typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::deriv_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
derivative_in_voxel(const coord_type &vid, const point_type &p) const
{
    deriv_type df;
    typedef traits = spurt::data_traits<value_type>;
    typedef ctraits = spurt::data_traits<coord_type>;
    coord_type order;
    ctraits::assign(order, 0);
    for (d = 0; d < dim; ++d)
    {
        traits::assign(df[d], 0);
        order[d] = 1;
        df[d] = convolve(p, vid, order);
        df[d] /= this->grid().spacing()[d];
        order[d] = 0;
    }
    return df;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, 
          typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::second_deriv_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
second_derivative_in_voxel(const coord_type &vid, const point_type &p) const
{
    second_deriv_type d2f;
    typedef traits = spurt::data_traits<value_type>;
    typedef ctraits = spurt::data_traits<coord_type>;
    coord_type order;
    ctraits::assign(order, 0);
    int offset=0;
    for (int d = 0; d < dim; ++d)
    {
        traits::assign(d2f[offset], 0);
        order[d] = 2; // diagonal
        d2f[offset] = convolve(p, vid, order);
        d2f[offset] /= this->grid().spacing()[d]*this->grid().spacing()[d];
        order[d] = 1;
        for (int dd=d+1; dd<dim; ++dd) {
            order[dd] = 1;
            offset += 1;
            d2f[offset] = convolve(p, vid, order);
            d2f[offset] /= this->grid().spacing()[d] * 
                           this->grid().spacing()[dd];
            order[dd] = 0;
        }
        order[d] = 0;
    }
    return d2f;
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, 
          typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::deriv_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel +>::
derivative(const point_type &p) const
{
    std::pair<coord_type, vec_type> r = this->m_grid.locate(p);
    return derivative_in_voxel(r.first, r.second);
}

template <typename Value_, size_t Dim_, typename Scalar_, typename Size_, 
          typename Kernel_>
typename differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::second_deriv_type
differential_image<Value_, Dim_, Scalar_, Size_, Kernel_>::
second_derivative(const point_type &p) const
{
    std::pair<coord_type, vec_type> r = this->m_grid.locate(p);
    return second_derivative_in_voxel(r.first, r.second);
}

} // namespace spurt
