#pragma once

#include <array>
#include <cmath>
#include <mutex>
#include <thread>

#include <math/types.hpp>
#include <data/raster_data.hpp>
#include <misc/progress.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
namespace spurt
{
namespace kernels
{
    // bleeding boundary handling
    template<typename Size_>
    inline int fix_index(Size_ a, Size_ s) {
        if (a < 0) return 0;
        else if (a >= s) return s-1;
        else return a;
    }

    struct Linear
    {
        static constexpr size_t size = 1;
        Linear() {}

        double operator()(double t, size_t dorder = 0) const
        {
            double u = (t < 0) ? -t : +t;
            if (u >= 1)
                return 0;
            else if (dorder == 0)
                return 1. - u;
            else if (dorder == 1)
                return (t < 0) ? 1. : -1.;
            else
                return 0.;
        }
    };

    struct MitchellNetravaliBC
    {
        static constexpr size_t size = 2;
        MitchellNetravaliBC(double B = 0, double C = 0.5) : m_B(B), m_C(C) {}

        double operator()(double t, size_t dorder = 0) const
        {
            double u = (t < 0) ? -t : t;
            if (u >= 2)
                return 0.;

            if (dorder == 0)
            {
                double usq = u * u;
                double ucu = usq * u;
                if (u < 1)
                {
                    return (2. - 1.5 * m_B - m_C) * ucu +
                           (-3 + 2 * m_B + m_C) * usq +
                           (1. - m_B / 3.);
                }
                else
                {
                    return (-m_B / 6. - m_C) * ucu +
                           (m_B + 5. * m_C) * usq +
                           (-2. * m_B - 8. * m_C) * u +
                           (4. * m_B / 3. + 4 * m_C);
                }
            }
            else if (dorder == 1)
            {
                // d|t|^3/dt: 3t^2 / -3*t^2
                // d|t|^2/dt: 2t
                // d|t|/dt: 1 / -1
                double dusq = 2 * t;
                double ducu = (t < 0) ? -3. * t * t : 3 * t * t;
                double du = (t < 0) ? -1. : 1.;
                if (u < 1.)
                    return (2. - 1.5 * m_B - m_C) * ducu +
                           (-3 + 2 * m_B + m_C) * dusq;
                else
                    return (-m_B / 6. - m_C) * ducu +
                           (m_B + 5. * m_C) * dusq +
                           (-2. * m_B - 8. * m_C) * du;
            }
            else if (dorder == 2)
            {
                // d2|t|^3/dt2: 6t / -6t
                // d2|t|^2/dt2: 2
                // d2|t|/dt2: 0
                double ddusq = 2;
                double dducu = (t < 0) ? -6. * t : 6. * t;
                if (u < 1.)
                    return (2. - 1.5 * m_B - m_C) * dducu +
                           (-3 + 2 * m_B + m_C) * ddusq;
                else
                    return (-m_B / 6. - m_C) * dducu +
                           (m_B + 5. * m_C) * ddusq;
            }
            else if (dorder == 3)
            {
                // d3|t|^3/dt3: 6 / -6
                double ddducu = (t < 0) ? -6. : 6.;
                if (u < 1.)
                    return (2. - 1.5 * m_B - m_C) * ddducu;
                else
                    return (-m_B / 6. - m_C) * ddducu;
            }
            else
                return 0.;
        }

        double m_B, m_C;
    };

    struct Gaussian {
        static constexpr double invsqrt2pi = 0.3989422804014327;
        size_t size;

        double _gauss(double t) const {
            return std::exp(-t * t / (2. * m_sig2));
        }

        Gaussian(double stddev=2, double cutoff=2)
            : m_sig(stddev), m_cutoff(cutoff), m_sig2(stddev * stddev)
        {
            m_sig4 = m_sig2 * m_sig2;
            m_sig6 = m_sig2 * m_sig4;
            size = std::ceil(m_sig * m_cutoff);
            m_norm = invsqrt2pi/m_sig;
        }

        double operator()(double t, size_t dorder=0) const {
            double u = (t < 0) ? -t : t;
            if (u >= size) return 0.;

            if (dorder == 0) {
                return m_norm * _gauss(t);
            }
            else if (dorder == 1) {
                return -m_norm * t / m_sig2 * _gauss(t);
            }
            else if (dorder == 2) {
                return m_norm * (t * t - m_sig2) / m_sig2 * _gauss(t);
            }
            else if (dorder == 3) {
                return -m_norm * t * (t*t - 3*m_sig2) / m_sig6 * _gauss(t);
            }
            else {
                throw std::runtime_error("invalid derivative order");
            }
        }

        double m_sig, m_sig2, m_sig4, m_sig6, m_cutoff, m_norm;
    };

} // namespace kernels

namespace differential
{
    template<typename T, size_t SpatialDim, typename Enable = void>
    struct derivatives {};
    
    template<typename T, size_t SpatialDim>
    struct derivatives<T, SpatialDim,
                       typename std::enable_if<std::is_scalar<T>::value>::type>
    {
        typedef T value_type;
        typedef T scalar_type;
        typedef small_vector<T, SpatialDim> first_derivative_type;
        typedef small_matrix<T, SpatialDim> second_derivative_type;
        
        static void set_dfdxi(int i, const value_type& v, first_derivative_type& df) {
            df[i] = v;
        }
        
        static void set_d2fdxij(int i, int j, const value_type& v, second_derivative_type& d2f) {
            d2f(i, j) = v;
        }
    };
   
    template<typename Storage, size_t SpatialDim>
    struct derivatives<small_vector_interface<Storage>, SpatialDim, void>
    {
        typedef typename Storage::value_type scalar_type;
        typedef small_vector_interface<Storage> value_type;
        typedef small_matrix<scalar_type, SpatialDim, SpatialDim> first_derivative_type;
        typedef small_tensor<scalar_type, SpatialDim, SpatialDim, SpatialDim> second_derivative_type;

        
        static void set_dfdxi(int i, const value_type& v, first_derivative_type& df) {
            df.column(i) = v;
        }
        // static const value_type& dfdxi(int i, const first_derivative_type& df) {
        //     return df.column(i);
        // }
        
        static void set_d2fdxij(int i, int j, const value_type& v, second_derivative_type& d2f) {
            d2f.column(i, j) = v;
        }
        // static const value_type& d2fdxij(int i, int j, const second_derivative_type& d2f) {
        //     return d2f.column(i, j);
        // }
        
    };
   
    template<typename Storage, size_t SpatialDim>
    struct derivatives<small_matrix_interface<Storage, SpatialDim, SpatialDim>, 
                       SpatialDim, void>
    {
        typedef typename Storage::value_type scalar_type;
        typedef small_matrix_interface<Storage, SpatialDim, SpatialDim> value_type;
        typedef small_tensor<scalar_type, SpatialDim, SpatialDim, SpatialDim> first_derivative_type;
        typedef small_tensor<scalar_type, SpatialDim, SpatialDim, SpatialDim*SpatialDim> second_derivative_type;
        
        static void set_dfdxi(int i, const value_type& v, 
                              first_derivative_type& df) {
            df.layer(i) = v;
        }
        
        static void set_d2fdxij(int i, int j, const value_type& v, 
                                second_derivative_type& d2f) {
            // there are N * N layers
            // J(i,j) = dfi/dxj
            // d2f/dxidxj = layers[i+j*N] 
            
            d2f.layer(i+j*SpatialDim) = v;
        }
        
    };
}

template <typename Size_, typename Scalar_, size_t Dim, typename Value_,
          typename Kernel_=kernels::Linear, typename Coord_ = small_vector<Size_, Dim>,
          typename Pos_ = small_vector<Scalar_, Dim> >
class image : public raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_>
{
public:
    typedef raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_> base_type;
    typedef typename base_type::size_type size_type;
    static const size_type dimension = base_type::dimension;

    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::pos_type pos_type;
    typedef typename base_type::coord_type coord_type;
    typedef typename base_type::bounds_type bounds_type;
    typedef typename base_type::grid_type grid_type;
    typedef typename base_type::value_type value_type;
    typedef data_traits<value_type> value_traits;

    typedef differential::derivatives<value_type, Dim> deriv_info_type;
    typedef typename deriv_info_type::first_derivative_type first_derivative_type;
    typedef typename deriv_info_type::second_derivative_type second_derivative_type;
    typedef typename deriv_info_type::value_type derivative_element;
    
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef image<Size_, Scalar_, Dim, Value_, Kernel_, Coord_, Pos_> self_type;

    template<typename OtherKernel_>
    using matching_image_type = image<Size_, Scalar_, Dim, Value_, OtherKernel_, Coord_, Pos_>;
    typedef Kernel_ kernel_type;

private:
    value_type interpolate(const pos_type &x, const coord_type &c) const {
        return convolve(x, c);
    }

    value_type convolve1D(const pos_type &u, const coord_type &c, 
                          const coord_type &dorder = coord_type(0)) 
                          const {  
        value_type r = value_traits::zero();
        const coord_type& res = this->grid().resolution();
        const int width = static_cast<int>(m_kernel.size);
        for (int i = 1-width; i <= width; ++i) {
            int ii = kernels::fix_index(c[0] + i, res[0]);
            scalar_type wi = m_kernel(u[0] - i, dorder[0]);
            r +=  wi * (*this)[ii];
        }
        return r;
    }

    value_type convolve2D(const pos_type &u, const coord_type &c, 
                          const coord_type &dorder = coord_type(0)) 
                          const {
        value_type r = value_traits::zero();
        const coord_type& res = this->grid().resolution();
        const int width = static_cast<int>(m_kernel.size);
        for (int j = 1 - width; j <= width; ++j)
        {
            int jj = kernels::fix_index(j + c[1], res[1]);
            scalar_type wj = m_kernel(u[1] - j, dorder[1]);
            for (int i = 1 - width; i <= width; ++i)
            {
                int ii = kernels::fix_index(i + c[0], res[0]);
                scalar_type wij = wj * m_kernel(u[0] - i, dorder[0]);
                r += wij * (*this)[this->m_grid.index(ii, jj)];
            }
        }
        return r;
    }

    value_type convolve3D(const pos_type &u, const coord_type &c, 
                          const coord_type &dorder = coord_type(0)) 
                          const {
        value_type r = value_traits::zero();
        const coord_type& res = this->grid().resolution();
        const int width = static_cast<int>(m_kernel.size);
        for (int k = 1-width; k <= width; ++k)
        {
            int kk = kernels::fix_index(k+c[2], res[2]);
            scalar_type wk = m_kernel(u[2] - k, dorder[2]);
            for (int j = 1-width; j <= width; ++j)
            {
                int jj = kernels::fix_index(j+c[1], res[1]);
                scalar_type wjk = wk * m_kernel(u[1] - j, dorder[1]);
                for (int i = 1-width; i <= width; ++i)
                {
                    int ii = kernels::fix_index(i+c[0], res[0]);
                    scalar_type wijk = wjk * m_kernel(u[0] - i, dorder[0]);
                    r += wijk * (*this)[this->m_grid.index(ii, jj, kk)];
                }
            }
        }
        return r;
    }

    value_type convolve4D(const pos_type &u, const coord_type &c,
                          const coord_type &dorder = coord_type(0)) 
                          const
    {
        value_type r = value_traits::zero();
        const coord_type& res = this->grid().resolution();
        const int width = static_cast<int>(m_kernel.size);
        for (int l = 1 - width; l <= width; ++l)
        {
            int ll = kernels::fix_index(l + c[3], res[3]);
            scalar_type wl = m_kernel(u[3] - l, dorder[3]);
            for (int k = 1 - width; k <= width; ++k)
            {
                int kk = kernels::fix_index(k + c[2], res[2]);
                scalar_type wkl = wl * m_kernel(u[2] - k, dorder[2]);
                for (int j = 1 - width; j <= width; ++j)
                {
                    int jj = kernels::fix_index(j + c[1], res[1]);
                    scalar_type wjkl = wkl * m_kernel(u[1] - j, dorder[1]);
                    for (int i = 1 - width; i <= width; ++i)
                    {
                        int ii = kernels::fix_index(i + c[0], res[0]);
                        scalar_type wijkl = wjkl * m_kernel(u[0] - i, dorder[0]);
                        r += wijkl * (*this)[this->m_grid.index(ii, jj, kk, ll)];
                    }
                }
            }
        }
        return r;
    }

    value_type convolve(const pos_type &u, const coord_type &c, 
                        const coord_type &dorder = coord_type(0),
                        size_type n = dimension) const {
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

public:
    image(const kernel_type &k = kernel_type()) : m_kernel(k) {}
    image(const self_type &other) = default;
    template<typename K_>
    image(const matching_image_type<K_>& other, 
          const kernel_type& k = kernel_type())
        : base_type(other), m_kernel(k) {}
    image(const grid_type &grid, const kernel_type &k = kernel_type())
        : base_type(grid), m_kernel(k) {}
    image(const grid_type &grid, const std::vector<value_type> &data,
          const kernel_type &k = kernel_type())
        : base_type(grid, data), m_kernel(k) {}
    template <typename _Iterator>
    image(const grid_type &grid, _Iterator begin, _Iterator end,
          const kernel_type &k = kernel_type())
        : base_type(grid, begin, end), m_kernel(k) {}
    image(base_type &data, const kernel_type &k = kernel_type())
        : base_type(data), m_kernel(k) {}

    ~image() {}

    // interpolation
    value_type value(const pos_type &p) const {
        std::pair<coord_type, pos_type> r = this->grid().locate(p);
        // std::cout << "position " << p << " was found in cell " << r.first << " with local coordinates " << r.second << '\n';
        return interpolate(r.second, r.first);
    }
    value_type value_in_voxel(const coord_type &vid, 
                              const pos_type &p) const {
        return interpolate(p, vid);
    }

    // derivative at arbitrary locations
    first_derivative_type derivative(const pos_type &p) const {
        std::pair<coord_type, pos_type> r = this->grid().locate(p);
        return derivative_in_voxel(r.first, r.second);
    }

    first_derivative_type derivative_in_voxel(const coord_type &vid, 
                                        const pos_type &p) const {
        first_derivative_type df;
        coord_type order(0);
        const pos_type& spc = this->grid().spacing();
        for (int d = 0; d < dimension; ++d)
        {
            order[d] = 1;
            auto v = convolve(p, vid, order)/spc[d];
            differential::derivatives<value_type, dimension>::set_dfdxi(d, v, df);
            order[d] = 0;
        }
        return df;
    }

    // second derivative at arbitrary locations
    second_derivative_type second_derivative(const pos_type &p) const {
        std::pair<coord_type, pos_type> r = this->m_grid.locate(p);
        return second_derivative_in_voxel(r.first, r.second);
    }

    second_derivative_type second_derivative_in_voxel(const coord_type &vid,
                                                      const pos_type &p) 
                                                      const {
        second_derivative_type d2f;
        coord_type order(0);
        const pos_type& spc = this->grid().spacing();
        for (int d = 0; d < dimension; ++d)
        {
            order[d] = 1;
            for (int d2=d; d2 < dimension; ++d2) {
                order[d2]+= 1;
                auto v = convolve(p, vid, order) / (spc[d]*spc[d2]);
                differential::derivatives<value_type, dimension>::set_d2fdxij(d, d2, v, d2f);
                if (d2>d)
                    differential::derivatives<value_type, dimension>::set_d2fdxij(d2, d, v, d2f);
                --order[d2];
            }
            order[d] = 0;
        }
        return d2f;
    };

private:
    kernel_type m_kernel;
};

template <typename Value, typename Kernel = kernels::Linear>
using image2d = image<long, double, 2, Value, Kernel>;
template <typename Value, typename Kernel = kernels::Linear>
using image3d = image<long, double, 3, Value, Kernel>;
template <typename Value, typename Kernel = kernels::Linear>
using image4d = image<long, double, 4, Value, Kernel>;

template <typename Value, typename Kernel = kernels::Linear>
using image2f = image<long, float, 2, Value, Kernel>;
template <typename Value, typename Kernel = kernels::Linear>
using image3f = image<long, float, 3, Value, Kernel>;
template <typename Value, typename Kernel = kernels::Linear>
using image4f = image<long, float, 4, Value, Kernel>;


template<typename Container>
void print_range(const Container& c, const std::string& what) {
    std::cout << "the range of " << what << " is: min: " << *std::min_element(c.begin(), c.end())
            << ", max: " << *std::max_element(c.begin(), c.end()) << '\n';
}

template<typename Kernel_, 
         typename Size_, typename Scalar_, size_t Dim, typename Value_,
         typename Coord_ = small_vector<Size_, Dim>,
         typename Pos_ = small_vector<Scalar_, Dim> >
raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_>
image_convolution(
    const raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_>& input, 
    const Kernel_& kernel, const Coord_& orders = Coord_(0))
{
    typedef raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_> raster_t;
    typedef Value_ value_t;
    typedef Pos_ pos_t;
    typedef Coord_ coord_t;
    typedef Scalar_ scalar_t;

    // sampled kernel
    std::vector<scalar_t> h(2*kernel.size+1);
    // "h[-k]" = h[kernel.size-k] 

    raster_t out1(input.grid(), input.begin(), input.end()); // copy input
    raster_t out2(input.grid(), value_t(0));

    raster_t *ptr1 = &out1;
    raster_t *ptr2 = &out2;
    // (X * K^N) = ( ( ( (X * K) * K ) * K ) ... ) * K

    auto res = input.grid().resolution();
    size_t npts = input.grid().size();

    std::atomic<size_t> tbb_progress_counter;
    std::mutex update_progress_mutex;

    // complexity Dim x npts x (2*kernel.size()+1)
    ProgressDisplay progress;
    progress.begin(npts*Dim, "Image convolution", 500);
    tbb_progress_counter = 0;
    long sz = kernel.size; // cast size to signed type for kernel indexing
    for (size_t d=0; d<Dim; ++d) {
        for (long i=-sz; i<=sz; ++i) {
            h[i+sz] = kernel(static_cast<scalar_t>(i), orders[d]);
        }
        tbb::parallel_for(tbb::blocked_range<size_t>(0, npts),
                          [&](tbb::blocked_range<size_t> r) {
            for (size_t n=r.begin(); n!=r.end(); ++n)
            {
                std::unique_lock<std::mutex> lock(update_progress_mutex, std::defer_lock);
                if (lock.try_lock())
                {
                    progress.update(tbb_progress_counter);
                }

                ++tbb_progress_counter;

                auto c = input.grid().coordinates(n);
                for (long k=-sz; k<=sz; ++k) {
                    coord_t u = c;
                    u[d] = std::min<Size_>(std::max<Size_>(u[d]-k, 0), res[d]-1);
                    (*ptr2)[n] += h[sz + k] * (*ptr1)(u);
                }
            }
        });
        // swap images
        std::swap(ptr1, ptr2);
        ptr2->initialize(value_t(0));
    }

    progress.end();
    return *ptr1;
}

} // namespace spurt
