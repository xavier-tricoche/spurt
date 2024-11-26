#pragma once

#include <array>
#include <cmath>
#include <math/types.hpp>
#include <data/raster_data.hpp>


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
                return (t < 0) ? -1. : 1.;
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

template <typename Size_, typename Scalar_, size_t Dim, typename Value_,
          typename Kernel_=kernels::Linear, typename Coord_ = small_vector<Size_, Dim>,
          typename Pos_ = small_vector<Scalar_, Dim> >
class image : public raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_>
{
public:
    typedef raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_> base_type;
    typedef typename base_type::size_type size_type;
    static const size_type dimension = base_type::dimension;
    static const size_type hess_dim = dimension*(dimension+1)/2;

    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::pos_type pos_type;
    typedef typename base_type::coord_type coord_type;
    typedef typename base_type::bounds_type bounds_type;
    typedef typename base_type::grid_type grid_type;
    typedef typename base_type::value_type value_type;
    typedef data_traits<value_type> value_traits;

    typedef std::array<value_type, dimension> derivative_type;
    typedef std::array<value_type, hess_dim>   second_derivative_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef image<Size_, Scalar_, Dim, Value_, Kernel_, Coord_, Pos_> self_type;
    typedef Kernel_ kernel_type;

private:
    value_type interpolate(const pos_type &x, const coord_type &c) const {
        return convolve(x, c);
    }

    value_type convolve1D(const pos_type &u, const coord_type &c, 
                          const coord_type &dorder = coord_type(0)) 
                          const {  
        value_type r = value_traits::zero();
        scalar_type allw = 0;
        const coord_type& res = this->grid().resolution();
        const coord_type& coord = c;
        const coord_type& order = dorder;
        const pos_type& pos = u;

        const int width = static_cast<int>(m_kernel.size);
        for (int i = 1-width; i <= width; ++i) {
            int ii = kernels::fix_index(coord[0] + i, res[0]);
            scalar_type wi = m_kernel(pos[0] - i, dorder[0]);
            allw += wi;
            r +=  wi * this->m_data[ii];
        }
        return r/allw;
    }

    value_type convolve2D(const pos_type &u, const coord_type &c, 
                          const coord_type &dorder = coord_type(0)) 
                          const {
        value_type r = value_traits::zero();
        scalar_type allw = 0;
        const coord_type& res = this->grid().resolution();
        const coord_type& coord = c;
        const coord_type& order = dorder;
        const pos_type& pos = u;
        const int width = static_cast<int>(m_kernel.size);
        for (int j = 1 - width; j <= width; ++j)
        {
            int jj = kernels::fix_index(j + coord[1], res[1]);
            scalar_type wj = m_kernel(pos[1] - j, order[1]);
            for (int i = 1 - width; i <= width; ++i)
            {
                int ii = kernels::fix_index(i + coord[0], res[0]);
                scalar_type wij = wj * m_kernel(pos[0] - i, order[0]);
                allw += wij;
                r += wij * this->m_data[this->m_grid.index(ii, jj)];
            }
        }
        return r / allw;
    }

    value_type convolve3D(const pos_type &u, const coord_type &c, 
                          const coord_type &dorder = coord_type(0)) 
                          const {
        bool verbose = false;
        if (verbose) std::cout << "convolve3d: u=" << u << ", c=" << c << ", dorder=" << dorder << '\n';
        if (verbose) std::cout << "width = " << m_kernel.size << '\n';
        value_type r = value_traits::zero();
        scalar_type allw = 0;
        const coord_type& res = this->grid().resolution();
        if (verbose) std::cout << "res = " << res << '\n';
        const coord_type& coord = c;
        const coord_type& order = dorder;
        const pos_type& pos = u;
        const int width = static_cast<int>(m_kernel.size);
        for (int k = 1-width; k <= width; ++k)
        {
            if (verbose) std::cout << "k=" << k << '\n';
            int kk = kernels::fix_index(k+coord[2], res[2]);
            if (verbose) std::cout << "kk=" << kk << '\n';
            scalar_type wk = m_kernel(pos[2] - k, order[2]);
            if (verbose) std::cout << "kernel(" << pos[2]-k << ", " << order[2] << ")=" << wk << '\n';
            for (int j = 1-width; j <= width; ++j)
            {
                if (verbose) std::cout << "j=" << j << '\n';
                int jj = kernels::fix_index(j+coord[1], res[1]);
                if (verbose) std::cout << "jj=" << jj << '\n';
                scalar_type wjk = wk * m_kernel(pos[1] - j, order[1]);
                if (verbose) std::cout << "kernel(" << pos[1]-j << ", " << order[1] << ") = " << m_kernel(pos[1]-j, order[1]) << '\n';
                for (int i = 1-width; i <= width; ++i)
                {
                    if (verbose) std::cout << "i=" << i << '\n';
                    int ii = kernels::fix_index(i+coord[0], res[0]);
                    if (verbose) std::cout << "ii=" << ii << '\n';
                    scalar_type wijk = wjk * m_kernel(pos[0] - i, order[0]);
                    if (verbose) std::cout << "kernel(" << pos[0]-i << ", " << order[0] << ")=" << wijk << '\n';
                    allw += wijk;
                    if (verbose) std::cout << "allw(" << i << ", " << j << ", " << k << ", " << order << ")=" << allw << '\n';
                    r += wijk * this->m_data[this->m_grid.index(ii, jj, kk)];
                    if (verbose) std::cout << "r=" << r << '\n';
                }
            }
        }
        if (verbose) {
            std::cout << "Returning " << r << " / " << allw << " = " << r/allw << '\n';
            if (std::abs(allw) < 1.0e-9) std::cout << "WARNING!!!\n";
        }
        return r; // / allw;
    }

    value_type convolve4D(const pos_type &u, const coord_type &c,
                          const coord_type &dorder = coord_type(0)) 
                          const
    {
        value_type r = value_traits::zero();
        scalar_type allw = 0;
        const coord_type&  res = this->grid().resolution();
        const coord_type&  coord = c;
        const coord_type& order = dorder;
        const pos_type&  pos = u;
        const int width = static_cast<int>(m_kernel.size);
        for (int l = 1 - width; l <= width; ++l)
        {
            int ll = kernels::fix_index(l + coord[3], res[3]);
            scalar_type wl = m_kernel(pos[3] - l, order[3]);
            for (int k = 1 - width; k <= width; ++k)
            {
                int kk = kernels::fix_index(k + coord[2], res[2]);
                scalar_type wkl = wl * m_kernel(pos[2] - k, order[2]);
                for (int j = 1 - width; j <= width; ++j)
                {
                    int jj = kernels::fix_index(j + coord[1], res[1]);
                    scalar_type wjkl = wkl * m_kernel(pos[1] - j, order[1]);
                    for (int i = 1 - width; i <= width; ++i)
                    {
                        int ii = kernels::fix_index(i + coord[0], res[0]);
                        scalar_type wijkl = wjkl * m_kernel(pos[0] - i, order[0]);
                        allw += wijkl;
                        r += wijkl * this->m_data[this->m_grid.index(ii, jj, kk, ll)];
                    }
                }
            }
        }
        return r / allw;
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
    image(const grid_type &grid, const kernel_type &k = kernel_type())
        : base_type(grid), m_kernel(k) {}
    image(const grid_type &grid, const std::vector<value_type> &data,
          const kernel_type &k = kernel_type())
        : base_type(grid, data), m_kernel(k) {}
    template <typename _Iterator>
    image(const grid_type &grid, _Iterator begin, _Iterator end,
          const kernel_type &k = kernel_type())
        : base_type(grid, begin, end), m_kernel(k) {}
    image(const base_type &data,
          const kernel_type &k = kernel_type())
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
    derivative_type derivative(const pos_type &p) const {
        std::pair<coord_type, pos_type> r = this->grid().locate(p);
        return derivative_in_voxel(r.first, r.second);
    }

    derivative_type derivative_in_voxel(const coord_type &vid, 
                                        const pos_type &p) const {
        derivative_type df;
        coord_type order(0);
        for (int d = 0; d < dimension; ++d)
        {
            order[d] = 1;
            df[d] = convolve(p, vid, order);
            df[d] /= this->grid().spacing()[d];
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
        int offset=0;
        for (int d = 0; d < dimension; ++d)
        {
            order[d] = 2; // diagonal
            d2f[offset] = convolve(p, vid, order);
            d2f[offset] /= this->grid().spacing()[d] * 
                           this->grid().spacing()[d];
            order[d] = 1;
            for (int dd=d+1; dd<dimension; ++dd) {
                order[dd] = 1;
                offset += 1;
                d2f[offset] = convolve(p, vid, order);
                d2f[offset] /= this->grid().spacing()[d]* 
                               this->grid().spacing()[dd];
                order[dd] = 0;
            }
            order[d] = 0;
        }
        return d2f;
    }

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

} // namespace spurt
