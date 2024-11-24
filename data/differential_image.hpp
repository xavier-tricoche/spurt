#ifndef __RASTER_HPP__
#define __RASTER_HPP__

#include <assert.h>
// STL
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <memory>
#include <array>
#include <raster.hpp>

namespace spurt {

namespace kernels {
struct Linear {
    const size_t size = 1;

    double operator()(double t, size_t dorder=0) const {
        double u = (t<0) ? -t : +t;
        if (u>=1) return 0;
        else if (dorder == 0) return 1.-u;
        else if (dorder == 1) return (t<0) ? -1. : 1.;
        else return 0.;
    }
};

struct MitchellNetravaliBC
{
    const size_t size = 2;
    MitchellNetravaliBC(double B=0, double C=0.5) : m_B(B), m_C(C) {}

    double operator()(double t, size_t dorder=0) const {
        double u = (t<0) ? -t : t;
        if (u>=2)
            return 0.;
        
        if (dorder == 0)
        {
            double usq = u * u;
            double ucu = usq * u;
            if (u<1) {
                return (2. - 1.5 * m_B - m_C)*ucu + 
                    (-3 + 2 * m_B + m_C)*usq + 
                    (1. - m_B/3.);
            }
            else {
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
        else return 0.;
    }

    double m_B, m_C;
};

} // namespace kernels

template <typename Value_, size_t Dim_, typename Type_ = double,
          typename Size_ = size_t, typename Kernel_ = kernels::Linear >
class differential_image : public raster_data<Value_, Dim_, Type_, Size_>
{
public:
    static const Size_ dim = Dim_;
    constexpr Size_ hessian_size = dim*(dim+1)/2;

    typedef raster_data<Value_, Dim_, Type_, Size_>  base_type;

    typedef size_t                               dim_type;
    typedef Type_                                scalar_type;
    typedef Value_                               value_type;
    typedef Size_                                size_type;
    typedef typename base_type::grid_type        grid_type;
    typedef typename base_type::point_type       point_type;
    typedef typename base_type::vec_type         vec_type;
    typedef typename base_type::coord_type       coord_type;
    typedef typename base_type::bounds_type      bounds_type;
    typedef std::array<value_type, dim>          deriv_type;
    typedef std::array<value_type, hessian_size> second_deriv_type;
    typedef typename base_type::iterator         iterator;
    typedef typename base_type::const_iterator   const_iterator;
    typedef Kernel_                              kernel_type;
    typedef differential_image<value_type, dim, scalar_type, size_type, kernel_type> self_type;

private:
    value_type interpolate(const vec_type& x, const coord_type& c) const;

    value_type convolve1D(const vec_type &u, const coord_type &c, const 
                          coord_type &dorder = coord_type::Zero()) const;

    value_type convolve2D(const vec_type &u, const coord_type &c, const 
                          coord_type &dorder = coord_type::Zero()) const;

    value_type convolve3D(const vec_type &u, const coord_type &c, const 
                          coord_type &dorder = coord_type::Zero()) const;

    value_type convolve4D(const vec_type &u, const coord_type &c, const 
                          coord_type &dorder = coord_type::Zero()) const;

    value_type convolve(const vec_type &u, const coord_type& c, const 
                        coord_type& dorder = coord_type::Zero(), 
                        size_type n=dim) const;

public:
    differential_image(const kernel_type& k = kernel_type()) : m_kernel(k) {}
    differential_image(const self_type& other) = default;
    differential_image(const grid_type &grid, 
                       const kernel_type &k = kernel_type()) 
        : base_type(grid), m_kernel(k) {}
    differential_image(const grid_type &grid, 
                       const std::vector<value_type> &data, 
                       const kernel_type &k = kernel_type())
        : base_type(grid, data), m_kernel(k) {}
    template <typename _Iterator>
    differential_image(const grid_type &grid, _Iterator begin, _Iterator end, 
                       const kernel_type &k = kernel_type())
        : base_type(grid, begin, end), m_kernel(k) {}
    differential_image(const base_type &data, 
                       const kernel_type &k = kernel_type()) 
        : base_type(data), m_kernel(k) {}

    ~differential_image() {}

    // interpolation
    value_type value(const point_type& p) const;
    value_type value_in_voxel(const coord_type& vid, const point_type& p) const;

    // derivative at arbitrary locations
    deriv_type derivative(const point_type& p) const ;
    deriv_type derivative_in_voxel(const coord_type &vid, const point_type &p) const;

    // second derivative at arbitrary locations
    second_deriv_type derivative(const point_type &p) const;
    second_deriv_type derivative_in_voxel(const coord_type &vid, 
                                          const point_type &p) const;

    kernel_type m_kernel;
};

template<typename Val_, typename Kernel=kernels::Linear>
using diff_image1d = differential_image<Val_, 1, double, Kernel>;
template <typename Val_, typename Kernel = kernels::Linear>
using diff_image2d = differential_image<Val_, 2, double, Kernel>;
template<typename Val_, typename Kernel=kernels::Linear>
using diff_image3d = differential_image<Val_, 3, double, Kernel>;
template<typename Val_, typename Kernel=kernels::Linear>
using diff_image4d = differential_image<Val_, 4, double, Kernel>;

template <typename Val_, typename Kernel = kernels::Linear>
using diff_image1f = differential_image<Val_, 1, float, Kernel>;
template<typename Val_, typename Kernel=kernels::Linear>
using diff_image2f = differential_image<Val_, 2, float, Kernel>;
template<typename Val_, typename Kernel=kernels::Linear>
using diff_image3f = differential_image<Val_, 3, float, Kernel>;
template <typename Val_, typename Kernel = kernels::Linear>
using diff_image4f = differential_image<Val_, 4, float, Kernel>;

} // namespace spurt

#include "detail/differential_image.hpp"

#endif
