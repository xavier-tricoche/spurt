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
#include <math/types.hpp>
#include <misc/meta_utils.hpp>

#include <image/nrrd_wrapper.hpp>

namespace spurt {

template< typename CoordArray_, typename SizeArray_=CoordArray_ >
long coord_to_index(const CoordArray_& coord, const SizeArray_& size);

template< typename SizeArray_, typename CoordArray_=SizeArray_ >
CoordArray_ index_to_coord(long idx, const SizeArray_& size);

// namespace spurt {
template<size_t Dim_, typename Size_>
class index_shifter
{
public:
    static const size_t dimension = Dim_;
    typedef Size_                      size_type;
    typedef Eigen::Vector<Size_, Dim_> array_type;
    typedef Eigen::Vector<Size_, Dim_> coord_type;

    index_shifter(const array_type& sizes = array_type(0));
    index_shifter(const index_shifter& other);

    // move index along a dimension
    size_type lift(size_type base, size_type dim) const;

    // return index of cell vertex, expressed in local coordinates {0,1}^n
    // 1D special case
    size_type operator()(size_type base, size_type i) const;
    // 2D special case
    size_type operator()(size_type base, size_type i, size_type j) const;
    // 3D special case
    size_type operator()(size_type base, size_type i, size_type j,
                         size_type k) const;
    // 4D special case
    size_type operator()(size_type base, size_type i, size_type j,
                         size_type k, size_type l) const;
    // general case
    size_type operator()(size_type base, const coord_type& c) const;
    size_type index(size_type base, size_type local_id) const;

    // per-cell index shift of all the vertices
    const std::vector<size_t>& cell_shifts(const size_type& base) const {
        return m_cell_shifts;
    }

    // map local index to local coordinates
    coord_type local_coords(size_t local_id) const;

private:
    array_type             m_sizes;       // raster dimensions
    array_type             m_dim_shifts;  // per-dimension index shift
    std::vector<size_type> m_cell_shifts; // inner-cell index shifts

    friend std::ostream& operator<<(std::ostream& os, const index_shifter& is) {
        os  << "[ [sizes: " << is.m_sizes
            << "]\n[dim shifts: "
            << is.m_dim_shifts
            << "]\n[cell shifts: ";
        std::copy(is.m_cell_shifts.begin(), is.m_cell_shifts.end(),
                  std::ostream_iterator<size_type>(os, ", "));
        os << "] ]";
        return os;
    }
};

template<size_t Dim_, typename Scalar_ = double,
         typename Size_ = size_t>
class raster_grid
{
public:
    static const size_t dim = Dim_;

    typedef Scalar_                                   scalar_type;
    typedef Size_                                     size_type;
    typedef size_t                                    dim_type;
    typedef index_shifter<dim, size_type>             shifter_type;
    typedef Eigen::Vector<Size_, Dim_>                coord_type;
    typedef Eigen::Vector<scalar_type, dim>           vec_type;
    typedef Eigen::Vector<scalar_type, dim>           point_type;
    typedef bounding_box<scalar_type, dim>            bounds_type;
    typedef raster_grid<dim, scalar_type, size_type>  self_type;

    // constructors
    raster_grid();
    raster_grid(const self_type& other) = default;
    raster_grid(const coord_type& resolution, const bounds_type& bounds,
                bool cell_based=false);
    raster_grid(const coord_type& resolution, const point_type& origin,
                const vec_type& spacing, bool cell_based=false);
    // vertex access
    point_type operator()(size_type i) const;
    point_type operator()(size_type i, size_type j) const;
    point_type operator()(size_type i, size_type j, size_type k) const;
    point_type operator()(size_type i, size_type j, size_type k,
                          size_type l) const;
    point_type operator()(const coord_type& ids) const;

    const coord_type& resolution() const { return m_res; }
    const bounds_type& bounds() const { return m_bounds; }
    const vec_type& spacing() const { return m_spacing; }
    const size_type& size() const { return m_npos; }

    // coordinates to index
    size_type index(size_type i) const;
    size_type index(size_type i, size_type j) const;
    size_type index(size_type i, size_type j, size_type k) const;
    size_type index(size_type i, size_type j, size_type k, size_type l) const;
    size_type index(const coord_type& ids) const;

    // index to coordinates
    coord_type coordinates(size_type idx) const;

    // index to vertex
    point_type operator[](size_type idx) const;

    // point to coordinates
    std::pair<coord_type, vec_type> locate(const point_type& x) const;

    // index manipulation helper
    const shifter_type& shifter() const { return m_shifter; }

private:
    coord_type   m_res;
    coord_type   m_cell_res;
    bounds_type  m_bounds;
    vec_type     m_spacing;
    size_type    m_npos;
    shifter_type m_shifter;
};

// Typical specializations
typedef raster_grid<1, double> rgrid1d;
typedef raster_grid<2, double> rgrid2d;
typedef raster_grid<3, double> rgrid3d;
typedef raster_grid<4, double> rgrid4d;

typedef raster_grid<1, float> rgrid1f;
typedef raster_grid<2, float> rgrid2f;
typedef raster_grid<3, float> rgrid3f;
typedef raster_grid<4, float> rgrid4f;


template<typename Value_, size_t Dim_, typename Scalar_ = double,
         typename Size_ = size_t>
class raster_data
{
public:
    static const Size_ dim = Dim_;

    typedef size_t                                             dim_type;
    typedef Scalar_                                            scalar_type;
    typedef Value_                                             value_type;
    typedef typename std::vector<value_type>::reference        reference_type;
    typedef typename std::vector<value_type>::const_reference  const_reference_type;
    typedef Size_                                              size_type;
    typedef raster_grid<dim, scalar_type, size_type>           grid_type;
    typedef typename grid_type::point_type                     point_type;
    typedef typename grid_type::vec_type                       vec_type;
    typedef typename grid_type::coord_type                     coord_type;
    typedef typename grid_type::bounds_type                    bounds_type;
    typedef typename std::vector<value_type>::iterator         iterator;
    typedef typename std::vector<value_type>::const_iterator   const_iterator;
    typedef raster_data<value_type, dim, scalar_type, size_type> self_type;

    // constructors
    raster_data() {}
	raster_data(const self_type& other) = default;
    raster_data(const grid_type& grid);
    raster_data(const grid_type& grid, const std::vector<value_type>& data);
    raster_data(const grid_type& grid, value_type init_val);
    template<typename _Iterator>
    raster_data(const grid_type& grid, _Iterator begin, _Iterator end);

    ~raster_data() {}

    const grid_type& grid() const { return m_grid; }

    iterator begin() { return m_data.begin(); }
    const_iterator begin() const { return m_data.begin(); }
    iterator end() { return m_data.end(); }
    const_iterator end() const { return m_data.end(); }

    const size_type& size() const { return m_grid.size(); }

    // lookup
    const_reference_type operator()(size_type i) const {
        return m_data[m_grid.index(i)];
    }
    const_reference_type operator()(size_type i, size_type j) const {
        return m_data[m_grid.index(i, j)];
    }
    const_reference_type operator()(size_type i, size_type j, size_type k) const {
        return m_data[m_grid.index(i, j, k)];
    }
    const_reference_type operator()(size_type i, size_type j, size_type k,
                                 size_type l) const {
        return m_data[m_grid.index(i, j, k, l)];
    }
    const_reference_type operator()(const coord_type& id) const {
        return m_data[m_grid.index(id)];
    }
    // array-type access for convenience
    const_reference_type operator[](size_type idx) const {
        return m_data[idx];
    }

    reference_type operator()(size_type i) {
        return m_data[m_grid.index(i)];
    }
    reference_type operator()(size_type i, size_type j) {
        return m_data[m_grid.index(i, j)];
    }
    reference_type operator()(size_type i, size_type j, size_type k) {
        return m_data[m_grid.index(i, j, k)];
    }
    reference_type operator()(size_type i, size_type j, size_type k,
                           size_type l) {
        return m_data[m_grid.index(i, j, k, l)];
    }
    reference_type operator()(const coord_type& id) {
        return m_data[m_grid.index(id)];
    }
    reference_type operator[](size_type idx) {
        return m_data[idx];
    }

    void save_as_nrrd(const std::string& filename) const;

protected:
    typedef typename grid_type::shifter_type shifter_type;
    grid_type                 m_grid;
    std::vector<value_type>   m_data;
};

template<typename _Value>
using raster1d = raster_data<_Value, 1, double>;
template<typename _Value>
using raster2d = raster_data<_Value, 2, double>;
template<typename _Value>
using raster3d = raster_data<_Value, 3, double>;
template<typename _Value>
using raster4d = raster_data<_Value, 4, double>;

template<typename _Value>
using raster1f = raster_data<_Value, 1, float>;
template<typename _Value>
using raster2f = raster_data<_Value, 2, float>;
template<typename _Value>
using raster3f = raster_data<_Value, 3, float>;
template<typename _Value>
using raster4f = raster_data<_Value, 4, float>;

template<typename Value_, size_t Dim_, typename Type_ = double,
         typename Size_ = size_t>
class image : public raster_data<Value_, Dim_, Type_, Size_>
{
public:
    static const Size_ dim = Dim_;

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
    typedef nvis::fixed_vector<value_type, dim>  deriv_type;
    typedef nvis::fixed_vector<value_type, 1>    val1;
    typedef nvis::fixed_vector<value_type, 2>    val2;
    typedef nvis::fixed_vector<value_type, 3>    val3;
    typedef nvis::fixed_vector<value_type, 4>    val4;
    typedef typename base_type::iterator         iterator;
    typedef typename base_type::const_iterator   const_iterator;
    typedef image<value_type, dim, scalar_type, size_type> self_type;

private:
    // interpolation
    // 1D case
    value_type linear(const vec_type& u, size_type a) const;
    // 2D case
    value_type bilinear(const vec_type& u, size_type a) const ;
    // 3D case
    value_type trilinear(const vec_type& u, size_type a) const;
    // 4D case
    value_type quadrilinear(const vec_type& u, size_type a) const;
    // General case
    value_type multilinear(const vec_type& u, size_type a,
                           size_t n=dim) const;
    value_type interpolate(const vec_type& x, const coord_type& c) const;

    // derivative within elements
    deriv_type dmultilinear(const vec_type& u, size_type a,
                            size_type n=dim) const;
    value_type dlinear_part(const vec_type& u, size_type a) const;
    value_type dbilinear_part(const vec_type& u, size_type a,
                              size_type d) const ;
    value_type dtrilinear_part(const vec_type& u, size_type a,
                               size_type d) const;
    value_type dquadrilinear_part(const vec_type& u, size_type a,
                                  size_type d) const;
    value_type dmultilinear_part(const vec_type& u, size_type a, size_type d,
                                 size_type n=dim) const;

    // derivative at vertices
    value_type central_diff(size_type i, size_type id, size_type dim) const;

    // helper functions
    value_type _value(size_type base, size_type loc_index) const {
        return this->m_data[this->_grid.shifter().index(base, loc_index)];
    }
    value_type _value(size_type base, const coord_type& c) const {
        return this->m_data[this->_grid.shifter()(base, c)];
    }
    size_type _lift(size_type base, size_type d) const {
        return this->m_grid.shifter().lift(base, d);
    }

public:
    image() {}
    image(const self_type& other) = default;
    image(const grid_type& grid) : base_type(grid) {}
    image(const grid_type& grid, const std::vector<value_type>& data)
        : base_type(grid, data) {}
    template<typename _Iterator>
    image(const grid_type& grid, _Iterator begin, _Iterator end)
        : base_type(grid, begin, end) {}
    image(const base_type& data) : base_type(data) {}

    ~image() {}

    // interpolation
    value_type value(const point_type& p) const;
    value_type value_in_voxel(const coord_type& vid, const point_type& p) const;

    // derivative at arbitrary locations
    deriv_type derivative(const point_type& p) const ;
    deriv_type derivative_in_voxel(const coord_type& vid, const point_type& p) const;

    // derivative at vertices
    val1 derivative(size_type i) const ;
    val2 derivative(size_type i, size_type j) const;
    val3 derivative(size_type i, size_type j, size_type k) const;
    val4 derivative(size_type i, size_type j, size_type k, size_type l) const;
    deriv_type derivative(const coord_type& id) const;
};

template<typename Val_>
using image1d = image<Val_, 1, double>;
template<typename Val_>
using image2d = image<Val_, 2, double>;
template<typename Val_>
using image3d = image<Val_, 3, double>;
template<typename Val_>
using image4d = image<Val_, 4, double>;

template<typename Val_>
using image1f = image<Val_, 1, float>;
template<typename Val_>
using image2f = image<Val_, 2, float>;
template<typename Val_>
using image3f = image<Val_, 3, float>;
template<typename Val_>
using image4f = image<Val_, 4, float>;


} // namespace spurt

#include "detail/raster.hpp"

#endif
