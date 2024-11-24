#pragma once

#include <assert.h>
// STL
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <memory>
#include <iomanip>
#include <math/types.hpp>
#include <misc/meta_utils.hpp>
#include <data/index_shifter.hpp>
#include <math/bounding_box.hpp>

namespace spurt {

template<typename Size_, typename Scalar_, size_t Dim, 
         typename Coord_=small_vector<Size_, Dim>, 
         typename Pos_=small_vector<Scalar_, Dim> >
class raster_grid
{
public:
    static const size_t dimension = Dim;

    typedef Scalar_                                           scalar_type;
    typedef Size_                                             size_type;
    typedef size_t                                            dim_type;
    typedef Coord_                                            coord_type;
    typedef index_shifter<Size_, Dim>                         shifter_type;
    typedef Pos_                                              pos_type;
    typedef Pos_                                              point_type;
    typedef bounding_box<pos_type>                            bounds_type;
    typedef raster_grid<Size_, Scalar_, Dim, Coord_, Pos_>    self_type;

    // constructors
    raster_grid() 
        : m_res(0), m_cell_res(0),
          m_bounds(), m_spacing(0), m_npos(0), m_shifter() {}
    raster_grid(const self_type& other) = default;

    raster_grid(const coord_type& resolution, const bounds_type& bounds,
                bool cell_based=false) {
        auto sz = bounds.max() - bounds.min();
        for (size_type i=0; i<dimension; ++i) {
            m_res[i] = resolution[i];
            if (cell_based)
                m_res[i] = m_res[i] + 1;
            m_cell_res[i] = m_res[i] - 1;
            m_spacing[i] = sz[i] / m_cell_res[i];
            m_npos *= m_res[i];
            m_bounds.min()[i] = bounds.min()[i];
            m_bounds.max()[i] = bounds.max()[i];
        }
    }
    raster_grid(const coord_type& resolution, const pos_type& origin,
                const pos_type& spacing, bool cell_based=false) {

        m_npos = 1;
        for (auto i=0; i<dimension; ++i) {
            auto n = resolution[i];
            auto x = origin[i];
            auto d = spacing[i];

            m_bounds.min()[i] = x;
            m_spacing[i] = d;
            if (cell_based) {
                m_res[i] = n + 1;
                m_cell_res[i] = n;
                m_bounds.max()[i] = x + n*d;
                m_npos *= n+1;
            }
            else {
                m_res[i] = n;
                m_cell_res[i] = n - 1;
                m_bounds.max()[i] = x + (n-1) * d;
                m_npos *= n;
            }
        }
    }

    // vertex access
    pos_type operator()(size_type i) const {
        assert(i < m_res[0]);
        pos_type r(m_bounds.min());
        r[i] +=
            static_cast<scalar_type>(i) * m_spacing[0];
        return r;
    }
    pos_type operator()(size_type i, size_type j) const {
        static_assert(dimension >= 2, "Invalid dimension (<2) in raster_grid::operator()");
        assert(i < m_res[0] && j < m_res[1]);
        pos_type r((*this)(i));
        r[1] += static_cast<scalar_type>(j) * m_spacing[1];
        return r;
    }
    pos_type operator()(size_type i, size_type j, size_type k) const {
        static_assert(dimension >= 3, "Invalid dimension (<3) in raster_grid::operator()");
        assert(i < m_res[0] && 
            j < m_res[1] && 
            k < m_res[2]);
        pos_type r((*this)(i, j));
        r[2] += static_cast<scalar_type>(k) * m_spacing[2];
        return r;
    }
    pos_type operator()(size_type i, size_type j, size_type k,
                        size_type l) const {
        static_assert(dimension >= 4, "Invalid dimension (<4) in raster_grid::operator()");
        assert(i < m_res[0] && 
            j < m_res[1] && 
            k < m_res[2] &&
            l < m_res[3]);
        pos_type r((*this)(i, j, k));
        r[3] += static_cast<scalar_type>(l) * m_spacing[3];
        return r;
    }
    pos_type operator()(const coord_type& ids) const {
        pos_type r;
        for (auto i = 0; i < dimension; ++i)
        {
            assert(ids[i] < m_res[i]);
            auto x = m_bounds.min()[i];
            auto d = m_spacing[i];
            auto n = ids[i];
            r[i] = x + n * d;
        }
        return r;
    }

    const coord_type& resolution() const { return m_res; }
    const bounds_type& bounds() const { return m_bounds; }
    const pos_type& spacing() const { return m_spacing; }
    const size_type& size() const { return m_npos; }

    // coordinates to index
    size_type index(size_type i) const {
        assert(i < m_res[0]);
        return i;
    }
    size_type index(size_type i, size_type j) const {
        assert(dimension >= 2);
        assert(i < m_res[0] &&
               j < m_res[1]);
        return m_shifter(0, i, j);
    }
    size_type index(size_type i, size_type j, size_type k) const {
        assert(dimension >= 3);
        assert(i < m_res[0] &&
               j < m_res[1] &&
               k < m_res[2]);
        return m_shifter(0, i, j, k);
    }
    size_type index(size_type i, size_type j, size_type k, size_type l) const {
        assert(dimension >= 4);
        assert(i < m_res[0] &&
               j < m_res[1] &&
               k < m_res[2] && 
               l < m_res[3]);
        return m_shifter(0, i, j, k, l);
    }
    size_type index(const coord_type& ids) const {
        return m_shifter(0, ids);
    }

    // index to coordinates
    coord_type coordinates(size_type idx) const {
        coord_type r;
        size_type aux;
        for (auto i = 0; i < dimension; ++i)
        {
            aux = idx / m_res[i];
            r[i] = idx % m_res[i];
            std::swap(idx, aux);
        }
        return r;
    }

    // index to vertex
    point_type operator[](size_type idx) const {
        return (*this)(coordinates(idx));
    }

    // point to coordinates
    std::pair<coord_type, pos_type> locate(const pos_type& x) const {
        eps_lexicographical_order eps_comp(1.0e-9);
    
        pos_type y = x - m_bounds.min();
        y /= m_spacing;
        pos_type q;
        coord_type c;
        scalar_type o, _max;
        for (int i = 0; i < dimension; ++i)
        {
            _max = static_cast<scalar_type>(m_cell_res[i]);
            if (y[i] > 0 && y[i] < _max)
            {
                o = floor(y[i]);
                q[i] = y[i] - o;
                c[i] = static_cast<size_type>(o);
            }
            else if (eps_comp.equal(y[i], _max))
            {
                c[i] = m_cell_res[i] - 1;
                q[i] = static_cast<scalar_type>(1);
            }
            else if (eps_comp.equal_zero(y[i]))
            {
                c[i] = 0;
                q[i] = static_cast<scalar_type>(0);
            }
            else
            {
                std::ostringstream os;
                os << std::setprecision(16)
                   << "invalid " << i+1 << "-th coordinate in raster_grid::locate()\n"
                   << "coordinate mapped from " << x[i] << " to " << y[i]
                   << " was tested against: [0, " << _max << "]\n";
                throw std::runtime_error(os.str());
            }
        }
        return std::make_pair(c, q);
    }

    // index manipulation helper
    const shifter_type& shifter() const { return m_shifter; }

private:
    coord_type   m_res;
    coord_type   m_cell_res;
    bounds_type  m_bounds;
    pos_type     m_spacing;
    size_type    m_npos;
    shifter_type m_shifter;
};

// Typical specializations
typedef raster_grid< long, double, 1 > rgrid1d;
typedef raster_grid< long, double, 2 > rgrid2d;
typedef raster_grid< long, double, 3 > rgrid3d;
typedef raster_grid< long, double, 4 > rgrid4d;

typedef raster_grid< long, float, 1 > rgrid1f;
typedef raster_grid< long, float, 2 > rgrid2f;
typedef raster_grid< long, float, 3 > rgrid3f;
typedef raster_grid< long, float, 4 > rgrid4f;

} // namespace spurt