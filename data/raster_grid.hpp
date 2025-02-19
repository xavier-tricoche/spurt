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
    typedef Pos_                                              pos_type;
    typedef Pos_                                              point_type;
    typedef bounding_box<pos_type>                            bounds_type;
    typedef raster_grid<Size_, Scalar_, Dim, Coord_, Pos_>    self_type;

    // constructors
    raster_grid() 
        : m_res(0), m_cell_res(0),
          m_bounds(), m_spacing(0), m_npos(0) {}
    raster_grid(const self_type& other) = default;

    raster_grid(const coord_type& resolution, const bounds_type& bounds,
                bool cell_based=false) : m_res(resolution) {
        auto sz = bounds.max() - bounds.min();
        m_npos = product(resolution);
        m_bounds.min() = bounds.min();
        m_bounds.max() = bounds.max();
        m_cell_res = m_res - 1;
        if (cell_based) m_cell_res += 1;
        m_spacing = sz / m_cell_res;

#ifdef SPURT_DEBUG
        std::cout << "grid constructor: npos = " << m_npos << ", bounds = " << m_bounds << ", cell res = " << m_cell_res << ", size = " << sz << ", spacing = " << m_spacing << '\n'; 
#endif
    }
    raster_grid(const coord_type& resolution, const pos_type& origin,
                const pos_type& spacing, bool cell_based=false)
        : m_res(resolution) {
        m_npos = product(m_res);
        m_bounds.min() = origin;
        m_spacing = spacing;
        m_cell_res = m_res-1;
        if (cell_based) m_cell_res += 1;
        m_bounds.max() = m_bounds.min() + m_spacing*m_cell_res;
    }

    // vertex access
    pos_type operator()(size_type i) const {
        assert(i < m_res[0]);
        pos_type r(m_bounds.min());
        r[0] += static_cast<scalar_type>(i) * m_spacing[0];
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
        assert(i < m_res[0] && j < m_res[1] && k < m_res[2]);
        pos_type r((*this)(i, j));
        r[2] += static_cast<scalar_type>(k) * m_spacing[2];
        return r;
    }
    pos_type operator()(size_type i, size_type j, size_type k, size_type l) const {
        static_assert(dimension >= 4, "Invalid dimension (<4) in raster_grid::operator()");
        assert(i < m_res[0] && j < m_res[1] && k < m_res[2] && l < m_res[3]);
        pos_type r((*this)(i, j, k));
        r[3] += static_cast<scalar_type>(l) * m_spacing[3];
        return r;
    }
    pos_type operator()(const coord_type& ids) const {
        pos_type r = m_bounds.min();
        for (auto i = 0; i < dimension; ++i)
        {
            assert(ids[i] < m_res[i]);
            r[i] += static_cast<scalar_type>(ids[i]) * m_spacing[i];
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
        return i + m_res[0]*j;
    }
    size_type index(size_type i, size_type j, size_type k) const {
        assert(dimension >= 3);
        assert(i < m_res[0] &&
               j < m_res[1] &&
               k < m_res[2]);
        return i + m_res[0]*(j + m_res[1]*k);
    }
    size_type index(size_type i, size_type j, size_type k, size_type l) const {
        assert(dimension >= 4);
        assert(i < m_res[0] &&
               j < m_res[1] &&
               k < m_res[2] && 
               l < m_res[3]);
        return i + m_res[0]*(j + m_res[1]*(k + m_res[2]*l));
    }
    size_type index(const coord_type& ids) const {
        size_type r = ids[dimension-1];
        for (int i=dimension-2; i>=0; --i) {
            r *= m_res[i];
            r += ids[i];
        }
        return r;
    }
    
    // index(i,j,k,l,m,n):
    // dim = 6
    // r = n
    // r = r * res[4] + m = m + res[4] * n
    // r = r * res[3] + l = l + res[3] * ( m + res[4] * n )
    // r = r * res[2] + k = k + res[2] * ( l + res[3]*( m + res[4] * n ) )
    // r = r * res[1] + j = j + res[1] * ( k + res[2]*( l + res[3]*( m + res[4] * n ) ) )
    // r = r * res[0] + i = i + res[0] * ( j + res[1]*( k + res[2]*( l + res[3]*(m + res[4] * n ) ) ) )

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
        const double eps = 1.0e-9;
        bool verbose = false;
        
        auto d = m_bounds.inf_distance(x);
        if (d > eps) {
            std::ostringstream os;
            os << "invalid query position " << x << '\n';
            os << "bounds are " << m_bounds << '\n';
            os << "outer distance is " << d << '\n';
            throw std::runtime_error(os.str());
        }
        else
        { 
            pos_type y = x - m_bounds.min();
            if (verbose) std::cout << "shifted position of " << x << " is " << y << '\n';
            y /= m_spacing;
            if (verbose) std::cout << "after normalization by spacing: " << y << '\n';
            coord_type c = floor(y);
            if (verbose) std::cout << "coordinates = " << c << '\n';
            pos_type q = y - c;
            if (verbose) std::cout << "local coordinates = " << q << '\n';
        
            for (int i=0; i<dimension; ++i) {
                if (c[i] == m_cell_res[i])
                {
                    if (verbose) std::cout << "cell coordinate in dim " << i << " is out of range (" << c[i] << ")\n";
                    c[i] = m_cell_res[i]-1;
                    if (verbose) std::cout << "it is corrected to " << c[i] << "\n";
                    q[i] = 1;
                }
                else if (y[i] < 0)
                {
                    if (verbose) std::cout << "position coordinate in dim " << i << " is negative (" << y[i] << ")\n";
                    c[i] = 0;
                    q[i] = 0;
                }
            }
            return std::make_pair(c, q);
        }
        // should not be reached...
        return std::make_pair(coord_type(), pos_type());
    }

private:
    coord_type   m_res;
    coord_type   m_cell_res;
    bounds_type  m_bounds;
    pos_type     m_spacing;
    size_type    m_npos;
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