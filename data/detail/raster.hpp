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

namespace xavier {

template< typename CoordArray_, typename SizeArray_ >
long coord_to_index(const CoordArray_& coord, const SizeArray_& size) {
    size_t N = size.size();

    // i + size0*(j + size1*(k + size2*(l + ....)))
    long idx=coord[N-1];
    for (long dim=N-2; dim>=0; --dim) {
        idx = coord[dim] + idx*size[dim];
    }
    return idx;
}

template< typename SizeArray_, typename CoordArray_ >
CoordArray_ index_to_coord(long idx, const SizeArray_& size) {
    CoordArray_ coord;
    typedef typename CoordArray_::value_type value_t;

    // i + size0*( j + size1*( k + size2*( l ) ) )
    size_t N = size.size();
    for (size_t i=0; i<N; ++i) {
        std::ldiv_t qr = std::div(idx, size[i]);
        coord[i] = static_cast<value_t>(qr.rem);
        idx = qr.quot;
    }
    return coord;
}


/*****************************************************************************

                                index_shifter

*****************************************************************************/

template<size_t Dim_, typename Size_>
inline typename index_shifter<Dim_, Size_>::coord_type
index_shifter<Dim_, Size_>::local_coords(size_t local_id) const {
    coord_type r;
    std::bitset<sizeof(size_type)> bits(local_id);
    for (size_t i=0 ; i<r.size() ; ++i) {
        r[i] = bits[i];
    }
    return r;
}

template<size_t Dim_, typename Size_>
index_shifter<Dim_, Size_>::index_shifter(const array_type& sizes)
    : m_sizes(sizes) {
    size_type shift = 1;
    for (size_t i=0 ; i<dimension ; ++i) {
        m_dim_shifts[i] = shift;
        shift *= m_sizes[i];
    }

    size_t points_per_cell = 1;
    points_per_cell = points_per_cell << dimension;
    m_cell_shifts.resize(points_per_cell);
    std::fill(m_cell_shifts.begin(), m_cell_shifts.end(), 0);

    for (size_t i=1 ; i<points_per_cell ; ++i) {
        std::bitset<sizeof(size_type)> bits(i);
        for (size_t j=0 ; j<dimension ; ++j) {
            if (bits[j]) m_cell_shifts[i] += m_dim_shifts[i];
        }
    }
}

template<size_t Dim_, typename Size_>
index_shifter<Dim_, Size_>::index_shifter(const index_shifter& other)
    : m_sizes(other.m_sizes), m_dim_shifts(other.m_dim_shifts),
      m_cell_shifts(other.m_cell_shifts) {}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::
operator()(size_type base, size_type i) const {
    return base + i;
}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::
operator()(size_type base, size_type i, size_type j) const {
    return base + i + j*m_dim_shifts[1];
}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::
operator()(size_type base, size_type i, size_type j, size_type k) const {
    return base + i + j*m_dim_shifts[1] + k*m_dim_shifts[2];
}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::
operator()(size_type base, size_type i, size_type j, size_type k,
           size_type l) const {
    return base + i + j*m_dim_shifts[1] + k*m_dim_shifts[2] + l*m_dim_shifts[3];
}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::
operator()(size_type base, const coord_type& c) const {
    size_type r(base);
    r += c[0];
    for (size_t i=1 ; i<dimension ; ++i)
        r += c[i]*m_dim_shifts[i];
    return r;
}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::lift(size_type base, size_type d) const {
    return base + m_dim_shifts[d];
}

template<size_t Dim_, typename Size_>
typename index_shifter<Dim_, Size_>::size_type
index_shifter<Dim_, Size_>::index(size_type base, size_type local_id) const {
    return base + m_cell_shifts[local_id];
}

/*****************************************************************************

                                raster_grid

*****************************************************************************/

template<size_t Dim_, typename Scalar_, typename Size_>
inline raster_grid<Dim_, Scalar_, Size_>::
raster_grid(const coord_type& resolution, const bounds_type& bounds,
            bool cell_based)
    : m_res(resolution), m_bounds(bounds), m_shifter(resolution)
{
    static_assert(dim != 0, "Invalid 0 dimension in raster_grid constructor");
    if (cell_based) m_res += coord_type(1);
    vec_type d = m_bounds.size();
    m_npos = 1;
    for (dim_type i = 0 ; i < dim ; ++i) {
        assert(m_res[i] > 1);
        m_spacing[i] = d[i] / static_cast<scalar_type>(m_res[i] - 1);
        m_npos *= m_res[i];
        m_cell_res[i] = m_res[i]-1;
    }
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline raster_grid<Dim_, Scalar_, Size_>::
raster_grid(const coord_type& resolution, const vec_type& origin,
            const vec_type& spacing, bool cell_based)
    : m_res(resolution), m_spacing(spacing), m_shifter(resolution)
{
    static_assert(dim != 0, "Invalid 0 dimension in raster_grid constructor");
    if (cell_based) m_res += coord_type(1);
    m_bounds.min() = origin;
    m_npos = 1;
    for (dim_type i = 0 ; i < dim ; ++i) {
        assert(m_res[i] > 1);
        m_bounds.max()[i] = origin[i] + static_cast<scalar_type>(m_res[i]-1)*m_spacing[i];
        m_npos *= m_res[i];
        m_cell_res[i] = m_res[i]-1;
    }
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::vec_type
raster_grid<Dim_, Scalar_, Size_>::operator()(size_type i) const
{
    assert(i < m_res[0]);
    vec_type r(m_bounds.min());
    r[0] += static_cast<scalar_type>(i)*m_spacing[0];
    return r;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::vec_type
raster_grid<Dim_, Scalar_, Size_>::operator()(size_type i, size_type j)
const
{
    static_assert(dim >= 2, "Invalid dimension (<2) in raster_grid::operator()");
    assert(i < m_res[0] && j < m_res[1]);
    vec_type r((*this)(i));
    r[1] += static_cast<scalar_type>(j)*m_spacing[1];
    return r;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::vec_type
raster_grid<Dim_, Scalar_, Size_>::operator()(size_type i, size_type j,
        size_type k) const
{
    static_assert(dim >= 3, "Invalid dimension (<3) in raster_grid::operator()");
    assert(i < m_res[0] && j < m_res[1] && k < m_res[2]);
    vec_type r((*this)(i,j));
    r[2] += static_cast<scalar_type>(k)*m_spacing[2];
    return r;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::vec_type
raster_grid<Dim_, Scalar_, Size_>::operator()(size_type i, size_type j,
        size_type k, size_type l) const
{
    static_assert(dim >= 4, "Invalid dimension (<4) in raster_grid::operator()");
    assert(i < m_res[0] && j < m_res[1] && k < m_res[2] && l < m_res[3]);
    vec_type r((*this)(i,j,k));
    r[3] += static_cast<scalar_type>(l)*m_spacing[3];
    return r;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::vec_type
raster_grid<Dim_, Scalar_, Size_>::operator()(const coord_type& ids) const
{
    for (dim_type i = 0 ; i < dim ; ++i) {
        assert(ids[i] < m_res[i]);
    }
    return m_bounds.min() + vec_type(ids)*m_spacing;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::size_type
raster_grid<Dim_, Scalar_, Size_>::index(size_type i) const
{
    assert(i < m_res[0]);
    return i;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::size_type
raster_grid<Dim_, Scalar_, Size_>::index(size_type i, size_type j) const
{
    static_assert(dim >= 2, "Invalid dimension (<2) in raster_grid::operator()");
    assert(i < m_res[0] && j < m_res[1]);
    return m_shifter(0, i, j);
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::size_type
raster_grid<Dim_, Scalar_, Size_>::index(size_type i, size_type j,
                                         size_type k) const
{
    static_assert(dim >= 3, "Invalid dimension (<3) in raster_grid::operator()");
    assert(i < m_res[0] && j < m_res[1] && k < m_res[2]);
    return m_shifter(0, i, j, k);
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::size_type
raster_grid<Dim_, Scalar_, Size_>::index(size_type i, size_type j,
                                         size_type k, size_type l) const
{
    static_assert(dim >= 4, "Invalid dimension (<4) in raster_grid::operator()");
    assert(i < m_res[0] && j < m_res[1] && k < m_res[2] && l < m_res[3]);
    return m_shifter(0, i, j, k, l);
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::size_type
raster_grid<Dim_, Scalar_, Size_>::index(const coord_type& ids) const
{
    return m_shifter(0, ids);
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::coord_type
raster_grid<Dim_, Scalar_, Size_>::coordinates(size_type idx) const
{
    coord_type r;
    size_type aux;
    for (dim_type i=0 ; i<dim ; ++i) {
        aux = idx / m_res[i];
        r[i] = idx % m_res[i];
        std::swap(idx, aux);
    }
    return r;
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline typename raster_grid<Dim_, Scalar_, Size_>::point_type
raster_grid<Dim_, Scalar_, Size_>::operator[](size_type idx) const
{
    return (*this)(coordinates(idx));
}

template<size_t Dim_, typename Scalar_, typename Size_>
inline std::pair<typename raster_grid<Dim_, Scalar_, Size_>::coord_type,
                 typename raster_grid<Dim_, Scalar_, Size_>::vec_type>
raster_grid<Dim_, Scalar_, Size_>::locate(const point_type& x) const
{
    typedef eps_compare< scalar_type > eps_comp_t;
    vec_type y = x - m_bounds.min();
    y /= m_spacing;
    vec_type q;
    coord_type c;
    scalar_type o, _max;
    for (int i=0; i<dim; ++i) {
        _max = static_cast<scalar_type>(m_cell_res[i]);
        if ( y[i] > 0 && y[i] < _max ) {
            o = floor(y[i]);
            q[i] = y[i] - o;
            c[i] = static_cast<size_type>(o);
        }
        else if ( eps_comp_t::equal(y[i], _max) ) {
            c[i] = m_cell_res[i]-1;
            q[i] = static_cast<scalar_type>(1);
        }
        else if ( eps_comp_t::equal_zero(y[i]) ) {
            c[i] = 0;
            q[i] = static_cast<scalar_type>(0);
        }
        else {
            std::ostringstream os;
            os << std::setprecision(16)
            << "invalid " << xavier::number_to_rank(i+1) << " coordinate in raster_grid::locate()\n"
            << "coordinate mapped from " << x[i] << " to " << y[i]
            << " was tested against: [0, " << _max << "]\n";
            throw std::runtime_error(os.str());
        }
    }
    return std::make_pair(c, q);
}

/*****************************************************************************

                                raster_data

*****************************************************************************/

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
raster_data<Value_, Dim_, Scalar_, Size_>::raster_data(const grid_type& grid)
    : m_grid(grid), m_data(grid.size())
{
    m_max_coord = vec_type(m_grid.resolution()) - vec_type(1.);
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
raster_data<Value_, Dim_, Scalar_, Size_>::
raster_data(const grid_type& grid, value_type init_val)
: m_grid(grid), m_data(grid.size(), init_val)
{
    m_max_coord = vec_type(m_grid.resolution()) - vec_type(1.);
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
raster_data<Value_, Dim_, Scalar_, Size_>::
raster_data(const grid_type& grid, const std::vector<value_type>& data)
    : m_grid(grid), m_data(data)
{
    assert(data.size() == grid.size());
    m_max_coord = vec_type(m_grid.resolution()) - vec_type(1.);
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
template<typename _Iterator>
raster_data<Value_, Dim_, Scalar_, Size_>::
raster_data(const grid_type& grid, _Iterator begin, _Iterator end)
    : m_grid(grid)
{
    std::copy(begin, end, std::back_inserter(m_data));
    assert(m_data.size() == m_grid.size());
    m_max_coord = vec_type(m_grid.resolution()) - vec_type(1.);
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
raster_data<Value_, Dim_, Scalar_, Size_>::raster_data(const self_type& other)
    : m_grid(other.m_grid), m_data(other.m_data)
{
    m_max_coord = other.m_max_coord;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
raster_data<Value_, Dim_, Scalar_, Size_>::raster_data(self_type&& other)
    : m_grid(other.m_grid), m_data(std::move(other.m_data))
{
    m_max_coord = other.m_max_coord;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
void raster_data<Value_, Dim_, Scalar_, Size_>::
save_as_nrrd(const std::string& filename) const
{
    std::string type_string = xavier::type2string<Scalar_>::type_name();
    size_t ncomp = xavier::data_traits<Value_>::size();

    std::ostringstream os;
    os << "NRRD0001\n";
    os << "# Complete NRRD file format specification at:\n";
    os << "# http://teem.sourceforge.net/nrrd/format.html\n";
    os << "type: " << type_string << '\n';
    os << "dimension: " << (ncomp > 1 ? Dim_+1 : Dim_) << '\n';
    os << "sizes:";
    if (ncomp > 1) os << " " << ncomp;
    for (int i=0; i<Dim_; ++i) {
        os << " " << m_grid.resolution()[i];
    }
    os << '\n';
    os << "spacings:";
    if (ncomp > 1) os << " nan";
    for (int i=0; i<Dim_; ++i) {
        os << " " << std::setprecision(17) << m_grid.spacing()[i];
    }
    os << '\n';
    os << "axis mins:";
    if (ncomp > 1) os << " nan";
    for (int i=0; i<Dim_; ++i) {
        os << " " << std::setprecision(17) << m_grid.bounds().min()[i];
    }
    os << '\n';
    os << "centerings:";
    if (ncomp > 1) os << " ???";
    for (int i=0; i<Dim_; ++i) {
        os << " node";
    }
    os << '\n';
    os << "endian: little\n";
    os << "encoding: raw\n";

    size_t sz = ncomp*m_data.size()*sizeof(Scalar_);
    std::cout << "exporting " << sz << " bytes" << '\n';
    std::cout << "grid res: " << m_grid.resolution() << '\n';
    std::cout << "grid size: " << m_grid.size() << '\n';

    std::ofstream out(filename.c_str(), std::ios::binary);
    out << os.str() << std::endl;
    out.write((char*)&m_data[0], sz);
    out.close();
}

/*****************************************************************************

                                image

*****************************************************************************/

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
inline typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
linear(const vec_type& u, size_type a) const
{
    return (1.-u[0])*this->m_data[a] + u[0]*this->m_data[a+1];
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
inline typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
bilinear(const vec_type& u, size_type a) const
{
    const typename base_type::shifter_type& shift = this->grid().shifter();
    const vec_type& hi = u;
    vec_type lo = vec_type(1.) - hi;
    return
        lo[0]*lo[1]*this->m_data[a           ] +
        hi[0]*lo[1]*this->m_data[shift(a,1)  ] +
        lo[0]*hi[1]*this->m_data[shift(a,0,1)] +
        hi[0]*hi[1]*this->m_data[shift(a,1,1)];
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
inline typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
trilinear(const vec_type& u, size_type a) const
{
    const typename base_type::shifter_type& shift = this->grid().shifter();
    const vec_type& hi = u;
    vec_type lo = vec_type(1.) - hi;

    return
        lo[0]*lo[1]*lo[2]*this->m_data[a             ] +
        hi[0]*lo[1]*lo[2]*this->m_data[shift(a,1)    ] +
        lo[0]*hi[1]*lo[2]*this->m_data[shift(a,0,1)  ] +
        hi[0]*hi[1]*lo[2]*this->m_data[shift(a,1,1)  ] +
        lo[0]*lo[1]*hi[2]*this->m_data[shift(a,0,0,1)] +
        hi[0]*lo[1]*hi[2]*this->m_data[shift(a,1,0,1)] +
        lo[0]*hi[1]*hi[2]*this->m_data[shift(a,0,1,1)] +
        hi[0]*hi[1]*hi[2]*this->m_data[shift(a,1,1,1)];
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
inline typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
quadrilinear(const vec_type& u, size_type a) const
{
    const typename base_type::shifter_type& shift = this->grid().shifter();
    const vec_type& hi = u;
    vec_type lo = vec_type(1.) - hi;
    return
        lo[0]*lo[1]*lo[2]*lo[3]*this->m_data[a] +
        hi[0]*lo[1]*lo[2]*lo[3]*this->m_data[shift(a,1)      ] +
        lo[0]*hi[1]*lo[2]*lo[3]*this->m_data[shift(a,0,1)    ] +
        hi[0]*hi[1]*lo[2]*lo[3]*this->m_data[shift(a,1,1)    ] +
        lo[0]*lo[1]*hi[2]*lo[3]*this->m_data[shift(a,0,0,1)  ] +
        hi[0]*lo[1]*hi[2]*lo[3]*this->m_data[shift(a,1,0,1)  ] +
        lo[0]*hi[1]*hi[2]*lo[3]*this->m_data[shift(a,0,1,1)  ] +
        hi[0]*hi[1]*hi[2]*lo[3]*this->m_data[shift(a,1,1,1)  ] +
        lo[0]*lo[1]*lo[2]*hi[3]*this->m_data[shift(a,0,0,0,1)] +
        hi[0]*lo[1]*lo[2]*hi[3]*this->m_data[shift(a,1,0,0,1)] +
        lo[0]*hi[1]*lo[2]*hi[3]*this->m_data[shift(a,0,1,0,1)] +
        hi[0]*hi[1]*lo[2]*hi[3]*this->m_data[shift(a,1,1,0,1)] +
        lo[0]*lo[1]*hi[2]*hi[3]*this->m_data[shift(a,0,0,1,1)] +
        hi[0]*lo[1]*hi[2]*hi[3]*this->m_data[shift(a,1,0,1,1)] +
        lo[0]*hi[1]*hi[2]*hi[3]*this->m_data[shift(a,0,1,1,1)] +
        hi[0]*hi[1]*hi[2]*hi[3]*this->m_data[shift(a,1,1,1,1)] ;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
multilinear(const vec_type& u, size_type a, size_t n) const
{
    switch (n) {
        case 0:
            throw std::runtime_error("invalid dimension in raster_data::multilinear()");
        case 1:
            return linear(u, a);
        case 2:
            return bilinear(u, a);
        case 3:
            return trilinear(u, a);
        case 4:
            return quadrilinear(u, a);
        default:
            scalar_type alpha = 1. - u[n-1];
            scalar_type beta = u[n-1];
            value_type v = alpha * multilinear(u, a, n-1);
            v += beta * multilinear(u, this->grid().shifter().lift(a, n-1), n-1);
            return v;
    }
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
interpolate(const vec_type& x, const coord_type& c) const
{
    return multilinear(x, this->grid().index(c));
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
dlinear_part(const vec_type& p, size_type a) const
{
    const size_type& lo = a;
    /*
        linear:
            L(u) = (1-u)f(u=0) + u*f(u=1)
            dL/du(u) = f(u=1) - f(u=0)
    */
    return this->m_data[lo+1] - this->m_data[lo];
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
dbilinear_part(const vec_type& p, size_type a, size_type d) const
{
    /*
        bilinear:
            B(u,v) = (1-v)*L(u,v=0) + v*L(u,v=1)
            dB/du(u,v) = (1-v)*dL/du(u,v=0) + v*dL/du(u,v=1)
            dB/dv(u,v) = L(u,v=1) - L(u,v=0)
    */
    const scalar_type& v = p[1];
    const size_type& lo = a;
    const size_type hi = _lift(a, 1);
    switch(d) {
        case 0: return
            (1-v)*dlinear_part(p, lo) + v*dlinear_part(p, hi);
        case 1: return
            linear(p, hi) - linear(p, lo);
    }
    throw std::runtime_error("invalid dimension in      raster_data::dbilinear_part()");
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
dtrilinear_part(const vec_type& p, size_type a, size_type d) const
{
    /*
        trilinear:
            T(u,v,w) = (1-w)*B(u,v,w=0) + w*B(u,v,w=1)
            dT/du(u,v,w) = (1-w)*dB/du(u,v,w=0) + w*dB/dv(u,v,w=1)
            dT/dv(u,v,w) = (1-w)*dB/dv(u,v,w=0) + w*dB/dv(u,v,w=1)
            dT/dw(u,v,w) = B(u,v,w=1) - B(u,v,w=0)
    */
    const scalar_type& w = p[2];
    const size_type& lo = a;
    const size_type hi = _lift(a, 2);
    switch (d) {
        case 0: return
            (1-w)*dbilinear_part(p, lo, 0) + w*dbilinear_part(p, hi, 0);
        case 1: return
            (1-w)*dbilinear_part(p, lo, 1) + w*dbilinear_part(p, hi, 1);
        case 2: return
            bilinear(p, hi) - bilinear(p, lo);
    }
    throw std::runtime_error("invalid dimension in image<>::dtrilinear_part()");
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
dquadrilinear_part(const vec_type& p, size_type a, size_type d) const
{
    /*
        quadrilinear:
            Q(u,v,w,x) = (1-x)*T(u,v,w,x=0) + x*T(u,v,w,x=1)
            dQ/du(u,v,w,x) = (1-x)*dT/du(u,v,w,x=0) + x*dT/dv(u,v,w,x=1)
            dQ/dv(u,v,w,x) = (1-x)*dT/dv(u,v,w,x=0) + x*dT/dv(u,v,w,x=1)
            dQ/dw(u,v,w,x) = (1-x)*dT/dw(u,v,w,x=0) + x*dT/dw(u,v,w,x=1)
            dQ/dx(u,v,w,x) = T(u,v,w,x=1) - T(u,v,w,x=0)
    */
    const scalar_type& x = p[3];
    const size_type& lo = a;
    const size_type hi = _lift(a, 3);
    switch (d) {
        case 0: return
            (1-x)*dtrilinear_part(p, lo, 0) + x*dtrilinear_part(p, hi, 0);
        case 1: return
            (1-x)*dtrilinear_part(p, lo, 1) + x*dtrilinear_part(p, hi, 1);
        case 2: return
            (1-x)*dtrilinear_part(p, lo, 2) + x*dtrilinear_part(p, hi, 2);
        case 3: return
            trilinear(p, hi) - trilinear(p, lo);
    }
    throw std::runtime_error("invalid dimension in image<>::dquadrilinear_part()");
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
dmultilinear_part(const vec_type& p, size_type a, size_type d,
    size_type n) const
{
    switch (n) {
        case 0:
            throw std::runtime_error("invalid dimension in " \
                "raster_data::dmultilinear_part()");
        case 1:
            return dlinear_part(p, a);
        case 2:
            return dbilinear_part(p, a, d);
        case 3:
            return dtrilinear_part(p, a, d);
        case 4:
            return dquadrilinear_part(p, a, d);
        default:
            if (n == d+1)
                return multilinear(p, _lift(a, n-1), n-1) - multilinear(p, a, n-1);
            else {
                const scalar_type& u = p[n-1];
                return (1-u)*dmultilinear_part(p, a, d, n-1) +
                    u*dmultilinear_part(p, _lift(a, n-1), d, n-1);
            }
    }
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::deriv_type
image<Value_, Dim_, Scalar_, Size_>::
dmultilinear(const vec_type& p, size_type a, size_type n) const
{
    deriv_type r;
    for (size_t d=0 ; d<dim ; ++d) {
        switch(n) {
            case 0: std::runtime_error("invalid dimension in raster_data::dmultilinear");
            case 1: return deriv_type(dlinear_part(p, a));
            case 2: r[d] = dbilinear_part(p, a, d);
            case 3: r[d] = dtrilinear_part(p, a, d);
            case 4: r[d] = dquadrilinear_part(p, a, d);
            default: r[d] = dmultilinear_part(p, a, d, n);
        }
    }
    return r;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
central_diff(size_type i, size_type id, size_type dim) const
{
    if (i > 0 && i < this->m_grid.resolution()[dim] - 1) {
        return 0.5*(this->_data[id+this->m_shift[dim]] - this->m_data[id-this->m_shift[dim]]) /
               this->m_grid.step()[dim];
    } else if (i == 0) {
        return (this->m_data[id+this->m_shift[dim]] - this->m_data[id]) / this->m_grid.step()[dim];
    } else if (i == this->m_grid.resolution()[dim] - 1) {
        return (this->m_data[id] - this->m_data[id-this->m_shift[dim]]) / this->m_grid.step()[dim];
    } else {
        throw std::runtime_error("invalid coordinates in raster_data::central_diff()");
    }
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::value_type
image<Value_, Dim_, Scalar_, Size_>::
value(const point_type& p) const
{
    std::pair<coord_type, vec_type> r = this->m_grid.locate(p);
    return interpolate(r.second, r.first);
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::deriv_type
image<Value_, Dim_, Scalar_, Size_>::
derivative(const vec_type& p) const
{
    std::pair<coord_type, vec_type> r = this->m_grid.locate(p);
    size_type id = this->m_grid.index(r.first);
    deriv_type d = dmultilinear(r.second, id);
    for (size_type i=0 ; i<dim ; ++i) d[i] /= this->m_grid.spacing()[i];
    return d;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::val1
image<Value_, Dim_, Scalar_, Size_>::
derivative(size_type i) const
{
    static_assert(dim == 1, "Invalid dimension (!=1) in image::derivative(Int)");
    val1 r;
    r[0] = central_diff(i, i, 0);
    return r;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::val2
image<Value_, Dim_, Scalar_, Size_>::
derivative(size_type i, size_type j) const
{
    static_assert(dim == 2, "Invalid dimension (!=2) in image::derivative(Int, Int)");
    val2 r;
    coord_type id = this->m_grid.idx(i, j);
    r[0] = central_diff(i, id, 0);
    r[1] = central_diff(j, id, 1);
    return r;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::val3
image<Value_, Dim_, Scalar_, Size_>::
derivative(size_type i, size_type j, size_type k) const
{
    static_assert(dim == 3, "Invalid dimension (!=3) in image::derivative(Int, Int, Int)");
    val3 r;
    coord_type id = this->m_grid.idx(i, j, k);
    r[0] = central_diff(i, id, 0);
    r[1] = central_diff(j, id, 1);
    r[2] = central_diff(k, id, 2);
    return r;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::val4
image<Value_, Dim_, Scalar_, Size_>::
derivative(size_type i, size_type j, size_type k, size_type l)
const
{
    static_assert(dim == 4, "Invalid dimension (!=4) in image::derivative(Int, Int, Int, Int)");
    vec_type r;
    coord_type id = this->m_grid.idx(i, j, k, l);
    r[0] = central_diff(i, id, 0);
    r[1] = central_diff(j, id, 1);
    r[2] = central_diff(k, id, 2);
    r[3] = central_diff(k, id, 3);
    return r;
}

template<typename Value_, size_t Dim_, typename Scalar_, typename Size_>
typename image<Value_, Dim_, Scalar_, Size_>::deriv_type
image<Value_, Dim_, Scalar_, Size_>::
derivative(const coord_type& id) const
{
    vec_type r;
    for (size_type i = 0 ; i < dim ; ++i) {
        r[i] = central_diff(id[i], this->m_grid.idx(id), i);
    }
    return r;
}
} // namespace xavier
