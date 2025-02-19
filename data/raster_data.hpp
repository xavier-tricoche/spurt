#pragma once

#include <assert.h>
// STL
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <memory>
#include <fstream>
#include <math/types.hpp>
#include <data/raster_grid.hpp>

namespace spurt
{

template <typename Size_, typename Scalar_, size_t Dim, typename Value_,
          typename Coord_ = small_vector<Size_, Dim>, 
          typename Pos_ = small_vector<Scalar_, Dim> >
class raster_data
{
public:
    typedef raster_grid<Size_, Scalar_, Dim, Coord_, Pos_> grid_type;
    typedef typename grid_type::size_type size_type;
    typedef typename grid_type::coord_type coord_type;
    typedef typename grid_type::pos_type pos_type;
    typedef typename grid_type::scalar_type scalar_type;
    typedef typename grid_type::bounds_type bounds_type;
    typedef Value_ value_type;
    typedef value_type& reference_type;
    typedef const value_type& const_reference_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    static const size_type dimension = grid_type::dimension;
    typedef raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_> self_type;

    // constructors
    raster_data() {}

    raster_data(const self_type &other) 
        : m_grid(other.m_grid), m_data(other.m_data) {}

    raster_data(const grid_type &grid, bool allocate=true) 
        : m_grid(grid), m_data() {
        if (allocate) {
            m_data.reset(new value_type[grid.size()]);
        }
    }

    raster_data(const grid_type& grid, value_type* data)
        : m_grid(grid), m_data(data) {}

    raster_data(const grid_type& grid, std::shared_ptr<value_type> data)
        : m_grid(grid), m_data(data) {}

    raster_data(const grid_type &grid, value_type init_val) 
        : m_grid(grid), m_data() {
        m_data.reset(new value_type[grid.size()]);
        std::fill(begin(), end(), init_val);
    } 

    raster_data(const grid_type &grid, const std::vector<value_type> &data) 
        : m_grid(grid), m_data() {
        assert(data.size() == grid.size());
        m_data.reset(new value_type[grid.size()]);
        std::copy(data.begin(), data.end(), begin());
    }

    template <typename _Iterator>
    raster_data(const grid_type &grid, _Iterator _begin, _Iterator _end) 
        : m_grid(grid), m_data() { 
        assert(std::distance(_begin, _end) == grid.size());
        m_data.reset(new value_type[grid.size()]);
        std::copy(_begin, _end, begin());
    }

    ~raster_data() {}

    void set_data(value_type* ptr) {
        m_data.reset(ptr);
    }

    void set_data(std::shared_ptr<value_type> sptr) {
        m_data = sptr;
    }

    void initialize(const value_type& v) {
         if (m_data.get() == nullptr) {
            m_data.reset(new value_type[m_grid.size()]);
         }
         std::fill(begin(), end(), v);
    }

    const grid_type &grid() const { return m_grid; }

    iterator begin() { return m_data.get(); }
    const_iterator begin() const { 
        return const_cast<const_iterator>(m_data.get()); 
    }
    iterator end() { 
        return m_data.get() + m_grid.size(); 
    }
    const_iterator end() const { 
        return const_cast<const_iterator>(m_data.get() + m_grid.size()); 
    }

    const size_type &size() const { return m_grid.size(); }

    // lookup
    const_reference_type operator()(size_type i) const
    {
        return m_data.get()[m_grid.index(i)];
    }
    const_reference_type operator()(size_type i, size_type j) const
    {
        return m_data.get()[m_grid.index(i, j)];
    }
    const_reference_type operator()(size_type i, size_type j, size_type k) const
    {
        return m_data.get()[m_grid.index(i, j, k)];
    }
    const_reference_type operator()(size_type i, size_type j, size_type k,
                                    size_type l) const
    {
        return m_data.get()[m_grid.index(i, j, k, l)];
    }
    const_reference_type operator()(const coord_type &id) const
    {
        return m_data.get()[m_grid.index(id)];
    }
    // array-type access for convenience
    const_reference_type operator[](size_type idx) const
    {
        return m_data.get()[idx];
    }

    reference_type operator()(size_type i)
    {
        return m_data.get()[m_grid.index(i)];
    }
    reference_type operator()(size_type i, size_type j)
    {
        return m_data.get()[m_grid.index(i, j)];
    }
    reference_type operator()(size_type i, size_type j, size_type k)
    {
        return m_data.get()[m_grid.index(i, j, k)];
    }
    reference_type operator()(size_type i, size_type j, size_type k,
                              size_type l)
    {
        return m_data.get()[m_grid.index(i, j, k, l)];
    }
    reference_type operator()(const coord_type &id)
    {
        return m_data.get()[m_grid.index(id)];
    }
    reference_type operator[](size_type idx)
    {
        return m_data.get()[idx];
    }
    
    value_type* data() { return m_data.get(); }
    const value_type* data() const { 
        return const_cast<const value_type*>(m_data.get()); 
    }

protected:
    grid_type m_grid;
    std::shared_ptr<value_type> m_data;
};

template <typename _Value>
using raster1d = raster_data<long, double, 1, _Value>;
template <typename _Value>
using raster2d = raster_data<long, double, 2, _Value>;
template <typename _Value>
using raster3d = raster_data<long, double, 3, _Value>;
template <typename _Value>
using raster4d = raster_data<long, double, 4, _Value>;

template <typename Size_, typename Scalar_, size_t Dim, typename Value_, typename Coord_, typename Pos_, typename AltValue_ = Value_>
void save_as_nrrd(const std::string &filename, 
                  const raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_>& data, const AltValue_& notused=AltValue_(0))
{
    typedef raster_data<Size_, Scalar_, Dim, Value_, Coord_, Pos_> data_type;
    typedef data_traits<AltValue_> dtraits;
    typedef AltValue_ value_type;
    typedef typename dtraits::value_type scalar_type;
    std::string type_string = spurt::type2string<scalar_type>::type_name();
    size_t ncomp = dtraits::size();

    std::ostringstream os;
    os << "NRRD0001\n";
    os << "# Complete NRRD file format specification at:\n";
    os << "# http://teem.sourceforge.net/nrrd/format.html\n";
    os << "type: " << type_string << '\n';
    os << "dimension: " << (ncomp > 1 ? Dim + 1 : Dim) << '\n';
    os << "sizes:";
    if (ncomp > 1)
        os << " " << ncomp;
    for (int i = 0; i < Dim; ++i)
    {
        os << " " << data.grid().resolution()[i];
    }
    os << '\n';
    os << "spacings:";
    if (ncomp > 1)
        os << " nan";
    for (int i = 0; i < Dim; ++i)
    {
        os << " " << std::setprecision(17) 
           << data.grid().spacing()[i];
    }
    os << '\n';
    os << "axis mins:";
    if (ncomp > 1)
        os << " nan";
    for (int i = 0; i < Dim; ++i)
    {
        os << " " << std::setprecision(17) 
           << data.grid().bounds().min()[i];
    }
    os << '\n';
    os << "centerings:";
    if (ncomp > 1)
        os << " ???";
    for (int i = 0; i < Dim; ++i)
    {
        os << " node";
    }
    os << '\n';
    os << "endian: little\n";
    os << "encoding: raw\n";

    size_t sz = ncomp * data.grid().size() * sizeof(scalar_type);
    std::cout << '\n' << filename << '\n';
    std::cout << "exporting "  << sz << " bytes" << '\n';
    std::cout << "grid res: "  << data.grid().resolution() << '\n';
    std::cout << "grid size: " << data.grid().size() << '\n';
    std::cout << "data type: " << type_string << '\n';
 
    std::ofstream out(filename.c_str(), std::ios::binary);
    out << os.str() << std::endl;
    out.write((char *)&data[0], sz);
    out.close();
}

} // namespace spurt