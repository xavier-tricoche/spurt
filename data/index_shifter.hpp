#pragma once

// STL
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <memory>
#include <math/types.hpp>

namespace spurt
{
    template<typename coord_type, typename sizes_type>
    long coord_to_index(const coord_type &coord, const sizes_type &sizes);

    template<typename coord_type, typename sizes_type>
    coord_type index_to_coord(long idx, const sizes_type &size);

    template <typename Size_, size_t Dim>
    class index_shifter
    {
    public:
        static const size_t dimension = Dim;
        typedef Size_ size_type;
        typedef small_vector<size_type, Dim> array_type;
        typedef small_vector<size_type, Dim> coord_type;

        template<typename T=size_type>
        index_shifter(const small_vector<T, Dim>& sizes=small_vector<T, Dim>(0))
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
        
        index_shifter(const index_shifter &other) 
            : m_sizes(other.m_sizes), m_dim_shifts(other.m_dim_shifts), 
              m_cell_shifts(other.m_cell_shifts) {}

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
        size_type operator()(size_type base, const coord_type &c) const;
        size_type index(size_type base, size_type local_id) const;

        // per-cell index shift of all the vertices
        const std::vector<size_t> &cell_shifts(const size_type &base) const
        {
            return m_cell_shifts;
        }

        // map local index to local coordinates
        coord_type local_coords(size_t local_id) const;

    private:
        array_type m_sizes;                   // raster dimensions
        array_type m_dim_shifts;              // per-dimension index shift
        std::vector<size_type> m_cell_shifts; // inner-cell index shifts

        friend std::ostream &operator<<(std::ostream &os, 
                                        const index_shifter &is)
        {
            os << "[ [sizes: " << is.m_sizes
               << "]\n[dim shifts: "
               << is.m_dim_shifts
               << "]\n[cell shifts: ";
            std::copy(is.m_cell_shifts.begin(), is.m_cell_shifts.end(),
                      std::ostream_iterator<size_type>(os, ", "));
            os << "] ]";
            return os;
        }
    };

    template <typename coord_type, typename sizes_type>
    long coord_to_index(const coord_type &coord, const sizes_type &sizes)
    {
        size_t N = coord.size();

        // i + size0*(j + size1*(k + size2*(l + ....)))
        auto idx = coord[N - 1];
        for (long dim = N - 2; dim >= 0; --dim)
        {
            idx = coord[dim] + idx * sizes[dim];
        }
        return idx;
    }

    template <typename coord_type, typename sizes_type>
    coord_type index_to_coord(long idx, const sizes_type &sizes)
    {
        size_t N = sizes.size();

        coord_type coord;
        typedef typename coord_type::value_type value_t;

        // i + size0*( j + size1*( k + size2*( l ) ) )
        for (size_t i = 0; i < N; ++i)
        {
            std::ldiv_t qr = std::div(idx, sizes[i]);
            coord[i] = static_cast<value_t>(qr.rem);
            idx = qr.quot;
        }
        return coord;
    }

    /***************************************************************************

                                    index_shifter

    ***************************************************************************/

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::coord_type
    index_shifter<Size_, N>::local_coords(size_t local_id) const
    {
        coord_type r;
        std::bitset<sizeof(size_type)> bits(local_id);
        for (size_t i = 0; i < r.size(); ++i)
        {
            r[i] = bits[i];
        }
        return r;
    }

    // template <typename Array>
    // index_shifter<Array>::index_shifter(const array_type &sizes)
    //     : m_sizes(sizes)
    // {
    //     std::cout << "shifter m_sizes=" << m_sizes << '\n';
    //     size_type shift = 1;
    //     for (size_t i = 0; i < dimension; ++i)
    //     {
    //         m_dim_shifts[i] = shift;
    //         shift *= m_sizes[i];
    //     }

    //     size_t points_per_cell = 1;
    //     points_per_cell = points_per_cell << dimension;
    //     std::cout << "shifter: points_per_cell=" << points_per_cell << '\n';
    //     m_cell_shifts.resize(points_per_cell);
    //     std::fill(m_cell_shifts.begin(), m_cell_shifts.end(), 0);

    //     for (size_t i = 1; i < points_per_cell; ++i)
    //     {
    //         std::bitset<sizeof(size_type)> bits(i);
    //         for (size_t j = 0; j < dimension; ++j)
    //         {
    //             if (bits[j])
    //                 m_cell_shifts[i] += m_dim_shifts[i];
    //         }
    //     }
    // }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::
    operator()(size_type base, size_type i) const
    {
        return base + i;
    }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::
    operator()(size_type base, size_type i, size_type j) const
    {
        return base + i + j * m_dim_shifts[1];
    }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::
    operator()(size_type base, size_type i, size_type j, size_type k) const
    {
        return base + i + j * m_dim_shifts[1] + k * m_dim_shifts[2];
    }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::
    operator()(size_type base, size_type i, size_type j, size_type k,
               size_type l) const
    {
        return base + i + j * m_dim_shifts[1] + k * m_dim_shifts[2] + l * m_dim_shifts[3];
    }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::
    operator()(size_type base, const coord_type &c) const
    {
        size_type r(base);
        r += c[0];
        for (size_t i = 1; i < dimension; ++i)
            r += c[i] * m_dim_shifts[i];
        return r;
    }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::lift(size_type base, size_type d) const
    {
        return base + m_dim_shifts[d];
    }

    template <typename Size_, size_t N>
    inline typename index_shifter<Size_, N>::size_type
    index_shifter<Size_, N>::index(size_type base, size_type local_id) const
    {
        return base + m_cell_shifts[local_id];
    }

} // namespace spurt
