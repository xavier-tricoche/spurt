#ifndef __raster_data_hpp__
#define __raster_data_hpp__

#include "grid.hpp"
#include <vector>
#include <assert.h>
#include <stdexcept>
#include <boost/static_assert.hpp>

namespace xavier
{

/// legacy wrapper for raster_data<>
template<typename T1, typename T2, int N>
class raster_data
{
public:
    typedef T1                                  data_type;
    typedef nvis::fixed_vector<data_type, N>    derivative_type;
    typedef T2                                  scalar_type;
    typedef grid<scalar_type, N>                grid_type;
    typedef typename grid_type::vec_type        vec_type;
    typedef typename grid_type::ivec_type       ivec_type;

    struct invalid_point : public std::runtime_error {
        invalid_point() : std::runtime_error("invalid point") {
        }
    };

    raster_data(const grid_type& grid, const data_type& init_value)
        : __data(grid.size(), init_value), __grid(grid), __verbose(false) {
    }

    raster_data(const grid_type& grid, const std::vector<data_type>& data)
        : __data(data), __grid(grid), __verbose(false) {
        assert(__data.size() == __grid.size());
    }

    template<typename iterator_type>
    raster_data(const grid_type& grid,
                const iterator_type& data_begin, const iterator_type& data_end)
        : __data(data_begin, data_end), __grid(grid), __verbose(false) {
        assert(__data.size() == __grid.size());
    }

    void verbose(bool v) const {
        __verbose = v;
    }

    const data_type& operator()(const ivec_type& c) const {
        // this will throw an exception if the coordinates are invalid
        return __data[__grid.index(c)];
    }

    data_type& operator()(const ivec_type& c) {
        // this will throw an exception if the coordinates are invalid
        return __data[__grid.index(c)];
    }

    data_type interpolate(const vec_type& x) const {
        // apply bilinear / trilinear interpolation
        BOOST_STATIC_ASSERT(N == 2 || N == 3);

        // an exception will be thrown if the position is outside the grid
        if (__verbose) __grid.verbose(true);
        std::pair<ivec_type, vec_type> tmp = __grid.local_coordinates(x);
        if (__verbose) __grid.verbose(false);
        ivec_type id = tmp.first;
        vec_type z = tmp.second;
        vec_type Z = vec_type(1) - z;

        if (N == 2) {
            // bilinear case
            return Z[0]*Z[1]*(*this)(id) +
                   z[0]*Z[1]*(*this)(id + ivec_type(1, 0)) +
                   z[0]*z[1]*(*this)(id + ivec_type(1, 1)) +
                   Z[0]*z[1]*(*this)(id + ivec_type(0, 1));
        }
        else {
            // trilinear case
            return Z[0]*Z[1]*Z[2]*(*this)(id) +
                   z[0]*Z[1]*Z[2]*(*this)(id + ivec_type(1, 0, 0)) +
                   z[0]*z[1]*Z[2]*(*this)(id + ivec_type(1, 1, 0)) +
                   Z[0]*z[1]*Z[2]*(*this)(id + ivec_type(0, 1, 0)) +
                   Z[0]*Z[1]*z[2]*(*this)(id + ivec_type(0, 0, 1)) +
                   z[0]*Z[1]*z[2]*(*this)(id + ivec_type(1, 0, 1)) +
                   z[0]*z[1]*z[2]*(*this)(id + ivec_type(1, 1, 1)) +
                   Z[0]*z[1]*z[2]*(*this)(id + ivec_type(0, 1, 1));
        }
    }

    derivative_type derivative(const vec_type& x) const {
        // apply bilinear / trilinear interpolation
        BOOST_STATIC_ASSERT(N == 2 || N == 3);

        // an exception will be thrown if the position is outside the grid
        if (__verbose) __grid.verbose(true);
        std::pair<ivec_type, vec_type> tmp = __grid.local_coordinates(x);
        if (__verbose) __grid.verbose(false);
        ivec_type id = tmp.first;
        vec_type z = tmp.second;
        vec_type Z = vec_type(1) - z;
        vec_type dz = 1. / __grid.spacing();
        vec_type dZ = -1*dz;

        derivative_type df;
        if (N == 2) {
            // bilinear case
            df[0] = dZ[0] * Z[1] * (*this)(id) +
                    dz[0] * Z[1] * (*this)(id + ivec_type(1, 0)) +
                    dz[0] * z[1] * (*this)(id + ivec_type(1, 1)) +
                    dZ[0] * z[1] * (*this)(id + ivec_type(0, 1));
            df[1] = Z[0] * dZ[1] * (*this)(id) +
                    z[0] * dZ[1] * (*this)(id + ivec_type(1, 0)) +
                    z[0] * dz[1] * (*this)(id + ivec_type(1, 1)) +
                    Z[0] * dz[1] * (*this)(id + ivec_type(0, 1));
        }
        else {
            // trilinear case
            df[0] = dZ[0] * Z[1] * Z[2] * (*this)(id) +
                    dz[0] * Z[1] * Z[2] * (*this)(id + ivec_type(1, 0, 0)) +
                    dz[0] * z[1] * Z[2] * (*this)(id + ivec_type(1, 1, 0)) +
                    dZ[0] * z[1] * Z[2] * (*this)(id + ivec_type(0, 1, 0)) +
                    dZ[0] * Z[1] * z[2] * (*this)(id + ivec_type(0, 0, 1)) +
                    dz[0] * Z[1] * z[2] * (*this)(id + ivec_type(1, 0, 1)) +
                    dz[0] * z[1] * z[2] * (*this)(id + ivec_type(1, 1, 1)) +
                    dZ[0] * z[1] * z[2] * (*this)(id + ivec_type(0, 1, 1));
            df[1] = Z[0] * dZ[1] * Z[2] * (*this)(id) +
                    z[0] * dZ[1] * Z[2] * (*this)(id + ivec_type(1, 0, 0)) +
                    z[0] * dz[1] * Z[2] * (*this)(id + ivec_type(1, 1, 0)) +
                    Z[0] * dz[1] * Z[2] * (*this)(id + ivec_type(0, 1, 0)) +
                    Z[0] * dZ[1] * z[2] * (*this)(id + ivec_type(0, 0, 1)) +
                    z[0] * dZ[1] * z[2] * (*this)(id + ivec_type(1, 0, 1)) +
                    z[0] * dz[1] * z[2] * (*this)(id + ivec_type(1, 1, 1)) +
                    Z[0] * dz[1] * z[2] * (*this)(id + ivec_type(0, 1, 1));
            df[2] = Z[0] * Z[1] * dZ[2] * (*this)(id) +
                    z[0] * Z[1] * dZ[2] * (*this)(id + ivec_type(1, 0, 0)) +
                    z[0] * z[1] * dZ[2] * (*this)(id + ivec_type(1, 1, 0)) +
                    Z[0] * z[1] * dZ[2] * (*this)(id + ivec_type(0, 1, 0)) +
                    Z[0] * Z[1] * dz[2] * (*this)(id + ivec_type(0, 0, 1)) +
                    z[0] * Z[1] * dz[2] * (*this)(id + ivec_type(1, 0, 1)) +
                    z[0] * z[1] * dz[2] * (*this)(id + ivec_type(1, 1, 1)) +
                    Z[0] * z[1] * dz[2] * (*this)(id + ivec_type(0, 1, 1));
        }
        return df;
    }

    const grid_type& get_grid() const {
        return __grid;
    }

    const std::vector<data_type>& get_data() const {
        return __data;
    }

    std::vector<data_type>& get_data() {
        return __data;
    }

    size_t size() const {
        return get_data().size();
    }

private:
    std::vector<data_type>  __data;
    const grid_type&        __grid;
    mutable bool            __verbose;
};

} // namespace xavier

#endif






















