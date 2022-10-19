#ifndef __NRRD_FIELD_HPP__
#define __NRRD_FIELD_HPP__

#include <iostream>
#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>
#include <data/raster.hpp>


namespace {
bool is_ok(double v)
{
    return !(std::isnan(v) || std::isinf(v));
}
}

namespace spurt {

namespace nrrd_utils {
template< typename T, int N >
void array2vector(std::vector<spurt::fixed_vector<double, N> >& out, const Nrrd* nrrd)
{
    int n = nrrd->axis[1].size * nrrd->axis[2].size * nrrd->axis[3].size;
    const T* data = (const T*)nrrd->data;
    out.resize(n);
    for (int i = 0 ; i < n ; ++i) {
        for (int j = 0 ; j < N ; ++j) {
            out[i][j] = data[N*i+j];
        }
    }
}

template< int N >
nvis::bounding_box<spurt::fixed_vector<double, N> > 
compute_bounds(const Nrrd* nrrd, bool verbose=false)
{
    typedef spurt::fixed_vector<double, N>     pos_type;
    typedef nvis::bounding_box<pos_type>        bbox_type;
    
    bbox_type bounds;
    pos_type lo, hi;
    for (int i = 0 ; i < N ; ++i) {
        const double& _min = nrrd->axis[i+1].min;
        const double& _max = nrrd->axis[i+1].max;
        const double& _spc = nrrd->axis[i+1].spacing;
        
        if (verbose) 
            std::cout << "axis " << i+1 
                      << ": min=" << _min 
                      << ", max=" << _max 
                      << ", step=" << _spc 
                      << std::endl;
        
        if (is_ok(_min) && is_ok(_max)) {
            lo[i] = _min;
            hi[i] = _max;
            if (verbose)
                std::cout << "min is ok = " << _min 
                          << " and max is ok = " << _max
                          << "spacing = " << _spc << '\n';
        }
        else if (is_ok(_spc)) {
            if (verbose) 
                std::cout << "spacing is ok: " << _spc << '\n';
            if (is_ok(_min)) {
                if (verbose) 
                    std::cout << "min is ok = " << _min << '\n';
                lo[i] = _min;
                hi[i] = _min + _spc * (nrrd->axis[i+1].size - 1);
            }
            else if (is_ok(_max)) {
                if (verbose)
                    std::cout << "max is ok = " << _max << '\n';
                hi[i] = _max;
                lo[i] = _max - _spc * (nrrd->axis[i+1].size - 1);
            }
            else {
                lo[i] = 0;
                hi[i] = _spc * (nrrd->axis[i+1].size - 1);
            }
        }
        else {
            if (is_ok(_min)) {
                if (verbose)
                    std::cout << "min is ok = " << _min << '\n';
                lo[i] = _min;
                hi[i] = _min + nrrd->axis[i+1].size - 1;
            }
            else if (is_ok(_max)) {
                if (verbose)
                    std::cout << "max is ok = " << _max << '\n';
                hi[i] = _max;
                lo[i] = _max - nrrd->axis[i+1].size - 1;
            }
            else {
                lo[i] = 0;
                hi[i] = nrrd->axis[i+1].size - 1;
            }
        }
        bounds.min()[i] = lo[i];
        bounds.max()[i] = hi[i];
    }

    return bounds;
}
}

template< int N >
class nrrd_field {
public:
    typedef fixed_vector<double, N>    value_type;
    typedef spurt::image<value_type, 3>     field_type;
    typedef typename field_type::grid_type   grid_type;
    
private:
    void initialize(const Nrrd* nrrd) {
        ivec3 res(nrrd->axis[1].size, nrrd->axis[2].size, nrrd->axis[3].size);
        _grid = new grid_type(res, nrrd_utils::compute_bounds<3>(nrrd));
        std::vector< value_type > data;
        switch (nrrd->type) {
        case nrrdTypeFloat:
            nrrd_utils::array2vector<float, N>(data, nrrd);
            break;
        case nrrdTypeDouble:
            nrrd_utils::array2vector<double, N>(data, nrrd);
            break;
        default:
            std::cerr << "nrrd_field::initialize(): " << nrrd->type << ": unsupported data type\n";
            throw;
        }
        _data = new field_type(*_grid, data);
    }

public:
    nrrd_field(const std::string& file_name) {
        Nrrd *nrrd = spurt::readNrrd(file_name);
        initialize(nrrd);
        nrrdNuke(nrrd);
    }

    nrrd_field(const Nrrd* nrrd) {
        initialize(nrrd);
    }

    ~nrrd_field() {
        delete _grid;
        delete _data;
    }

    value_type operator()(const vec3& x) const {
        return interpolate(x);
    }
    
    value_type operator()(const ivec3& c) const {
        return value(c[0], c[1], c[2]);
    }

    value_type value(int n) const {
        return (*_data)(_grid->coord(n));
    }

    value_type value(int u, int v, int w) const {
        return (*_data)(u, v, w);
    }

    value_type interpolate(const vec3& p) const {
        return _data->value(p);
    }

    const grid_type& grid() const {
        return *_grid;
    }

    const bbox3& bounds() const {
        return _grid->bounds();
    }

    const ivec3& resolution() const {
        return _grid->resolution();
    }

    const vec3& spacing() const {
        return _grid->spacing();
    }

private:
    grid_type    *_grid;
    field_type   *_data;
    bbox3  _bounds;
};

} // namespace spurt

#endif
































