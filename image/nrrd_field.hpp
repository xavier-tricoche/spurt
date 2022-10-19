#ifndef __NRRD_FIELD_HPP__
#define __NRRD_FIELD_HPP__

#include <iostream>
#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>
#include <data/raster.hpp>
#include <misc/meta_utils.hpp>


namespace xavier {

template< typename Value_, int N, typename Scalar_=double >
class nrrd_field {
public:
    static const size_t dim = N;
    static const size_t valsize = data_traits<Value_>::size();
    typedef Scalar_                          scalar_type; // a number
    typedef Value_                           value_type;  // an attribute
    typedef xavier::image<value_type, N>     field_type;  // a continuous function
    typedef typename field_type::grid_type   grid_type;   // a spatial data structure
    typedef typename field_type::point_type  point_type;  // a position in space
    typedef typename field_type::coord_type  coord_type;  // discrete coordinates
    typedef typename field_type::bounds_type bounds_type; // physical bounds
    typedef typename field_type::vec_type    vector_type; // a vector 
    
private:
    void initialize(const Nrrd* nrrd) {
        coord_type res;
        
        size_t offset = (valsize == 1 ? 0 : 1);
        size_t required_dim = dim + offset;
        
        if (nrrd->dim < required_dim) {
            throw std::runtime_error("Invalid field dimension exceeds Nrrd dimension");
        }
        
        for (int i=0; i<dim; ++i) res[i] = nrrd->axis[i+offset].size;
        
        m_grid = new grid_type(res, nrrd_utils::get_bounds<dim>(nrrd, offset));
        std::vector< value_type > data;
        nrrd_utils::to_vector(data, nrrd);
        m_data = new field_type(*m_grid, data);
    }

public:
    nrrd_field(const std::string& file_name, bool is_scalar=true) {
        Nrrd *nrrd = nrrd_utils::readNrrd(file_name);
        initialize(nrrd);
        nrrdNuke(nrrd);
    }

    nrrd_field(const Nrrd* nrrd) {
        initialize(nrrd);
    }

    ~nrrd_field() {
        delete m_grid;
        delete m_data;
    }

    value_type operator()(const point_type& x) const {
        return m_data->value(x);
    }

    value_type value(int n) const {
        return (*m_data)(m_grid->coord(n));
    }

    value_type value(int u, int v, int w=0) const {
        switch (dim) {
            case 2:return (*m_data)(u, v);
            case 3: return (*m_data)(u, v, w);
            default: throw std::runtime_error("Unsupported grid dimension");
        }
    }
    
    value_type value(const coord_type& c) const {
        return field_type::value(c);
    }

    const grid_type& grid() const {
        return *m_grid;
    }

    const bounds_type& bounds() const {
        return m_grid->bounds();
    }

    const coord_type& resolution() const {
        return m_grid->resolution();
    }

    const coord_type& spacing() const {
        return m_grid->spacing();
    }

private:
    grid_type    *m_grid;
    field_type   *m_data;
    bounds_type  m_bounds;
};

} // namespace xavier

#endif
































