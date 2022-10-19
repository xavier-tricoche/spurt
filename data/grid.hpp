#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <boost/static_assert.hpp>
#include <stdexcept>
#include <data/metric.hpp>
#include <data/raster.hpp>


namespace spurt {
    
const ivec4 voxel_faces[6] = {
    ivec4(0, 1, 2, 3), ivec4(4, 5, 6, 7),
    ivec4(0, 1, 5, 4), ivec4(1, 2, 6, 5),
    ivec4(2, 3, 7, 6), ivec4(3, 0, 4, 7)
};

const ivec3 voxel_vertices[8] = {
    ivec3(0, 0, 0), ivec3(1, 0, 0),
    ivec3(1, 1, 0), ivec3(0, 1, 0),
    ivec3(0, 0, 1), ivec3(1, 0, 1),
    ivec3(1, 1, 1), ivec3(0, 1, 1)
};

/// Legacy wrapper for raster_grid
template<typename Type, int N>
class grid : public raster_grid<N, Type>
{
    static Type _modulo(Type a, Type b) {
        Type r = fmod(a, b - 1);
        if (!r) {
            return 0;
        } else {
            return a >= 0 ? r : b + r - 1;
        }
    }
    
public:
    typedef raster_grid<N, Type>             base_type;
    typedef grid<Type, N>                    self_type;
    typedef metric<Type, N>                  metric_type;
    typedef typename base_type::size_type    index_type;
    typedef typename base_type::coord_type   ivec_type;
    typedef typename base_type::vec_type     vec_type;
    typedef typename base_type::vec_type     vertex_type;
    typedef typename base_type::bounds_type  bounds_type;
    typedef fixed_vector<bool, N>      bvec_type;
    
    grid(const base_type& other) : base_type(other), 
        m_periodic(false), m_verbose(false) {
        m_metric.bounds() = other.bounds();
        m_metric.periodic() = m_periodic;
    }

    grid(const ivec_type& size, const vec_type& spacing)
        : base_type(size, vec_type(0), spacing), 
          m_periodic(false), m_verbose(false) {
        m_metric.bounds() = this->bounds();
        m_metric.periodic() = m_periodic;
    }
    
    grid(const ivec_type& size, const bounds_type& _bounds)
        : base_type(size, _bounds), m_periodic(false), m_verbose(false) {
        m_metric.bounds() = this->bounds();
        m_metric.periodic() = m_periodic;
    }
    
    grid(const ivec_type& size, const vec_type& spacing, const bvec_type& periodic)
        : base_type(size, vec_type(0), spacing), m_periodic(periodic),
          m_verbose(false) {
        m_metric.bounds() = this->bounds();
        m_metric.periodic() = m_periodic;
    }
    
    grid(const ivec_type& size, const vec_type& spacing, const vec_type& origin)
        : base_type(size, origin, spacing), m_periodic(false),
          m_verbose(false) {
        m_metric.bounds() = this->bounds();
        m_metric.periodic() = m_periodic;
    }
    
    grid(const ivec_type& size, const vec_type& spacing, 
         const vec_type& origin, const bvec_type& periodic)
        : base_type(size, origin, spacing), m_periodic(periodic),
          m_verbose(false) {
        m_metric.bounds() = this->bounds();
        m_metric.periodic() = m_periodic;
    }
    
    void verbose(bool v) const {
        m_verbose = v;
    }
    
    const ivec_type& dimensions() const { return this->resolution(); }

    std::pair<ivec_type, vec_type> local_coordinates(const vec_type& x) const {
        return this->locate(x);
    }
    
    ivec_type imodulo(const ivec_type& c) const {
        const ivec_type _size = this->resolution();
        if (!any(m_periodic)) return c;
        ivec_type pc(c);
        for (int i = 0 ; i < N ; ++i) {
            if (m_periodic[i]) {
                pc[i] = _modulo(c[i], _size[i]);
                // std::cerr << pc[i] << " = _modulo(" << c[i] << ", " << _size[i] << ")\n";
            }
            else if (c[i] < 0 || c[i] >= _size[i]) {
                // std::cerr << "c[" << i << "] = " << c[i] << " is out of range\n";
                throw std::runtime_error("grid::imodulo(): invalid coordinates");
            }
        }
        return pc;
    }
    
    vec_type dmodulo(const vec_type& x) const {
        vec_type y = m_metric.modulo(x);
        if (!this->bounds().inside(y))
            throw std::runtime_error("grid::dmodulo(): invalid coordinates");
        return y;
    }
    
    const metric_type& getm_metric() const {
        return m_metric;
    }
    
    bool on_boundary(const ivec_type& c) const {
        const ivec_type& _size = this->resolution();
        for (int i = 0 ; i < N ; ++i)
            if (!m_periodic[i] && (c[i] == 0 || c[i] == _size[i] - 1))
                return true;
        return false;
    }
    
    bvec_type which_boundary(const ivec_type& c) const {
        const ivec_type& _size = this->resolution();
        bvec_type r(false);
        for (int i = 0 ; i < N ; ++i)
            if (!m_periodic[i] && (c[i] == 0 || c[i] == _size[i] - 1))
                r[i] = true;
        return r;
    }
    
    std::list<index_type> neighbors(index_type idx) const {
        ivec_type c = this->coordinates(idx);
        std::list<index_type> r;
        for (int i = 0 ; i < N ; ++i) {
            index_type l = c[i];
            ivec_type previous(c), next(c);
            --previous[l];
            ++next[l];
            try {
                previous = modulo(previous);
                r.push_back(index(previous));
            } catch (...) {}
            try {
                next = modulo(next);
                r.push_back(index(next));
            } catch (...) {}
        }
        return r;
    }
    
private:
    bvec_type         m_periodic;
    metric_type       m_metric;
    mutable    bool   m_verbose;
};

}

#endif