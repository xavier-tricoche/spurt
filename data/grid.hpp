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
#include <string>


namespace spurt { namespace grid {
    
template<unsigned int N>
struct cell_type {
    typedef std::array<size_t, N> ids_type;
    static constexpr unsigned int nb_vertices = N;
    
    size_t index(unsigned int i) const {
        if (i >= N) throw std::runtime_error("Cell vertex index out of bounds: " + std::to_string(i));
        return m_ids[i];
    }
    
    size_t& index(unsigned int i) {
        if (i >= N) throw std::runtime_error("Cell vertex index out of bounds: " + std::to_string(i));
        return m_ids[i];
    }
    
    ids_type m_ids;
};

struct quadrilateral : public cell_type<4> {
};

struct tetrahedron : public cell_type<4> {
};

struct pyramid : public cell_type<5> {
};

struct prism : public cell_type<6> {
};

struct hexahedron : public cell_type<8> {
};

    
const nvis::ivec4 voxel_faces[6] = {
    nvis::ivec4(0, 1, 2, 3), nvis::ivec4(4, 5, 6, 7),
    nvis::ivec4(0, 1, 5, 4), nvis::ivec4(1, 2, 6, 5),
    nvis::ivec4(2, 3, 7, 6), nvis::ivec4(3, 0, 4, 7)
};

const nvis::ivec3 voxel_vertices[8] = {
    nvis::ivec3(0, 0, 0), nvis::ivec3(1, 0, 0),
    nvis::ivec3(1, 1, 0), nvis::ivec3(0, 1, 0),
    nvis::ivec3(0, 0, 1), nvis::ivec3(1, 0, 1),
    nvis::ivec3(1, 1, 1), nvis::ivec3(0, 1, 1)
};

/// Legacy wrapper for raster_grid
template<typename Type, int N>
class uniform_grid : public raster_grid<N, Type>
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
    typedef uniform_grid<Type, N>            self_type;
    typedef metric<Type, N>                  metric_type;
    typedef typename base_type::size_type    index_type;
    typedef typename base_type::coord_type   ivec_type;
    typedef typename base_type::vec_type     vec_type;
    typedef typename base_type::vec_type     vertex_type;
    typedef typename base_type::bounds_type  bounds_type;
    typedef nvis::fixed_vector<bool, N>      bvec_type;
    
    uniform_grid(const base_type& other) : base_type(other), 
        _periodic(false), _verbose(false) {
        _metric.bounds() = other.bounds();
        _metric.periodic() = _periodic;
    }

    uniform_grid(const ivec_type& size, const vec_type& spacing)
        : base_type(size, vec_type(0), spacing), 
          _periodic(false), _verbose(false) {
        _metric.bounds() = this->bounds();
        _metric.periodic() = _periodic;
    }
    
    uniform_grid(const ivec_type& size, const bounds_type& _bounds)
        : base_type(size, _bounds), _periodic(false), _verbose(false) {
        _metric.bounds() = this->bounds();
        _metric.periodic() = _periodic;
    }
    
    uniform_grid(const ivec_type& size, const vec_type& spacing, const bvec_type& periodic)
        : base_type(size, vec_type(0), spacing), _periodic(periodic),
          _verbose(false) {
        _metric.bounds() = this->bounds();
        _metric.periodic() = _periodic;
    }
    
    uniform_grid(const ivec_type& size, const vec_type& spacing, const vec_type& origin)
        : base_type(size, origin, spacing), _periodic(false),
          _verbose(false) {
        _metric.bounds() = this->bounds();
        _metric.periodic() = _periodic;
    }
    
    uniform_grid(const ivec_type& size, const vec_type& spacing, 
                 const vec_type& origin, const bvec_type& periodic)
        : base_type(size, origin, spacing), _periodic(periodic),
          _verbose(false) {
        _metric.bounds() = this->bounds();
        _metric.periodic() = _periodic;
    }
    
    void verbose(bool v) const {
        _verbose = v;
    }
    
    const ivec_type& dimensions() const { return this->resolution(); }

    std::pair<ivec_type, vec_type> local_coordinates(const vec_type& x) const {
        return this->locate(x);
    }
    
    ivec_type imodulo(const ivec_type& c) const {
        const ivec_type _size = this->resolution();
        if (!nvis::any(_periodic)) return c;
        ivec_type pc(c);
        for (int i = 0 ; i < N ; ++i) {
            if (_periodic[i]) {
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
        vec_type y = _metric.modulo(x);
        if (!this->bounds().inside(y))
            throw std::runtime_error("grid::dmodulo(): invalid coordinates");
        return y;
    }
    
    const metric_type& get_metric() const {
        return _metric;
    }
    
    bool on_boundary(const ivec_type& c) const {
        const ivec_type& _size = this->resolution();
        for (int i = 0 ; i < N ; ++i)
            if (!_periodic[i] && (c[i] == 0 || c[i] == _size[i] - 1))
                return true;
        return false;
    }
    
    bvec_type which_boundary(const ivec_type& c) const {
        const ivec_type& _size = this->resolution();
        bvec_type r(false);
        for (int i = 0 ; i < N ; ++i)
            if (!_periodic[i] && (c[i] == 0 || c[i] == _size[i] - 1))
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
    bvec_type         _periodic;
    metric_type       _metric;
    mutable    bool   _verbose;
};

} // namespace grid
} // namespace spurt


#endif






























