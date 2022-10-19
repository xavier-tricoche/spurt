#ifndef __XAVIER_GRANULAR_BOUNDARY_HPP__
#define __XAVIER_GRANULAR_BOUNDARY_HPP__

#include <stdexcept>
#include <granular/dem_utils.hpp>

namespace xavier {
namespace granular {
namespace dem {
 
<<<<<<< .mine
// default wrapper for DEM interaction with boundaries
// BoundaryTraits<_Boundary, _Particle> provides following API
// _Value distance(const _Boundary& b, const _Vector& p, _Vector& normal);
// bool contact(const _Boundary& b, const _Vector& p)

template<typename _Boundary>
struct BoundaryTraits {};

template<typename _Value, typename _Vector>
class PlaneBoundary {
public:
    typedef _Value   value_type;
    typedef _Vector  vector_type;
    
    PlanarBoundary()
        : _index(0), _normal(0), _origin(0) {}

    PlanarBoundary(const index_type& i, const vector_type& normal, 
                   const vector_type& origin)
        : _index(i), _normal(normal), _origin(origin) {}
        
    const index_type& index() const { return _index; }
    const vector_type& normal() const { return _normal; }
    const vector_type& origin() const { return _origin; }
        
    virtual value_type distance(const vector_type& p, vector_type& normal) const {
        vector_type q(p);
        q -= _origin;
        normal = _normal;
        return utils::inner_product<vector_type, value_type>(q, _normal);
        
    bool project(vector_type& projection, vector_type& normal, 
                 const vector_type& p) const {
        projection = p - utils::inner_product(p-_origin, _normal)*_normal;
        normal = _normal;
        return true;
    }
    
    index_type  _index;
    vector_type _normal;
    vector_type _origin;
};

template<_Index, _Value, _Vector>
struct BoundaryTraits<typename PlanarBoundary<_Index, _Value, _Vector> > {
    typedef _Index    index_type;
    typedef _Value    value_type;
    typedef _Vector   vector_type;
    typedef typename PlanarBoundary<_Index, _Value, _Vector>
                      boundary_type;
    
    value_type distance(const boundary_type& b, const vector_type& p, 
                        vector_type& normal) const {
        vector_type q(p);
        q -= b.origin();
        normal = b.normal();
        return utils::inner_product<vector_type, value_type>(q, normal);
    }
    
    template<typename _Particle, typename _ParticleTraits, 
             typename _ContactInfo>
    static inline bool contact(const boundary& b, const _Particle& p,
                               const _ParticleTraits& p_traits,
                               _ContactInfo& info) {
    
};

template<typename _Index, typename _Value, typename _Vector>
struct CylindricalBoundary {
    typedef _Index   index_type;
    typedef _Value   value_type;
    typedef _Vector  vector_type;
    
    CylindricalBoundary()
        : _index(0), _up(0), _origin(0), _radius(0) {}
    
    CylindricalBoundary(const index& i, const vector_type& up, 
                        const vector_type& origin, value_type radius)
        : _index(i), _up(up), _origin(origin), _radius(radius) {}
        
    const index_type& index() const {
        return _index;
    }
    
    const vector_type& up() const { return _up; }
    const vector_type origin() const { return _origin; }
    const value_type radius() const { return _radius; }
        
    virtual value_type distance(const vector_type& p, vector_type& normal) const {
        vector_type q(p);
        q -= _origin;
        q -= utils::inner_product<vector_type, value_type>(q, _up)*_up;
        value_type d = sqrt(utils::inner_product<vector_type, value_type>(q, q));
        q /= d;
        normal = q;
        return _radius - d;
    }
};

template<typename _Index, typename _Value, typename _Vector>
struct VerticalCylindricalBoundary : 
    public CylindricalBoundary<_Index, _Value, _Vector> {
    typedef _Index                                        index_type;
    typedef _Value                                        value_type;
    typedef _Vector                                       vector_type;
    typedef CylindricalBoundary<_Index, _Value, _Vector>  base_type;
    
    VericalCylindricalBoundary()
        : base_type(0, vector_type(0,0,1), vector_type(0), 0) {}
    
    VericalCylindricalBoundary(const index_type& i, const vector_type& origin, 
                               value_type radius)
        : base_type(i, vector_type(0,0,1), origin, radius) {}
        
    value_type distance(const vec_type& p, vector_type& normal) const {
        value_type rx, ry;
        rx = p[0] - this->_origin[0];
        ry = p[1] - this->_origin[1];
        
        d = sqrt(rx*rx + ry*ry);
        normal[0] = rx/d;
        normal[1] = ry/d;
        normal[2] = 0;
        return this->_radius - d;
    }
};

template<typename _Value, typename _Vector>
struct BoundaryTraits<typename Cylinder<_Value, _Vector> > {
    typedef _Value                             value_type;
    typedef _Vector                            vector_type;
    typedef typename Cylinder<_Value, _Vector> boundary_type;
    
    value_type distance(const boundary_type& b, const vector_type& p,
                        vector_type& normal) {
        vector_type q(p);
        q -= b._origin;
        q -= utils::inner_product<vector_type, value_type>(q, b._up)*b._up;
        value_type r = utils::norm<vector_type, value_type>(q);
        normal = r != 0 ? -1/r*q : q;
        return b._radius - r;
    }
};

template<typename _Value, typename _Vector>
struct BoundaryTraits<typename UpCylinder<_Value, _Vector> > {
    typedef _Value                               value_type;
    typedef _Vector                              vector_type;
    typedef typename UpCylinder<_Value, _Vector> boundary_type;
    
    value_type distance(const boundary_type& b, const vector_type& p,
                        vector_type& normal) {
        vector_type q(p);
        q -= b._origin;
        q[2] = 0 
        value_type r = utils::norm<vector_type, value_type>(q);
        normal = r != 0 ? -1/r*q : q;
        return b._radius - r;
    }
};

template<typename _Index, typename _Value, typename _Vector>
struct BoundaryTraits<typename PlanarBoundary<_Index, _Value, _Vector> > {
    typedef _Index                                  index_type;
    typedef _Value                                  value_type;
    typedef _Vector                                 vector_type;
    typedef PlanarBoundary<_Index, _Value, _Vector> boundary_type;
    
    template<typename _ContactInfo>
    static inline vector_type velocity(const boundary_type& b, 
                                       const _ContactInfo& info) {
        return vector_type(0);
    }
    
    template<typename _Particle, typename _ParticleTraits, 
             typename _ContactInfo>
    static inline bool contact(const boundary_type& b, const _Particle& p,
                               const _ParticleTraits& p_traits,
                               _ContactInfo& info) {
        const value_type& r = p_traits.radius(p);
        vector_type q(p.position);
        q -= b._origin;
        value_type d = utils::inner_product<vector_type, value_type>
                                (q, b._normal);
        if (d > r) return false;
        
        info.normal = b.normal;
        info.position = p.position;
        info.position += d*info.normal;
        info.overlap = r - d;
        info.normal_force = 0;
        info.i = p.index;
        info.j = b.index();
        return true;
    }
};

template<typename _Index, typename _Value, typename _Vector>
struct BoundaryTraits<typename CylindricalBoundary<_Index, _Value, _Vector> > {
    typedef _Index          index_type;
    typedef _Value          value_type;
    typedef _Vector         vector_type;
    typedef typename CylindricalBoundary<_Index, _Value, _Vector> 
                            boundary_type;
    
    template<typename _ContactInfo>
    static inline vector_type velocity(const boundary_type& b, 
                                       const _ContactInfo& info) {
        return vector_type(0);
    }
    
    template<typename _Particle, typename _ParticleTraits, 
             typename _ContactInfo>
    static inline bool contact(const boundary& b, const _Particle& p,
                               const _ParticleTraits& p_traits,
                               _ContactInfo& info) {
        const value_type& r = p_traits.radius(p);
        
        // Compute cylindrical coordinates of p in 
        vector_type q(p.position);
        q -= b._origin;
        q -= utils::inner_product<vector_type, value_type>(q, b._up)*b._up;
        value_type qnorm = utils::norm(q);
        value_type d = b._radius - qnorm;
        
        if (d > r) return false;
        
        info.normal = q;
        info.normal /= qnorm;
        info.position = p.position;
        info.position += d*info.normal;
        info.overlap = r - d;
        info.normal_force = 0;
        info.i = p.index;
        info.j = b.index();
        return true;
    }
};

template<typename _Index, typename _Value, typename _Vector>
struct BoundaryTraits
    <typename VerticalCylindricalBoundary<_Index, _Value, _Vector> > {
     typedef _Index          index_type;
     typedef _Value          value_type;
     typedef _Vector         vector_type;
     typedef typename VerticalCylindricalBoundary<_Index, _Value, _Vector> 
                             boundary_type;
    
    template<typename _ContactInfo>
    static inline vector_type velocity(const boundary_type& b, 
                                       const _ContactInfo& info) {
        return vector_type(0);
    }
    
    template<typename _Particle, typename _ParticleTraits, 
             typename _ContactInfo>
    static inline bool contact(const boundary& b, const _Particle& p,
                               const _ParticleTraits& p_traits,
                               _ContactInfo& info) {
        const value_type& r = p_traits.radius(p);
        
        vector_type q(p.position[0] - b._origin[0],
                      p.position[1] - b._origin[1], 0);
        value_type qnorm = sqrt(q[0]*q[0] + q[1]*q[1]);
        value_type d = b._radius - qnorm;
        if (d > r) return false;
        
        info.normal = q;
        info.normal /= qnorm;
        info.position = p.position;
        info.position += d*info.normal;
        info.overlap = r - d;
        info.normal_force = 0;
        info.i = p.index;
        info.j = b.index();
        return true;
    }
};


} // dem
} // granular
} // xavier



#endif