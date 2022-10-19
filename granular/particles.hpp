#ifndef __XAVIER_GRANULAR_PARTICLES_HPP__
#define __XAVIER_GRANULAR_PARTICLES_HPP__


namespace xavier {
namespace granular {
namespace dem {

template<typename _Particle>
class ParticleTraits {};

template<typename _Index, typename _Value, typename _Vector>
struct Spherical {
    typedef _Index  index_type;
    typedef _Value  value_type;
    typedef _Vector vector_type;
    
    Spherical() 
        : _radius(0), _mass(0), _k_load(0), _k_unload(0),
          _k0(0), _gamma(0), _mu(0), 
          index(0), position(0), velocity(0), angular_velocity(0), 
          force(0), torque(0), displacement(0) {}
          
    Spherical(const index& _index, const value_type& _radius,
              const value_type& _mass, const value_type& _k_load,
              const value_type& _k_unload, const value_type& _k_ratio,
              const value_type& _gamma, const value_type& mu,
              const vector_type& _position, 
              const vector_type& _velocity = vector_type(0),
              const vector_type& _angular_velocity = vector_type(0),
              const vector_type& _force = vector_type(0),
              const vector_type& _torque = vector_type(0)) 
        : _radius(_radius), _mass(_mass), _k_load(_k_load), 
          _k_unload(_k_unload), _k_ratio(_k_ratio), _gamma(_gamma), 
          _mu(_mu),
          index(_index), position(_position), velocity(_velocity),
          angular_velocity(_angular_velocity), force(_force),
          torque(_torque), displacement(0) {}
    
    // Specific physical and geometric attributes of the particle
    value_type _radius, _mass, _k_load, _k_unload, 
               _k_ratio, _gamma, _mu;
    
    // These values describe the identity and state of a particle
    index_type  index;
    vector_type position;
    vector_type velocity;
    vector_type angular_velocity;
    vector_type force;
    vector_type torque;
    vector_type displacement;
};

template<typename _Index, typename _Value, typename _Vector>
struct GenericSpherical {
    typedef _Index  index_type;
    typedef _Value  value_type;
    typedef _Vector vector_type;
    
    GenericSpherical() 
        : index(0), position(0), velocity(0), angular_velocity(0), 
          force(0), torque(0), displacement(0) {}
              
    GenericSpherical(const index_type& _index, 
                     const vector_type& _position, 
                     const vector_type& _velocity = vector_type(0),
                     const vector_type& _angular_velocity = vector_type(0),
                     const vector_type& _force = vector_type(0),
                     const vector_type& _torque = vector_type(0)) 
        : index(_index), position(_position), velocity(_velocity),
          angular_velocity(_angular_velocity), force(_force),
          torque(_torque), displacement(0) {}

    // These values describe the specific state of any given particle
    index_type  index;
    vector_type position;
    vector_type velocity;
    vector_type angular_velocity;
    vector_type force;
    vector_type torque;
    vector_type displacement;
};

template<typename _Index, typename _Value, typename _Vector>
class ParticleTraits<typename Spherical<_Index, _Value, _Vector> > {
public:
    typedef _Index                                       index_type;
    typedef _Value                                       value_type;
    typedef _Vector                                      vector_type;
    typedef typename Spherical<_Index, _Value, _Vector>  particle;
    
    static const value_type& radius(const particle& p) {
        return p._radius;
    }
    static value_type rcut(const particle& p1, const particle& p2) {
        return p1._radius + p2._radius;
    }
    static value_type rcut_square(const particle& p1, const particle& p2) {
        value_type rc = rcut(p1, p2);
        return rc*rc;
    }
    static const value_type& mass(const particle& p) {
        return p._mass;
    }
    static const value_type& k_load(const particle& p) {
        return p._k_load;
    }
    static const value_type& k_unload(const particle& p) {
        return p._k_unload;
    }
    static const value_type& k_ratio(const particle& p) {
        return p._k_ratio;
    }
    static const value_type& gamma(const particle& p) {
        return p._gamma;
    }
    static const value_type& mu(const particle& p) {
        return p._mu;
    }
    static value_type inertia(const particle& p) {
        const value_type& m = p._mass;
        const value_type& r = p._radius;
        return m*r*r;
    }
};

template<typename _Index, typename _Value, typename _Vector>
class ParticleTraits<typename GenericSpherical<_Index, _Value, _Vector> > {
public:
    typedef _Index                                             index_type;
    typedef _Value                                             value_type;
    typedef _Vector                                            vector_type;
    typedef typename GenericSpherical<_Index, _Value, _Vector> particle;
    
    ParticleTraits<particle>() 
        : _radius(0), _mass(0), _k_load(0), _k_unload(0), _k0(0),
          _gamma(0), _mu(0) {}
            
    ParticleTraits<particle>(const value_type& radius, 
                             const value_type& mass, 
                             const value_type& k_load, 
                             const value_type& k_unload, 
                             const value_type& k_ratio, 
                             const value_type& gamma, 
                             const value_type& mu)
        : _radius(radius), _mass(mass), _k_load(k_load),
          _k_unload(k_unload), _k_ratio(k_ratio), _gamma(gamma),
          _mu(mu), _inertia(mass*radius*radius), _rcut(2*radius),
          _rcut_square(4*radius*radius) {}
                        
    const value_type& radius(const particle&) const {
        return _radius;
    }
    const value_type& rcut(const particle& p1, const particle& p2) const {
        return _rcut;
    }
    const value_type& rcut_square(const particle& p1, 
                                  const particle& p2) const {
        return _rcut_square;
    }
    const value_type& mass(const particle&) {
        return _mass;
    }
    const value_type& k_load(const particle&) {
        return _inertia;
    }
    const value_type& k_unload(const particle&) {
        return _k_unload;
    }
    const value_type& k_ratio(const particle&) {
        return _k_ratio;
    }
    const value_type& gamma(const particle&) {
        return _gamma;
    }
    const value_type& mu(const particle&) {
        return _mu;
    }
    const value_type& inertia(const particle& p) {
        return _inertia;
    }

private:
    value_type _density;    // "\rho"
    value_type _radius;     // "r"
    value_type _mass;       // "m" (=4/3\pi r^3 \rho)
    value_type _inertia;    // moment of inertia "I_0" (=mr^2)

    // Normal force parameters
    value_type _k_load;     // loading spring constant
    value_type _k_unload;   // unloading spring constant

    // Tangential force (friction) parameters
    value_type _k_ratio;    // ratio between tangential and normal stiffness
    value_type _gamma;      // Mindlin's elastic friction parameter
    value_type _mu;         // coefficient of friction
    
    // Auxiliary values
    value_type _rcut;       // minimum distance for contact
    value_type _rcut_square;// minimum distance for contact squared
};

    
} // dem
} // granular
} // xavier

#endif