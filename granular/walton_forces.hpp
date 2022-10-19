#ifndef __XAVIER_GRANULAR_WALTON_FORCES_HPP__
#define __XAVIER_GRANULAR_WALTON_FORCES_HPP__

#include <granular/dem_utils.hpp>

namespace xavier {
namespace granular {
namespace walton {

template<typename _Value, typename _Vector>
struct ContactHistory {
    typedef _Vector  vector_type;
    typedef _Value   value_type;
    
    ContactHistory() 
        : T_old(0), alpha0(0), fn_old(0), t_star(0), k0(0) {}
    
    vector_type   T_old;  // last friction vector
    value_type    alpha0; // zero unloading force displacement
    value_type    fn_old; // last normal force magnitude
    value_type    t_star; // reference friction value
    value_type    k0;     // initial tangential stiffness
};

template<typename _Index, typename _Value, typename _Vector>
struct ContactInformation {
    typedef _Index  index_type;
    typedef _Vector vector_type;
    typedef _Value  value_type;
    
    ContactInformation() 
        : i(0), j(0), position(0), normal(0), overlap(0), normal_force(0) {}
    
    index_type   i, j;
    vector_type  position;
    vector_type  normal;
    value_type   overlap;
    value_type   normal_force;
};

// default wrapper for DEM interaction with Walton's force model
template<typename _Index, typename _Value, typename _Vector,
         typename _Particle, typename _ParticleTraits,
         typename _Boundary, typename _BoundaryTraits,
         typename _History>
struct ForceModel {
    typedef _Index                                       index_type;
    typedef _Value                                       value_type;
    typedef _Vector                                      vector_type;
    typedef _Particle                                    particle;
    typedef _ParticleTraits                              particle_traits;
    typedef _Boundary                                    boundary;
    typedef _BoundaryTraits                              boundary_traits;
    typedef ContactHistory<_Value, _Vector>              history;
    typedef ContactInformation<_Index, _Value, _Vector>  contact_information;
    
    static inline
    bool contact(const particle& pi, const particle& pj,
                 contact_information& info, const particle_traits& traits) {
        const value_type& ri = p_traits.radius(pi);
        const value_type& rj = p_traits.radius(pj);
        value_type rcut   = p_traits.rcut(pi, pj);
        value_type rcutsq = p_traits.rcut_square(pi, pj);
        
        for (int n=0 ; n<3 ; ++n) {
            info.normal[n] = pj[n] - pi[n];
            if (fabs(info.normal[n]) > rcut) return false;
        }
        value_type dist_sq = 
            utils::norm_sq<vector_type, value_type>(info.normal);
        if (dist_sq > rcutsq) return false;
        
        value_type dist = sqrt(dist_sq);
        info.normal *= 1/dist;
        info.overlap = rcut - dist;
        info.position = 0.5*(pi.position + pj.position); // only true if monodisperse
        info.i = pi.index;
        info.j = pj.index;
        info.normal_force = 0;
        info.k0 = p_traits.k_ratio(pi) * p_traits.k_load(pi);
        return true;
    }
    
    static inline
    bool contact_boundary(const particle& p, const boundary& b,
                          contact_information& info, 
                          const particle_traits& p_traits,
                          const boundary_traits& b_traits) {
        const value_type& r = p_traits.radius(p);
        const value_type rcut = r;
        
        vector_type normal;
        value_type dist = b.distance(p.position, normal);
        if (dist > rcut) return false;
        
        // bp is tangential contact point between p and b
        // (i.e., normal projection of p onto surface)
        vector_type bp = p - normal;
        
        info.normal *= 1/dist;
        info.overlap = rcut - dist;
        info.position = bp;
        info.i = pi.index;
        info.j = b.index;
        info.normal_force = 0;
        info.k0 = p_traits.k_ratio(p)*p_traits.k_load(p);
        return true;
    }
    
    static inline value_type
    normal_force(const particle& pi, const particle& pj, 
                 contact_information& info, const particle_traits& p_traits,
                 history& h) {
        info.fn = walton::normal_force<vector_t, value_t>
                            (info.overlap, h.alpha0, info.k0,
                             p_traits.k_load(pi), p_traits.k_unload(pi),
                             p_traits.k_ratio(pi));
        return info.fn;
    }

    static inline value_type
    normal_force_boundary(const particle& p, const boundary& bj,
                          contact_information& info,
                          const particle_traits& p_traits,
                          const boundary_traits& b_traits, 
                          history& h) {
        info.fn = walton::normal_force<vector_t, value_t>
                            (info.overlap, h.alpha0, info.k0,
                             p_traits.k_load(p), p_traits.k_unload(p),
                             p_traits.k_ratio(p));
        return info.fn;
    }
    
    static inline vector_type
    tangential_force(const particle& pi, const particle& pj, 
                     contact_information& info,
                     const particle_traits& p_traits, 
                     history& h) {
        return walton::tangential_force<vector_t, value_t>
                (h.T_old, info.normal, (pj.velocity - pi.velocity), 
                 0.5*(p_traits.radius(pi)*pi.angular_velocity + 
                      p_traits.radius(pj)*pj.angular_velocity),
                 info.normal_force, h.fn_old, info.k0, 
                 p_traits.mu(pi), p_traits.gamma(pi), h.tstar);
                                
    }
    
    static inline vector_type 
    tangential_force_boundary(const particle& p, const boundary& bj,
                              contact_information& info,
                              const particle_traits& p_traits,
                              const boundary_traits& b_traits, 
                              history& h) {
        return walton::tangential_force<vector_t, value_t>
                (h.T_old, info.normal, 
                 (b_traits.velocity(bj, contact_info) - pi.velocity), 
                 p_traits.radius(p)*p.angular_velocity,
                 info.normal_force, h.fn_old, info.k0, p_traits.mu(p),
                 p_traits.gamma(p), h.tstar);
    }
};

// Walton & Braun (1986) normal force model
// F_N = K_load * overlap               % loading
// F_N = K_unload * (overlap - alpha0)  % unloading
// with alpha0 = max_overlap*(1-e^2) chosen to ensure continuity 
// at turning point
template<typename _Vector, typename _Value = typename _Vector::value_type>
inline _Value normal_force(_Value overlap, _Value& alpha0, _Value& k0,
                           const _Value& kn_load, const _Value& kn_unload,
                           const _Value& kt_over_kn) {
    if (overlap > alpha0) { // loading phase
        alpha0 = overlap; // initially alpha0 stores overlap history
        return kn_load*overlap;
    }
    else { // unloading phase
        alpha0 *= (1-kn_load/kn_unload);
        return kn_unload*(overlap-alpha0);
    }
}

// Walton (1992) tangential force model
// Step 1: Initialize friction as norm-preserving projection of last 
//         tangential force vector onto contact tangent plane
template<typename _Vector, typename _Value = typename _Vector::value_type>
inline _Vector T_init(const _Vector& T_old, 
                      const _Vector& Normal) {
    _Vector T_new(project_on_plane(T_old, Normal));
    _Value t_new_normsq = utils::norm_sq(T_new);
    _Value t_old_normsq = utils::norm_sq(T_old);
    T_new *= sqrt(t_old_normsq/t_new_normsq);
    return T_new;
}

// Step 2: Project relative surface displacement onto contact tangent plane
template<typename _Vector, typename _Value = typename _Vector::value_type>
inline _Vector Delta_s(const _Vector& Normal, 
                       const _Vector& Vij, 
                       const _Vector& ROmega_avg, 
                       _Value dt) {
    _Vector Ds(project_on_plane(Vij, Normal));
    Ds += utils::outer_product(ROmega_avg, Normal);
    Ds *= dt;
    return Ds;                                          
}

// Step 3: update T* in expression of tangential stiffness
template<typename _Value>
inline _Value scaled_T_star(_Value t_star, _Value fn_old, _Value fn) {
    return t_star*fabs(fn/fn_old);
}

// Step 4: Compute tangential stiffness Kt
template<typename _Value>
inline _Value k_tang(_Value t, _Value t_star, 
                     _Value k0, _Value mu, 
                     _Value fn, 
                     _Value gamma,
                     bool increasing) {
    if (increasing) {
        return k0*pow(1. - (t - t_star)/(mu*fn - t_star), gamma);
    }
    else {
        return k0*pow(1. - (t_star - t)/(mu*fn + t_star), gamma);
    }
}

// Step 5: Compute component of friction force parallel to old friction force
template<typename _Vector, typename _Value = typename _Vector::value_type>
inline _Vector compute_T_parallel(const _Vector& T, const _Value& kt, 
                                  const _Vector& Delta_s,
                                  _Vector& Delta_s_orthogonal,
                                  bool& sign_changed) {
    _Value t_norm = utils::norm(T);
    _Value delta_s_parallel = utils::inner_product(Delta_s, T)/t_norm;
    Delta_s_orthogonal = Delta_s - (delta_s_parallel/t_norm)*T;
    value_type t_parallel = t_norm + kt*delta_s_parallel;
    if (delta_s_parallel < 0 && t_parallel < 0) sign_changed = true;
    return (t_parallel/t_norm)*T;
}

// Complete tangential force model
template<typename _Vector, typename _Value = typename _Vector::value_type>
_Vector tangential_force(const _Vector& T_old, const _Vector& Normal,
                         const _Vector& Vij, const _Vector& ROmega_avg,
                         _Value fn, _Value fn_old, _Value k0, _Value mu, 
                         _Value gamma, _Value dt, _Value& t_star) {
    
    // compute initial tangential friction (T)
    _Vector T(compute_T_init(T_old, Normal));
    
    // compute relative surface displacement (delta_s)
    _Vector Delta_s(compute_delta_s(Normal, Vij, ROmega_avg, dt));
    
    // compute effective tangential stiffness in direction parallel to
    // existing friction force (Kt)
    t_star = update_T_star(t_star, fn_old, fn);
    _Value t = utils::norm(T);
    _Value t_old = utils::norm(T_old);
    _Value Kt = compute_Kt(T, T_star, K0, mu, Fn, gamma, (t >= t_old));
    
    // compute friction force component parallel to existing friction force
    bool signed_changed = false;
    _Vector Delta_s_orthogonal;
    _Vector T_new(compute_T_parallel(T, kt, Delta_s, Delta_s_orthogonal, 
                                    sign_changed));
    if (sign_changed) t_star *= -1;
    
    // add friction force component orthogonal to existing friction force
    T_new += k0*Delta_s_orthogonal;
    _Value t_new = utils::norm(T_new);
    _Value ratio = mu*fn/t_new;
    
    // clamp magnitude at friction limit
    if (ratio < 1) return T_new *= ratio;
    return T_new;
}

} // walton
} // granular
} // xavier

#endif