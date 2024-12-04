#ifndef __XAVIER_GRANULAR_FORCES_HPP__
#define __XAVIER_GRANULAR_FORCES_HPP__

#include <granular/dem_utils.hpp>

namespace spurt {
namespace granular {
namespace dem {

/*
    Computation at each step
    
    In parallel
    -----------
    1) Loop over all particles i
       1a) identify all the contacts (old and new) between particle i and 
           particle j > i
       1b) compute associated normal and tangential forces and store in
           place for particle i
       1c) history of existing contact is updated in place
       1d) history of new contact (or contact removal) is created and
           stored per-thread
       1e) store force i<->j for particle j per thread
       
    In serial
    ---------
    2) Update global list of contact history by consolidating 
       per-thread history updates
    3) Compute particle force information by consolidating per-thread
       i > j forces
       
    In parallel
    -----------
    4) Update the position of each particle using Newton's second law.
    
    Optional 
    --------
    5) If necessary, recompute the neighbor list of each particle
    6) Update boundary position / geometry
    7) Do a binary dump of current configuration
*/

/** Computation of the forces exerted upon a particle #i by all its neighbors
    in contact and by the boundaries of the domain. Note that to avoid 
    redundant computations in parallel, only forces caused by particles with
    indices higher than i are computed at this stage. Specifically, the forces
    and torques computed here are stored directly in the i-th entry of the 
    supplied _ParticleContainer object while the reciprocal forces are stored
    in what must be a thread-safe data structure.

    Explanations about expected template types:
    - Index is an index type to reference particles and boundaries.
    - _BoundaryContainer is an iterable container of boundaries 
      (e.g., std::list<_Boundary>).
    - _ParticleContainer is a random access iterable container of particles
      (e.g., std::vector<_Particle>).
    - _HistoryContainer is a (_Pair, _History) symbol table where _Pair is a 
      pair of indices corresponding to two objects (particle or boundary) 
      in contact and _History describes the history of their contact force 
      (e.g., std::map<std::pair<_Index, _Index>, _History, std::less<_Pair> >).
    - _ForceModel is a type that supports the computation of normal and 
      tangential force between any two objects in contact (particle or 
      boundary).
    - _ForceContainer is a (_Pair, _Force) symbol table where _Pair is as
      before and _Force is an object that stores both forces and torques 
      between two particles. It is used to store reciprocal forces.
      (e.g., std::map<std::pair<_Index, _Index>, std::pair<_Vector, _Vector> >).
    - _BoundaryTraits supports the determination of the contact between
      particles and boundaries and provides all the information needed to 
      compute the corresponding forces and torques (with supplied history 
      information).
    - _ParticleTraits supports the determination of the contact between
      particles and provides all the information needed to compute the 
      corresponding forces and torques (with supplied history information).
 */

template<typename _BoundaryContainer,
         typename _BoundaryTraits,
         typename _ContactInformation,
         typename _ForceModel,
         typename _ForceContainer,
         typename _HistoryContainer,
         typename _Index,
         typename _NeighborList,
         typename _ParticleContainer,
         typename _ParticleTraits>
void compute_forces(_ForceContainer& forces,
                    const _Index& i,
                    const _BoundaryContainer& boundaries,
                    const _BoundaryTraits& b_traits,
                    const _ParticleContainer& particles,
                    const _ParticleTraits& p_traits,
                    const _ForceModel& model,
                    const _NeighborList& neighbors,
                    _HistoryContainer& old_history,
                    _HistoryContainer& new_history,
                    _HistoryContainer& old_boundary_history,
                    _HistoryContainer& new_boundary_history) 
{
    typedef _BoundaryContainer                          boundary_container;
    typedef _ContactInformation                         contact_information;
    typedef _ForceContainer                             force_container;
    typedef _HistoryContainer                           history_container;
    typedef _Index                                      index;
    typedef _NeighborList                               neighbor_list;
    typedef _ParticleContainer                          particle_container;
    typedef typename _BoundaryContainer::data_type      boundary;
    typedef typename _ForceContainer::data_type         force;
    typedef typename _ForceContainer::key_type          id_pair;
    typedef typename _ForceContainer::value_type        force_info;
    typedef typename _HistoryContainer::data_type       history;
    typedef typename _ParticleContainer::data_type      particle;
    typedef typename particle::vector_type              vector_t;
    typedef typename vector_t::value_type               value_t;
    typedef typename boundary_container::const_iterator const_boundary_iterator;
    typedef typename history_map::const_iterator        const_history_iterator;
    typedef typename history_map::iterator              history_iterator;
    typedef typename neighborlist::const_iterator       const_neighbor_iterator;
    
    particle& pi = particles[i];
    pi.force = 0;
    pi.torque = 0;
    const vector_t& xi = pi.position;
    
    value_t overlap;
    vector_t normal;
    history_iterator hist_it;
        
    // Particle-particle forces
    for (const_neighbor_iterator nit=neighbors.begin() ;
         nit!=neighbors.end() ; ++nit) {
        const index& j = *nit;
        const particle& pj = particles[j];
        contact_information contact_info;
         
        // Quickly rule out irrelevant cases
        // 1. only consider neighbors with higher index
        if (j < i) continue;
        // 2. is there an actual contact?
        if (!model::contact(pi, pj, contact_info, p_traits)) continue;

        id_pair ij(i, j), ji(j, i);
        // is this an existing contact (w/ history)?
        hist_it = old_history.find(ij);
        if (hist_it == old_history.end()) {
            // if not, create a new entry for future reference
            std::pair<history_iterator, bool> r = new_history.insert(ij);
            hist_it = r.first;
        }
        history& hist = *hist_it;
        
        // compute normal force
        value_t fn = model::normal_force(pi, pj, contact_info, p_traits, hist);
        vector_t F = -fn*normal;
        
        // compute tangential force and associated torque
        const particle& pj = particles[j];
        vector_t T = model::tangential_force(pi, pj, contact_info, p_traits, 
                                             hist);
        vector_t _torque = utils::outer_product(normal, T);
        
        // update forces for particle i
        F += T;
        pi.force += F;
        pi.torque += _torque;
        
        // save reciprocal force for particle i acting on particle j:
        // opposite forces, same torque
        forces.insert(force_info(ji, force(-F, _torque)));
    }
    
    // particle-boundary forces
    for (const_boundary_iterator bound_it=boundaries.begin() ; 
         bound_it != boundaries.end() ; ++bound_it) {
        const boundary& bj = *bound_it;
        index j = bj.index();
        contact_information contact_info;
        
        if (!model::contact_boundary(pi, bound, contact_info, p_traits, 
                                     b_traits)) 
            continue;
        
        id_pair ij(i, j), ji(j, i);
        
        hist_it = old_boundary_history.find(ij);
        if (hist_it == old_boundary_history.end()) {
            std::pair<history_iterator, bool> r = 
                new_boundary_history.insert(ij);
            hist_it = r.first;
        }
        history& hist = *hist;
        
        // compute normal force
        value_t fn = model::normal_force_boundary(pi, bj, contact_info,
                                                  p_traits, b_traits, 
                                                  hist);
        vector_t F = -fn*normal
        
        // compute tangential force
        const particle& pj = particles[j];
        vector_t T = model::tangential_force_boundary(pi, bj, contact_info,
                                                      p_traits, b_traits,
                                                      hist);
        F += T;
        pi.force += F;
        pi.torque += utils::outer_product(normal, T);
        
        forces.insert(force_info(ji, force(-F, _torque)));
        
        // NB: no need to save reciprocal force and torque if boundary
        // is assumed unaffected by its interaction with particles and
        // its motion is entirely controled by arbitrary external forces
        // prescribed by the user.
    }
}

template<typename _Index,
         typename _ParticleContainer,
         typename _ParticleTraits,
         typename _Vector>
inline double leapfrog_step(const _Index& i, _ParticleContainer& particles,
                            const _ParticleTraits& particle_traits,
                            const _Vector& gravity, double dt) {
    typedef _Index                                  index;
    typedef _ParticleContainer                      particle_container;
    typedef typename _ParticleContainer::data_type  particle;
    typedef _Vector                                 vector_t;
    typedef typename _ParticleTraits::value_type    value_t;
    
    const value_t& mass    = particle_traits.mass(i);
    const value_t& inertia = particle_traits.inertia(i);
    particle& pi = particles[i];
    pi.velocity += (dt/mass)*pi.forces;
    pi.velocity += dt*gravity;
    vector_t dp = dt*pi.velocity;
    pi.position += dp;
    pi.displacement += dp;
    pi.angular_velocity += (dt/inertia)*pi.torque;
    
    return utils::norm(pi.displacement);
}

    
} // dem
} // granular
} // spurt