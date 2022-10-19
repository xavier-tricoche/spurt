#ifndef __XAVIER_GRANULAR_NEIGHBORS_HPP__
#define __XAVIER_GRANULAR_NEIGHBORS_HPP__

#include <set>
#include <vector>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace xavier {
namespace granular {
namespace dem {

template<typename _Index, typename _Particle>
struct CellGrid {
    typedef _Index            index;
    typedef std::set<_Index>  bucket;
    typedef _Particle         particle;
    
    size_t coord_to_index(const nvis::uvec3& c) const {
        return c[0] + resolution[0]*(c[1] + resolution[2]*c[2]);
    }
    
    nvis::uvec3 pos_to_coord(const nvis::vec3& x) const {
        const nvis::vec3& min = bounds.min();
        nvis::vec3 p = x - min;
        p /= step;
        return nvis::uvec3((size_t)floor(p[0]);
                           (size_t)floor(p[1]);
                           (size_t)floor(p[2]));
    }
    
    bool check_coord(int& i, int dim) const {
        const bool& per = periodic[dim];
        if (i>=0 && i<resolution[dim]) return true;
        else if (!periodic[dim]) return false;
        else if (i<0) i = resolution[dim]-1;
        else i = 0;
        return true;
    }
    
    CellGrid(const nvis::bbox3& _bounds, const nvis::bvec3& _periodic,
             const nvis::uvec3& _resolution, 
             const std::vector<particle>& particles) 
        : bounds(_bounds), periodic(_periodic), resolution(_resolution) {
        
        const nvis::vec3& min = bounds.min();
        buckets.resize(resolution[0]*resolution[1]*resolution[2]);
        
        step = bounds.size() / nvis::vec3(resolution - nvis::uvec3(1));
        
        std::vector<particle>::const_iterator it;
        for (it=particles.begin() ; it!=particles.end() ; ++it) {
            size_t n = coord_to_index(pos_to_coord(it->position));
            buckets[n].insert(it->index);
        }
    }
    
    void neighbor_candidates(std::list<index>& candidates,
                             const nvis::vec3& p) const {
        candidates.clear();
        
        nvis::uvec3 c = pos_to_coord(p);
        for (int dz=-1; dz<2 ; ++dz) {
            int z = (int)c[2] + dz;
            if (!check_coord(z, 2)) continue;
            for (int dy=-1 ; dy<2 ; ++dy) {
                int y = (int)c[1] + dy;
                if (!check_coord(y, 1)) continue;
                for (int dx=-1 ; dx<2 ; ++dx) {
                    int x = (int)c[0] + dx;
                    if (!check_coord(x, 0)) continue;
                    nvis::uvec3 d(x, y, z);
                    n = coord_to_index(d);
                    std::copy(buckets[n].begin(), buckets[n].back(), 
                              std::back_inserter(candidates)); 
                }
            }
        }
    }
    
    template<typename Container>
    void remove_particles(const Container& removed, 
                          const nvis::uvec3& bucket_id) {
         typedef typename Container::const_iterator const_iterator;
         for (const_iterator it=removed.begin() ; it!=removed.end() ; ++it) {
             buckets[bucket_id].erase(*it);
         }
    }
    
    template<typename Container>
    void add_particles(const Container& added, 
                       const nvis::uvec3& bucket_id) {
         typedef typename Container::const_iterator const_iterator;
         for (const_iterator it=added.begin() ; it!=added.end() ; ++it) {
             buckets[bucket_id].insert(*it);
         }
    }
    
    nvis::uvec3 resolution;
    nvis::vec3  step;
    nvis::bbox3 bounds;
    nvis::bvec3 periodic;
    std::vector<bucket> buckets;
};

template<typename _Index, typename _Particle>
struct NeighborList : public std::list<_Index> {
    typedef _Index                      index;
    typedef std::list<_Index>           base_type;
    typedef base_type::const_iterator   const_iterator;
    typedef base_type::iterator         iterator;
    typedef _Particle                   particle;
    
    NeighborList(const index& i, 
                 const typename CellGrid<index, particle>& cell_grid,
                 const std::vector<particle>& particles, 
                 double rm_sq) {
        const nvis::vec3& xi = particles[i].position;
        cell_grid.neighbors_candidates(*this, xi);
        for (iterator it=this->begin() ; it!=this->end() ; ) {
            nvis::vec3 rij = particles[*it].position;
            rij -= xi;
            if (nvis::inner(rij, rij) > rm_sq) {
                it = this->erase(it);
            }
            else ++it;
        }
    }
};

} // dem
} // granular
} // xavier

#endif 