#ifndef __ADAPTIVE_MAP_SAMPLING_HPP__
#define __ADAPTIVE_MAP_SAMPLING_HPP__

#include <string>
#include <list>
#include <vector>
#include <math/fixed_vector.hpp>
#include "maps_lib/definitions.hpp"
#include "maps_lib/orbits.hpp"
#include "maps_lib/index.hpp"
#include "poincare/fixpoints.hpp"

namespace map_analysis {
typedef std::vector<nvis::vec2>                         separatrix_type;
typedef std::pair<int, std::vector<separatrix_type> >   separatrix_connection_type;

struct adaptive_map_sampling_params {

    adaptive_map_sampling_params()
        : hdf_name("/Users/xmt/data/NIMROD/CDXUB7/h5/cdxub7n.h5"), time_step("1000"),
          out_name("/Users/xmt/tmp/cdxub7"), eps(1.0e-6), n(5000), m(50), dq(0.05),
          mr(2.), mt(150000), ma(0.01), Ma(1.), mp(20), pm(1), err(0.1), dx(1.0e-4),
          dtheta(0.75*M_PI), close_d(1.0) {}
          
    // hdf_name:    input HF5 file name
    // time_step:   corresponding time step
    // out_name:    output file name
    // eps:         integration precision
    // dq:          period assessment precision (will be enforced)
    // mr:          max triangle aspect ratio
    // ma/Ma        min/max triangle area
    // err:         max relative approximation error
    // n:           number of sampled positions in phase portrait (!= number of sampled orbits)
    // m:           number of iterations of the map
    // mt:          max number of triangles per period
    // mp:          max considered period in analysis
    // pm:          min considered period in analysis
    // dx:          smallest allowable sampling distance in index computation
    // dtheta:      maximum angular distance allowed in tracking
    // close_d:     close distance between chains of same type
    
    std::string hdf_name, time_step, out_name;
    double eps, dq, mr, ma, Ma, err, dx, dtheta, close_d;
    unsigned int n, m, mt, mp, pm;
};

struct adaptive_map_sampling_output {

    // base_mesh:           initial triangulation from which all others are derived
    // p_orbits:            index of p-relevant orbits within central_map_orbits
    // p_meshes:            triangulations resulting from the adaptive sampling of the p-map
    // p_cand_tris:         index of triangles that have a non-zero naive (linear) poincare index for period p
    // p_sing_tris:         index of triangles that actually contain a p-fixed point (currently subset of p_cand_tris)
    // p_prob_tris:         index of triangles that posed problems
    // p_index_vectors:     vectors computed in the course of determining the Poincare index of a triangle
    // p_chains:            chains of fixed points of period p
    // p_prob_edges:        index of edges whose angle coordinates in p-map could not be tracked
    // p_degenerate_points: points where index computation failed
    // p_separatrices:      separatrices computed along saddle chains
    // p_rejected_separatrices: separatrices that were rejected during the construction
    
    mesh_type                                               base_mesh;
    std::vector<mesh_type>                                  p_meshes;
    std::vector<std::list<unsigned int> >                   p_orbits;
    std::vector<std::list<unsigned int> >                   p_cand_tris;
    std::vector<std::list<unsigned int> >                   p_sing_tris;
    std::vector<std::list<unsigned int> >                   p_prob_tris;
    std::vector<std::list<step_type> >                      p_index_vectors;
    std::vector<std::list<std::list<xavier::fixpoint> > >   p_chains;
    std::vector<std::list<xavier::Edge> >                   p_prob_edges;
    std::vector<std::list<nvis::vec2> >                     p_degenerate_points;
    std::vector<std::list<separatrix_connection_type> >     p_separatrices;
    std::vector<std::list<separatrix_type> >                p_rejected_separatrices;
};

void adaptive_map_sampling(adaptive_map_sampling_output& output,
                           const adaptive_map_sampling_params& params);
                           
} // map_analysis




#endif






