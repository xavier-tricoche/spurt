#ifndef __TOKATOPO_HPP__
#define __TOKATOPO_HPP__

#include "maps/triangulation.hpp"
#include "maps/basic_definitions.hpp"

#include <poincare/fixpoints.hpp>
#include <poincare/metric.hpp>

namespace tokatopo {
typedef std::pair<double, double>                                       pair_type;
typedef std::vector<nvis::vec2>::const_iterator                         const_iter_type;
typedef spurt::point_data                                              point_data;
typedef spurt::triangulation<point_data, spurt::point_locator>        mesh_type;
typedef mesh_type::index_type                                           index_type;
typedef mesh_type::triangle_type                                        triangle_type;
typedef mesh_type::point_type                                           point_type;
typedef boost::rational<int>                                            rational_type;
typedef spurt::interval<double>                                        interval_type;

extern int res, maxp, maxit, n1, n2, it, m, maxround, max_nb_triangles, max_period, period;
extern char* outs, *file, *ts;
extern double h, dq, min_area, max_area, max_ratio, close_enough;

extern mesh_type* sampling_mesh;
extern std::vector< std::vector< spurt::fixpoint > > all_chains;
extern spurt::map_metric metric;

extern std::vector<std::pair<unsigned int, nvis::vec2> > probed_seeds;

extern unsigned int verbose_level;

double min_norm(unsigned int& pid, unsigned int oid, const spurt::map_metric& metric, unsigned int p);

void initialize(int argc, char* argv[]);

int compute_fixed_points();
}

inline double tokatopo::
min_norm(unsigned int& pid, unsigned int oid, const spurt::map_metric& metric, unsigned int p)
{
    if (oid>=spurt::__map_orbits.size()) {
        throw std::runtime_error("min_norm: requested min norm on invalid orbit");
    }
    const spurt::orbit& orbit = spurt::__map_orbits[oid];
    if (p >= orbit.size()) {
        return std::numeric_limits<double>::max();
    }
    
    double min = std::numeric_limits<double>::max();
    for (int id = 0 ; id < orbit.size() ; ++id) {
        double norm = nvis::norm(spurt::vector_value(oid, id, p, metric));
        if (norm<min) {
            min=norm;
            pid = id;
        }
    }
    
    return min;
}


#endif





