#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <math.h>
#include <sstream>
#include <fstream>

// boost
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>


// #include <dfftw.h>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>
#include <math/bounding_box.hpp>

// spurt
#include <math/math.hpp>
#include <poincare/metric.hpp>

// christoph
#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <tokamak/poincare_map.hpp>

// teem
#include <teem/nrrd.h>

// OpenMP
#if _OPENMP
#include <omp.h>
#endif

// meshing
#include <maps/basic_definitions.hpp>
#include <maps/triangulation.hpp>
#include <maps/adaptive_triangulation.hpp>
#include <maps/quality_control.hpp>
#include <maps/IO.hpp>

#include <maps/period.hpp>
#include <data/quadtree.hpp>

#define __REFINE_ZEROS__

int res, maxp, maxit, n1, n2, it, m, maxround, max_nb_triangles, max_period, period;
char* outs, *file, *ts;
double h, max_diff, min_area, max_area, max_ratio;

struct point_counter {
    typedef spurt::quadtree<int, nvis::vec2>   quadtree_type;
    
    point_counter() : n(0) {}
    
    void operator()(const quadtree_type& q) {
        n += q.data().size();
    }
    
    int n;
};

struct leaf_display {
    typedef spurt::quadtree<int, nvis::vec2>   quadtree_type;
    
    leaf_display() : os() {}
    void operator()(const quadtree_type& q) {
        os << q << '\n';
    }
    
    void reset() {
        os.str("");
        os.clear();
    }
    
    std::ostringstream os;
};

struct point_locator {
    typedef spurt::quadtree<int, nvis::vec2>   quadtree_type;
    typedef quadtree_type::data_type            data_type;
    
    point_locator(const nvis::bbox2& bounds, size_t max_depth, size_t max_nb_pts)
        : qt(bounds, max_depth, max_nb_pts) {}
        
    point_locator(const point_locator& pl) : qt(pl.qt) {}
    
    int find_close_point(const nvis::vec2& x, int id) {
        quadtree_type& leaf = qt.find(x);
        int res = (leaf.data().size() ? leaf.data()[0].second : 0);
        leaf.insert(x, id);
        return res;
    }
    
    void display() const {
        std::string str("");
        display(str);
    }
    
    void display(const std::string& preamble) const {
        leaf_display display;
        qt.for_each_leaf(display);
        std::cerr << preamble << "quadtree contains: " << display.os.str() << std::endl;
    }
    
    quadtree_type qt;
};

template<typename I>
struct orbit_integrator {
    typedef I   integrator_type;
    
    orbit_integrator(const integrator_type& integrator, size_t nsteps,
                     const spurt::map_metric& metric)
        : __nsteps(nsteps), __integ(integrator), __metric(metric) {}
        
    void operator()(const nvis::vec2& x0,
                    std::vector<nvis::vec2>& points,
                    std::vector<spurt::point_data>& data) const {
        const integrator_type* pmap = __integ.clone();
        points.clear();
        data.clear();
        std::vector<nvis::vec2> __steps;
        try {
            pmap->map(x0, __steps, __nsteps);
        } catch (...) {
            std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
            return;
        }
        if (__steps.empty()) {
            std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
            return;
        }
        
        std::pair<double, double> sf = spurt::period_x_periodic(__steps, x0, __metric);
        
        points.resize(__steps.size() + 1);
        data.resize(__steps.size() + 1);
        size_t orbit_id = spurt::__map_orbits.size();
        points[0] = __metric.modulo(x0);
        data[0] = spurt::point_data(orbit_id, 0);
        for (int i = 0 ; i < __steps.size() ; ++i) {
            points[i+1] = __metric.modulo(__steps[i]);
            data[i+1] = spurt::point_data(orbit_id, i + 1);
        }
        spurt::__map_orbits.push_back(spurt::orbit(points, sf.first));
    }
    
    size_t                  __nsteps;
    const integrator_type&  __integ;
    spurt::map_metric      __metric;
};


struct vector_value_function {
    vector_value_function(int period, double max_q_dist, const spurt::map_metric& metric)
        : __period(period), __max_q_dist(max_q_dist), __metric(metric) {}
        
    int order() const {
        return 1;
    }
    
    std::string name() const {
        std::ostringstream os;
        os << __period << "-map";
        return os.str();
    }
    
    bool is_valid(const spurt::point_data& d) const {
        return true;
    }
    
    std::string value_string(const spurt::point_data& p) const {
        double q = p.period();
        double mind = std::numeric_limits<double>::max();
        for (int i = 1 ; i < __period ; ++i) {
            mind = std::min(fabs((double)__period / (double)i - q), mind);
        }
        nvis::vec2 v;
        if (mind > __max_q_dist)
            // out of range. set 0 (invalid) value
        {
            v = nvis::vec2(0, 0);
        } else {
            v = vector_value(p, __period, __metric);
        }
        std::ostringstream os;
        os << v[0] << " " << v[1] << " 0.";
        return os.str();
    }
    
    int                 __period;
    double              __max_q_dist;
    spurt::map_metric  __metric;
};


void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "f",      "file",                 airTypeString,  1, 1, &file,                NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",                 airTypeString,  1, 1, &ts,                  NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &outs,                NULL,       "output nrrd file name");
    hestOptAdd(&hopt, "eps",    "eps",                  airTypeDouble,  0, 1, &h,                   "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "n1",     "# 1D samples",         airTypeInt,     0, 1, &n1,                  "100",      "number of 1D samples");
    // hestOptAdd(&hopt, "n2",  "# 2D samples",         airTypeInt,     0, 1, &n2,                  "2000",     "number of 2D samples");
    hestOptAdd(&hopt, "m",      "# rotations",          airTypeInt,     0, 1, &m,                   "50",       "number of rotations");
    hestOptAdd(&hopt, "d",      "safety factor step",   airTypeDouble,  0, 1, &max_diff,            "0.01",     "maximum discrepancy between consecutive safety factors");
    // hestOptAdd(&hopt, "mr",  "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,           "3.",       "maximum triangle aspect ratio");
    // hestOptAdd(&hopt, "nr",  "max # rounds",         airTypeInt,     0, 1, &maxround,            "10",       "max number of 1D approximation iterations");
    hestOptAdd(&hopt, "mt",     "max # triangles",      airTypeInt,     0, 1, &max_nb_triangles,    "500000",   "max number of triangles in adaptive sampling");
    hestOptAdd(&hopt, "ma",     "min triangle area",    airTypeDouble,  1, 1, &min_area,            NULL,       "min triangle area in adaptive sampling");
    hestOptAdd(&hopt, "Ma",     "max triangle area",    airTypeDouble,  1, 1, &max_area,            NULL,       "max triangle area in adaptive sampling");
    // hestOptAdd(&hopt, "mp",      "max period",           airTypeInt,     0, 1, &max_period,          "25",       "max considered period in adaptive refinement");
    hestOptAdd(&hopt, "p",      "period",               airTypeInt,     1, 1, &period,              NULL,       "targeted period");
    
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Compute approximative q-profile using magnetic field integration",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef std::pair<double, double>                                       pair_type;
typedef std::vector<nvis::vec2>::const_iterator                         const_iter_type;
typedef spurt::point_data                                              point_data;
typedef spurt::triangulation<point_data, point_locator>                mesh_type;
typedef mesh_type::index_type                                           index_type;
typedef mesh_type::triangle_type                                        triangle_type;
typedef mesh_type::point_type                                           point_type;

std::vector<boost::rational<int> > ref_values;

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(file), std::string(ts));
    field->periodic_coordinates(false);
    
    bool per[2] = {true, false};
    spurt::map_metric metric(field->bounds(), per);
    
    poincare_map pmap(field);
    pmap.precision(h);
    orbit_integrator<poincare_map> intg(pmap, m, metric);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double width = diagonal[0];
    double height = diagonal[1];
    double dh = height / (double)(n1 - 1);
    double x_center = bounds.center()[0];
    
    typedef boost::rational<int>    rational_type;
    std::set<rational_type> rationals;
    for (int den = 1 ; den <= period ; ++den) {
        rationals.insert(rational_type(period, den));
    }
    
    ref_values.clear();
    std::copy(rationals.begin(), rationals.end(),
              std::back_inserter(ref_values));
    std::cerr << "there are " << ref_values.size() << " interesting values\n";
    for (int i = 0 ; i < ref_values.size() ; ++i) {
        std::cerr << ref_values[i] << " = " << spurt::value(ref_values[i]) << '\n';
    }
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    srand48(time(0));
    std::vector< std::pair<point_type, spurt::point_data> > all_points;
    
    // reset central orbit repository
    spurt::__map_orbits.clear();
    
    // 1D sampling along symmetry axis
#pragma openmp parallel
    {
        std::vector<nvis::vec2> steps;
        std::vector<spurt::point_data> data;
        #pragma omp for schedule(dynamic,1)
        for (int i = 0 ; i < n1 ; ++i) {
#if _OPENMP
            const int thread = omp_get_thread_num();
            if (!thread) {
                std::cout << i << " from " << n1 << " values computed\n";
            }
#endif
            nvis::vec2 x0 = bounds.min() + nvis::vec2(0.5 * diagonal[0], i * dh);
            x0[0] += 0.05 * dh * (-1. + 2.*drand48());
            intg(x0, steps, data);
            if (steps.empty()) {
                continue;
            }
            
            for (int i = 0 ; i < steps.size() ; ++i) {
                all_points.push_back(std::make_pair(metric.modulo(steps[i]),
                                                    data[i]));
            }
        }
    }
    std::cout << all_points.size() << " points sampled after 1D seeding\n";
    
    std::cerr << "building initial triangulation after 1D sampling\n";
    // form the boundary
    std::pair<nvis::vec2, spurt::point_data> boundary[4];
    std::vector<nvis::vec2> corner(4);
    corner[0] = bounds.min();
    corner[1] = bounds.min() + nvis::vec2(bounds.size()[0], 0);
    corner[2] = bounds.max();
    corner[3] = bounds.min() + nvis::vec2(0, bounds.size()[1]);
    size_t orbit_id = spurt::__map_orbits.size();
    spurt::__map_orbits.push_back(spurt::orbit(corner, 0));
    for (int i = 0 ; i < 4 ; ++i) {
        boundary[i].first = corner[i];
        boundary[i].second = spurt::point_data(orbit_id, i);
    }
    
    mesh_type mesh(boundary, point_locator(bounds, 10, 10));
    
    // std::cerr << "before triangulating: " << std::endl;
    // const std::vector<spurt::orbit>& all_orbits = spurt::__map_orbits;
    // for (int i = 0 ; i < all_orbits.size() ; ++i) {
    //  std::cerr << "orbit #" << i
    //  << " contains " << all_orbits[i].size() << " points starting at " << all_orbits[i][0]
    //  << " and its period is " << all_orbits[i].period() << std::endl;
    // }
    // for (int i = 0 ; i < all_points.size() ; ++i) {
    //  std::cerr << "point at " << all_points[i].first << " has orbit #"
    //  << all_points[i].second.orbit_id() << std::endl;
    // }
    
    mesh.set_export_basename("1D_sampling_triangulation");
    try {
        mesh.insert_points(all_points);
    } catch (...) {
        export_VTK(mesh, "bug.vtk", "this is faulty triangulation", true, false);
        return -1;
    }
    
    // std::cerr << "after triangulating: " << std::endl;
    // for (int i = 0 ; i < all_orbits.size() ; ++i) {
    //  std::cerr << "orbit #" << i
    //  << " contains " << all_orbits[i].size() << " points starting at " << all_orbits[i][0]
    //  << " and its period is " << all_orbits[i].period() << std::endl;
    // }
    // for (int i = 0 ; i < mesh.get_nb_vertices() ; ++i) {
    //  std::cerr << "point at " << mesh.get_vertex(i) << " has orbit #"
    //  << mesh.get_data(i).orbit_id() << std::endl;
    // }
    
    
    std::cerr << "triangulation initially contains " << mesh.get_nb_triangles()
              << " triangles\n";
              
    // spurt::export_VTK(mesh, "initial_triangulation.vtk", "none", false);
    
    std::cerr << "refining triangulation to achieve prescribed max area\n";
    mesh.set_export_basename("area_triangulation");
    bool valid_mesh = spurt::refine(mesh, intg,
                                     spurt::triangle_area_controller(max_area),
                                     max_nb_triangles);
    std::cerr << "refinement was " << (valid_mesh ? "successful" : "not successful") << '\n';
    std::cout << mesh.get_nb_vertices() << " points sampled after area-driven 2D refinement\n";
    
    for (std::set<rational_type>::const_iterator iter = rationals.begin() ; iter != rationals.end() ; ++iter) {
        spurt::interval<double> valid(spurt::value(*iter) - max_diff,
                                       spurt::value(*iter) + max_diff);
        std::ostringstream os;
        os << "q=" << iter->numerator() << "/" << iter->denominator()
           << "_triangulation";
        mesh.set_export_basename(os.str());
        std::cerr << "refining triangulation around q=" << spurt::value(*iter) << '\n';
        valid_mesh = spurt::refine(mesh, intg,
                                    spurt::triangle_interval_inclusion_controller(valid, min_area),
                                    max_nb_triangles);
        std::cerr << "refinement was " << (valid_mesh ? "successful" : "not successful") << '\n';
        std::cout << mesh.get_nb_vertices() << " points sampled after value-driven 2D refinement\n";
    }
    
#ifdef __REFINE_ZEROS__
    {
        std::ostringstream os;
        os << "zero_period=" << period << "_triangulation";
        mesh.set_export_basename(os.str());
        std::cerr << "refining triangulation for period " << period << '\n';
        valid_mesh = spurt::refine(mesh, intg,
                                    spurt::triangle_contains_fixed_point_controller(min_area, period, metric),
                                    max_nb_triangles);
        std::cerr << "refinement was " << (valid_mesh ? "successful" : "not successful") << '\n';
        std::cout << mesh.get_nb_vertices() << " points sampled after zero-driven 2D refinement\n";
    }
#endif
    
    std::ostringstream os_name, os_comment;
    os_name << outs << "-vectors-n1=" << n1 << "-n2=" << n2
            << "-m=" << m << "-p=" << period << ".vtk";
    os_comment << "Sample points from " << file << " at t=" << ts << " with estimated safety factor";
    spurt::export_VTK<mesh_type, vector_value_function>(mesh, os_name.str(), os_comment.str(),
            vector_value_function(period, max_diff, metric));
}




















































































































































































































