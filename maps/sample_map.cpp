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

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>
#include <math/bounding_box.hpp>

// map analysis
#include <math/math.hpp>
#include <poincare/metric.hpp>
#include <poincare/newton.hpp>
#include <maps/period.hpp>
#include <maps/regular_mesh.hpp>
#include <poincare/topology.hpp>
#include <poincare/fixpoints.hpp>
#include <poincare/macros.hpp>
#include <maps/sampling.hpp>
#include <maps/temporary.hpp>

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

#include <maps/tokatopo.hpp>

using namespace tokatopo;

namespace {
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
            // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
            return;
        }
        if (__steps.empty()) {
            // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
            return;
        }
        
        std::pair<double, double> sf = spurt::period_x_periodic(__steps, x0, __metric);
        
        // std::cerr << "found period = " << sf.first << std::endl;
        
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
    
    void precision(double h) {
        __integ.precision(h);
    }
    
    size_t                  __nsteps;
    const integrator_type&  __integ;
    spurt::map_metric      __metric;
};
}

void init(int argc, char* argv[])
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
    hestOptAdd(&hopt, "m",      "# rotations",          airTypeInt,     0, 1, &m,                   "50",       "number of rotations");
    hestOptAdd(&hopt, "dq",     "safety factor step",   airTypeDouble,  0, 1, &dq,                  "0.01",     "maximum discrepancy between consecutive safety factors");
    hestOptAdd(&hopt, "mr",     "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,           "2.",       "maximum triangle aspect ratio");
    hestOptAdd(&hopt, "mt",     "max # triangles",      airTypeInt,     0, 1, &max_nb_triangles,    "500000",   "max number of triangles in adaptive sampling");
    hestOptAdd(&hopt, "ma",     "min triangle area",    airTypeDouble,  1, 1, &min_area,            NULL,       "min triangle area in adaptive sampling");
    hestOptAdd(&hopt, "Ma",     "max triangle area",    airTypeDouble,  1, 1, &max_area,            NULL,       "max triangle area in adaptive sampling");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0, 1, &max_period,          "25",       "max considered period in fixed point search");
    // hestOptAdd(&hopt, "p",       "period",               airTypeInt,     1, 1, &period,              NULL,       "targeted period");
    
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Compute approximative q-profile using magnetic field integration",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

std::vector<rational_type> ref_values;

template<typename T>
inline bool has_zero(int period, const T& mesh, int cellid)
{
    typedef T                                               triangulation_type;
    typedef typename triangulation_type::triangle_type      triangle_type;
    typedef typename triangulation_type::index_type         index_type;
    typedef typename triangulation_type::point_type         point_type;
    typedef typename triangulation_type::data_type          data_type;
    
    const triangle_type& tri = mesh.get_triangle_vertices(cellid);
    nvis::vec2 v[3];
    for (int i = 0 ; i < 3 ; ++i) {
        const data_type& d = mesh.get_data(tri[i]);
        v[i] = spurt::vector_value(d, period, spurt::__default_metric);
    }
    
    nvis::vec2 col[2] = { v[0] - v[2], v[1] - v[2] };
    double denom, beta[3];
    denom = col[0][0] * col[1][1] - col[1][0] * col[0][1];
    if (denom == 0) {
        return false;
    }
    beta[0] = (-v[2][0] * col[1][1] + v[2][1] * col[1][0]) / denom;
    beta[1] = (-v[2][1] * col[0][0] + v[2][0] * col[0][1]) / denom;
    beta[2] = 1. - beta[0] - beta[1];
    return (beta[0] >= 0 && beta[1] >= 0 && beta[2] >= 0);
}

inline double min_q_dist(double q, const std::vector<rational_type>& r)
{
    double _min = std::numeric_limits<double>::max();
    for (int i = 0 ; i < r.size() ; ++i) {
        _min = std::min(fabs(q - spurt::value(r[i])), _min);
    }
    
    return _min;
}

template<typename T, typename F>
struct Lt_point_norm {
    typedef std::pair<T, F> pair_type;
    
    bool operator()(const pair_type& p0, const pair_type& p1) {
        return p0.second < p1.second;
    }
    
};

struct vector_func {
    vector_func(unsigned int period, double minq, double maxq)
        : __period(period), __minq(minq), __maxq(maxq) {}
        
    int order() const {
        return 1;
    }
    
    std::string name() const {
        std::ostringstream os;
        os << __period << "-vector";
        return os.str();
    }
    
    bool is_valid(const spurt::point_data& d) const {
        if (d.period() < __minq || d.period() > __maxq) {
            return false;
        }
        // else {
        //  std::cerr << "period " << d.period() << " is between " << __minq << " and " << __maxq << '\n';
        // }
        
        nvis::vec2 v = spurt::vector_value(d, __period, spurt::__default_metric);
        return (v[0] != std::numeric_limits<double>::max());
    }
    
    std::string value_string(const spurt::point_data& d) const {
        if (d.period() < __minq && d.period() > __maxq) {
            return "0 0 0";
        }
        
        nvis::vec2 v = spurt::vector_value(d, __period, spurt::__default_metric);
        std::ostringstream os;
        os << v[0] << " " << v[1] << " 0";
        return os.str();
    }
    
    unsigned int __period;
    double __minq, __maxq;
};

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(file), std::string(ts));
    field->periodic_coordinates(false);
    
    bool per[2] = {true, false};
    spurt::map_metric metric(field->bounds(), per);
    
    spurt::__default_metric = metric;
    
    poincare_map pmap(field);
    pmap.precision(h);
    orbit_integrator<poincare_map> intg(pmap, m, metric);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double width = diagonal[0];
    double height = diagonal[1];
    double x_center = bounds.center()[0];
    
    std::set<rational_type> tmp_factors;
    for (int i = 1 ; i < max_period ; ++i) {
        for (int j = 1 ; j <= i ; ++j) {
            tmp_factors.insert(rational_type(i, j));
        }
    }
    std::vector<double> factors;
    std::cerr << "interesting rational periods = \n";
    for (std::set<rational_type>::const_iterator it = tmp_factors.begin() ; it != tmp_factors.end() ; ++it) {
        factors.push_back(spurt::value(*it));
        std::cerr << factors.back() << "\n";
    }
    double min_dq = std::numeric_limits<double>::max();
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    srand48(time(0));
    std::vector< std::pair<point_type, spurt::point_data> > all_points;
    
    // reset central orbit repository
    spurt::__map_orbits.clear();
    
#if 0
    int nsamples[] = { 2*max_period, 100 };
    nvis::vec2 dp(1, 1);
    spurt::sample_on_raster(bounds, pmap, m, nsamples, dp);
#else
    {
        unsigned int res[2];
        double ratio = bounds.size()[1] / bounds.size()[0];
        res[0] = ceil(sqrt((double)n1 / ratio));
        res[1] = ceil(res[0] * ratio);
        spurt::experimental::sample_on_raster(spurt::__map_orbits,
                                               pmap, metric, res, m);
        std::cerr << "Done\n";
    }
#endif
    
    // assign period to each orbit and enter data points in array
    for (int i = 0 ; i < spurt::__map_orbits.size() ; ++i) {
        spurt::orbit& obt = spurt::__map_orbits[i];
        double q = (spurt::period_x_periodic(obt.points(), metric)).first;
        for (int j = 0 ; j < obt.size() ; ++j) {
            obt[j] = metric.modulo(obt[j]);
        }
        obt.period() = q;
        if (q < 0) {
            std::cerr << "Returned a negative period!\n";
            std::cerr << "points were:\n";
            for (int j = 0 ; j < obt.size() ; ++j) {
                std::cerr << j << ": " << obt[j] << '\n';
            }
            assert(false);
        }
        
        for (int j = 0 ; j < obt.points().size() ; ++j) {
            all_points.push_back(std::make_pair(obt.points()[j], spurt::point_data(i, j)));
        }
    }
    
    
    
    std::cerr << "building initial triangulation after sampling\n";
    // form the boundary
    double ref_l = nvis::norm(bounds.size());
    nvis::vec2 dx(0.005*ref_l, 0.005*ref_l);
    std::pair<nvis::vec2, spurt::point_data> boundary[4];
    std::vector<nvis::vec2> corner(4);
    corner[0] = bounds.min() - dx;
    corner[1] = bounds.min() + nvis::vec2(bounds.size()[0] + dx[0], -dx[1]);
    corner[2] = bounds.max() + dx;
    corner[3] = bounds.min() + nvis::vec2(-dx[0], bounds.size()[1] + dx[1]);
    size_t orbit_id = spurt::__map_orbits.size();
    spurt::__map_orbits.push_back(spurt::orbit(corner, 0));
    for (int i = 0 ; i < 4 ; ++i) {
        boundary[i].first = corner[i];
        boundary[i].second = spurt::point_data(orbit_id, i);
    }
    
    sampling_mesh = new mesh_type(boundary, spurt::point_locator());
    mesh_type& mesh = *sampling_mesh;
    try {
        mesh.insert_points(all_points);
    } catch (...) {
        export_VTK(mesh, "bug.vtk", "this is a faulty triangulation", true, false);
        return -1;
    }
    
    {
        std::ostringstream os_name, os_comment;
        os_name << outs << "-initial-mesh.vtk";
        os_comment << "Sample points from " << file << " at t=" << ts << " with estimated safety factor";
        spurt::export_VTK(mesh, os_name.str(), os_comment.str(), true);
    }
    
    
    bool valid_mesh;
    nvis::bbox2 interior_bounds(bounds.min() + 0.01*diagonal,
                                bounds.min() + 0.99*diagonal);
                                
    std::cerr << "before refinement, mesh statistics are:\n";
    spurt::value_stats<double> area_stats, edge_stats;
    spurt::statistics(mesh, area_stats, edge_stats, interior_bounds);
    std::cerr << "area stats: mean=" << area_stats.mean << ", variance=" << area_stats.var
              << ", min=" << area_stats.min << ", max=" << area_stats.max << std::endl;
    std::cerr << "edge stats: mean=" << edge_stats.mean << ", variance=" << edge_stats.var
              << ", min=" << edge_stats.min << ", max=" << edge_stats.max << std::endl;
              
    // std::cerr << "refinement to achieve prescribed max area is... " << std::flush;
    // valid_mesh = spurt::refine(mesh, intg,
    //                             spurt::triangle_area_controller(max_area),
    //                             max_nb_triangles);
    std::cerr << "refinement to achieve prescribed max edge length is... " << std::flush;
    valid_mesh = spurt::refine(mesh, intg,
                                spurt::triangle_edge_length_controller(edge_stats.mean + sqrt(edge_stats.var)),
                                max_nb_triangles);
    std::cerr << (valid_mesh ? "successful" : "not successful") << " (" << mesh.get_nb_vertices()
              << " / " << mesh.get_nb_triangles() << ")\n";
              
    {
        std::ostringstream os_name, os_comment;
        os_name << outs << "-area-control-only.vtk";
        os_comment << "Sample points from " << file << " at t=" << ts << " with estimated safety factor";
        spurt::export_VTK(mesh, os_name.str(), os_comment.str(), true);
    }
    typedef boost::rational<int>                rational_type;
    typedef spurt::point_data                  data_type;
    typedef std::pair<nvis::vec2, data_type>    data_point_type;
    
    spurt::MarkedPos marked_positions(0.5);
    std::vector< std::vector< spurt::fixpoint > > chains;
    all_chains.clear();
    
    spurt::map_debug::verbose_level = 1;
    
    for (int p = 1 ; p <= max_period ; ++p) {
        std::cerr << "\nprocessing period " << p << "...\n";
        
        std::cerr << "there are currently " << spurt::__map_orbits.size() << " orbits on record\n";
        
        // determine what rationals should be considered
        std::set<rational_type> rationals;
        for (int den = 1 ; den <= p ; ++den) {
            rational_type r(p, den);
            if (r.numerator() != p) {
                continue;
            }
            rationals.insert(r);
        }
        ref_values.clear();
        std::copy(rationals.begin(), rationals.end(),
                  std::back_inserter(ref_values));
        std::cerr << "there are " << ref_values.size() << " interesting values\n";
        for (int i = 0 ; i < ref_values.size() ; ++i) {
            std::cerr << ref_values[i] << " = " << spurt::value(ref_values[i]) << '\n';
        }
        
        // refine triangulation around interesting rationals
        for (std::set<rational_type>::const_iterator iter = rationals.begin() ; iter != rationals.end() ; ++iter) {
            double r = spurt::value(*iter);
            std::cerr << "refinement around q=" << r << " was... " << std::flush;
            interval_type range((p == 1 ? 0 : r - dq), r + dq);
            valid_mesh = spurt::refine(mesh, intg,
                                        spurt::triangle_interval_inclusion_controller(range, min_area),
                                        max_nb_triangles);
            std::cerr << (valid_mesh ? "successful" : "not successful") << '\n';
            std::cout << mesh.get_nb_vertices() << " points sampled after value-driven 2D refinement\n";
            
            // save the resulting mesh as a vector field
            std::ostringstream os;
            os << outs << "-vector_field-p=" << p << "-q=" << r << ".vtk";
            vector_func functor(p, range.__min, range.__max);
            std::ostringstream os2;
            os2 << p << "-vector for q between " << range.__min << " and " << range.__max;
            spurt::export_VTK(mesh, os.str(), os2.str(), functor, true);
            
        } // for all rational safety factors matching the prescribed period
        
        // loop over all current vertices
        
        
        typedef std::pair<unsigned int, unsigned int>   pair_ui;
        std::vector<std::pair<pair_ui, double> > point_to_norm;
        
#ifdef I_LOVE_TO_BE_INEFFICIENT
        std::set<size_t> close_orbits;
        for (int i = 0 ; i < mesh.get_nb_vertices() ; ++i) {
            const data_type& data = mesh.get_data(i);
            if (close_orbits.find(data.orbit_id()) != close_orbits.end()) {
                continue;
            }
            
            // check that we are close to valid safety factor
            if (min_q_dist(data.period(), ref_values) > _dq) {
                continue;
            }
            close_orbits.insert(data.orbit_id());
        }
        std::cerr << close_orbits.size() << " candidates within prescribed q-distance for period "
                  << p << std::endl;
                  
        // identify position with smallest p-norm for all found chains      std::set<size_t>::const_iterator iter;
        for (iter = close_orbits.begin() ; iter != close_orbits.end() ; ++iter) {
            int orbit_id = *iter;
            size_t sz = spurt::__map_orbits[orbit_id].size();
            double min_norm = std::numeric_limits<double>::max();
            int min_id = -1;
            for (int i = 0 ; i < sz ; ++i) {
                double norm = nvis::norm(vector_value(orbit_id, i, p, metric));
                if (norm < min_norm) {
                    min_id = i;
                    min_norm = norm;
                }
            }
            
            point_to_norm.push_back(std::make_pair(pair_ui(orbit_id, min_id), min_norm));
        }
#else
        double max_norm = 0;
        int orbit_count = 0;
        for (int i = 0 ; i < spurt::__map_orbits.size() ; ++i) {
            const spurt::orbit& obt = spurt::__map_orbits[i];
            if (min_q_dist(obt.period(), ref_values) < 0.1) { // 0.1 as in "way out of whack"
                ++orbit_count;
                unsigned int pid;
                double norm = tokatopo::min_norm(pid, i, metric, p);
        
                std::cerr << "orbit #" << i << " has min norm at " << obt[pid] << " = " << norm << std::endl;
        
                if (norm < close_enough) {
                    point_to_norm.push_back(std::make_pair(pair_ui(i, pid), norm));
                    max_norm = std::max(max_norm, norm);
                }
            }
        }
        std::cerr << "selected " << point_to_norm.size() << " points from "
                  << orbit_count << " orbits lying in 0.1 q-vicinity. "
                  << "Max norm of smallest step is " << max_norm << '\n';
#endif
        
        std::sort(point_to_norm.begin(), point_to_norm.end(),
                  Lt_point_norm<pair_ui, double>());
                  
        std::vector<nvis::vec2> seeds(point_to_norm.size());
        for (int i = 0 ; i < point_to_norm.size() ; ++i) {
            pair_ui pui = point_to_norm[i].first;
            unsigned int orbit_id = pui.first;
            unsigned int p_id = pui.second;
            const spurt::orbit& orbit = spurt::__map_orbits[orbit_id];
            seeds[i] = orbit[p_id];
            probed_seeds.push_back(std::pair<unsigned int, nvis::vec2>(orbit_id, seeds[i]));
        }
        
        pmap.precision(h / (4.*(double)p));
        find_fixed_points(pmap, marked_positions, seeds, p,
                          max_period, 0.5, (double)p*10.*h, chains,
                          tokatopo::verbose_level);
        for (int n = 0 ; n < chains.size() ; ++n) {
            for (int k = 0 ; k < chains[n].size() ; ++k) {
                marked_positions.add(chains[n][k].pos, chains[n][k].K);
            }
        }
        std::cerr << "for period p=" << p << " fixed point search returned " << chains.size() << " chains ";
        std::copy(chains.begin(), chains.end(), std::back_inserter(all_chains));
        std::cerr << "and there are a total of " << all_chains.size() << " chains so far\n";
        //
        // std::cerr << "minimum per-orbit " << p << "-norms\n";
        // for (int i = 0 ; i < std::min(point_to_norm.size(), (size_t)10) ; ++i) {
        //  std::cerr << i + 1 << ": " << point_to_norm[i].first << " has norm " << point_to_norm[i].second
        //            << '\n';
        // }
        
        std::ostringstream os_name, os_comment;
        os_name << outs << "-refined_for_up_to_p=" << p << ".vtk";
        os_comment << "Sample points from " << file << " at t=" << ts << " with estimated safety factor";
        spurt::export_VTK(mesh, os_name.str(), os_comment.str(), true);
        
        std::cerr << "upon completing fixed points search there are " << spurt::__map_orbits.size()
                  << " orbits available\n";
                  
    } // for all periods
}





















































