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


#define __REFINE_SAFETY_FACTOR__

namespace tokatopo {
int res, maxp, maxit, n1, n2, it, m, maxround, max_nb_triangles, max_period, period;
char* outs, *file, *ts;
double h, dq, min_area, max_area, max_ratio;

double close_enough = 0;

mesh_type* sampling_mesh;
std::vector< std::vector< xavier::fixpoint > > all_chains;
xavier::map_metric metric;
unsigned int verbose_level = 1;
std::vector<std::pair<unsigned int, nvis::vec2> > probed_seeds;
}

using namespace tokatopo;


void tokatopo::initialize(int argc, char* argv[])
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
    hestOptAdd(&hopt, "dq",     "safety factor step",   airTypeDouble,  0, 1, &dq,                  "0.01",     "maximum discrepancy between consecutive safety factors");
    hestOptAdd(&hopt, "mr",     "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,           "2.",       "maximum triangle aspect ratio");
    // hestOptAdd(&hopt, "nr",  "max # rounds",         airTypeInt,     0, 1, &maxround,            "10",       "max number of 1D approximation iterations");
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

inline double min_q_dist(double q, const std::vector<rational_type>& r)
{
    double _min = std::numeric_limits<double>::max();
    for (int i = 0 ; i < r.size() ; ++i) {
        _min = std::min(fabs(q - xavier::value(r[i])), _min);
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

void shuffle(std::vector<size_t>& v)
{
    int N = v.size();
    std::vector<size_t> tmp;
    srand48(time(0));
    std::vector<bool> used(N, false);
    for (int n = 0 ; n < N ; ++n) {
        int i = (int)floor(drand48() * N);
        while (used[i]) {
            i = (i + 1) % N;
        }
        tmp.push_back(v[i]);
        used[i] = true;
    }
    v.swap(tmp);
}

void test_error_estimate(tokamak_nimrod_parametric* field)
{
    std::fstream output("test_error_estimate.txt", std::ios::out);
    poincare_map pmap(field);
    nvis::bbox2 bounds = field->bounds();
    bounds.min() += nvis::vec2(1, 1);
    bounds.max() -= nvis::vec2(1, 1);
    nvis::vec2 diagonal = bounds.size();
    
    std::vector<nvis::vec2> seed(1000);
    for (int i = 0 ; i < 1000 ; ++i) {
        nvis::vec2 x = bounds.min();
        x[0] += drand48() * diagonal[0];
        x[1] += drand48() * diagonal[1];
        seed[i] = x;
    }
    
    for (double eps = 1.0e-10 ; eps <= 1.0e-2 ; eps *= 2.0) {
        pmap.precision(eps);
        double mean = 0;
        double meann = 0;
        std::vector<double> error, nsteps;
        map2d::value_type val;
        for (int i = 0 ; i < seed.size() ; ++i) {
            try {
                val = pmap.map_complete(seed[i], 1);
            } catch (...) {
                continue;
            }
            
            double err = nvis::norm(val.err);
            
            if (std::isnan(err) || std::isinf(err) || err > 100. || err == 0.) {
                std::cerr << "invalid error value for eps=" << eps << ": " << err << " at " << seed[i]
                          << " after " << val.nsteps << " successful steps and "
                          << val.total_nsteps << " total steps" << std::endl;
                continue;
            }
            mean += err;
            error.push_back(err);
            
            // std::cerr << "error(" << eps << ") at " << seed[i] << " = " << err << " after " << val.nsteps << " steps\n";
            
            double n = val.nsteps;
            meann += n;
            nsteps.push_back(n);
        }
        mean /= (double)error.size();
        meann /= (double)nsteps.size();
        
        double var = 0;
        double varn = 0;
        for (int i = 0 ; i < error.size() ; ++i) {
            var += (error[i] - mean) * (error[i] - mean);
            varn += (nsteps[i] - meann) * (nsteps[i] - meann);
        }
        var /= error.size();
        varn /= nsteps.size();
        
        double min = *std::min_element(error.begin(), error.end());
        double minn = *std::min_element(nsteps.begin(), nsteps.end());
        double max = *std::max_element(error.begin(), error.end());
        double maxn = *std::max_element(nsteps.begin(), nsteps.end());
        
        output << eps << " \t" << mean << " \t" << var << " \t" << min << " \t" << max << " \t"
               << meann << " \t" << varn << " \t" << minn << " \t" << maxn << '\n';
    }
    output.close();
}

int tokatopo::compute_fixed_points()
{
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(file), std::string(ts));
    field->periodic_coordinates(false);
    
    probed_seeds.clear();
    
    // test_error_estimate(field);
    
    bool per[2] = {true, false};
    metric = xavier::map_metric(field->bounds(), per);
    
    xavier::__default_metric = metric;
    
    poincare_map pmap(field);
    pmap.precision(h);
    xavier::orbit_integrator<poincare_map> intg(pmap, m, metric);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double width = diagonal[0];
    double height = diagonal[1];
    double x_center = bounds.center()[0];
    
    if (close_enough <= 0) {
        close_enough = 0.01 * nvis::norm(diagonal);    // 1% is _very_ conservative
    }
    
    std::set<rational_type> tmp_factors;
    for (int i = 1 ; i < max_period ; ++i) {
        for (int j = 1 ; j <= i ; ++j) {
            tmp_factors.insert(rational_type(i, j));
        }
    }
    std::vector<double> factors;
    std::cerr << "interesting rational periods = \n";
    for (std::set<rational_type>::const_iterator it = tmp_factors.begin() ; it != tmp_factors.end() ; ++it) {
        factors.push_back(xavier::value(*it));
        std::cerr << factors.back() << "\n";
    }
    double min_dq = std::numeric_limits<double>::max();
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    srand48(time(0));
    std::vector< std::pair<point_type, xavier::point_data> > all_points;
    
    // reset central orbit repository
    xavier::__map_orbits.clear();
    
    int nsamples[] = { 2*max_period, 100 };
    nvis::vec2 dp(1, 1);
    nvis::timer intg_timer;
    xavier::sample_on_raster(bounds, pmap, m, nsamples, dp);
    std::cerr << "integration of " << nsamples[0]*nsamples[1] << " orbits over " << m << " iterations took " << intg_timer.elapsed() << " s.\n";
    
    // assign period to each orbit and enter data points in array
    for (int i = 0 ; i < xavier::__map_orbits.size() ; ++i) {
        xavier::orbit& obt = xavier::__map_orbits[i];
        double q = (xavier::period_x_periodic(obt.points(), metric)).first;
        for (int j = 0 ; j < obt.size() ; ++j) {
            obt[j] = metric.modulo(obt[j]);
        }
        obt.period() = q;
        if (q < 0) {
            std::cerr << "Returned a negative period!\n";
            assert(false);
        }
        
        for (int j = 0 ; j < obt.points().size() ; ++j) {
            all_points.push_back(std::make_pair(obt.points()[j], xavier::point_data(i, j)));
        }
    }
    
    std::cerr << "building initial triangulation after sampling\n";
    // form the boundary
    double ref_l = nvis::norm(bounds.size());
    nvis::vec2 dx(0.005*ref_l, 0.005*ref_l);
    std::pair<nvis::vec2, xavier::point_data> boundary[4];
    std::vector<nvis::vec2> corner(4);
    corner[0] = bounds.min() - dx;
    corner[1] = bounds.min() + nvis::vec2(bounds.size()[0] + dx[0], -dx[1]);
    corner[2] = bounds.max() + dx;
    corner[3] = bounds.min() + nvis::vec2(-dx[0], bounds.size()[1] + dx[1]);
    size_t orbit_id = xavier::__map_orbits.size();
    xavier::__map_orbits.push_back(xavier::orbit(corner, 0));
    for (int i = 0 ; i < 4 ; ++i) {
        boundary[i].first = corner[i];
        boundary[i].second = xavier::point_data(orbit_id, i);
    }
    
    sampling_mesh = new mesh_type(boundary, xavier::point_locator());
    mesh_type& mesh = *sampling_mesh;
    nvis::timer triangulation_time;
    try {
        mesh.insert_points(all_points);
    } catch (...) {
        export_VTK(mesh, "bug.vtk", "this is a faulty triangulation", true, false);
        return -1;
    }
    std::cerr << "triangulation of " << all_points.size() << " initial points took " << triangulation_time.elapsed() << " s.\n";
    //
    //
    //
    // // ensure a minimal sampling density of the map
    //
    // int res_y = n1 / 2;
    // int res_x = (int)floor(width / height * (double)res_y);
    //
    // std::cerr << "resolution of density control mesh is " << res_x
    //           << " x " << res_y << std::endl;
    // xavier::regular_mesh raster(bounds, res_x, res_y);
    //
    // // ** shuffle seeding locations to avoid degenerate triangulation
    //
    // std::vector<bool> valid_id(res_x*res_y, true);
    // for (int i = 0 ; i < all_points.size() ; ++i) {
    //  int n = raster.idx(all_points[i].first);
    //  if (n < 0) continue;
    //  else valid_id[n] = false;
    // }
    // std::vector<size_t> seed_ids;
    // for (int i = 0 ; i < valid_id.size() ; ++i) {
    //  if (valid_id[i]) {
    //      seed_ids.push_back(i);
    //  }
    // }
    // shuffle(seed_ids);
    //
    // int nb_seeds = seed_ids.size();
    // std::cerr << nb_seeds << " cells out of " << res_x * res_y << " are still empty\n";
    // for (int i = 0 ; i < seed_ids.size() ; ++i) {
    //  int n = seed_ids[i];
    //  if (!valid_id[n]) continue;
    //  nvis::vec2 x0 = raster.random_point(n);
    //  std::vector<nvis::vec2> steps;
    //  std::vector<xavier::point_data> data;
    //  intg(x0, steps, data);
    //  valid_id[n] = false;
    //
    //  std::vector<int> modified;
    //  mesh.insert_points(steps, data, modified);
    //
    //  for (int i = 0 ; i < steps.size() ; ++i) {
    //      valid_id[raster.idx(metric.modulo(steps[i]))] = false;
    //  }
    // }
    // std::cerr << mesh.get_nb_vertices() << " points after density control\n";
    
    {
        std::ostringstream os_name, os_comment;
        os_name << outs << "-initial-mesh.vtk";
        os_comment << "Sample points from " << file << " at t=" << ts << " with estimated safety factor";
        xavier::export_VTK(mesh, os_name.str(), os_comment.str(), true);
    }
    
    
    bool valid_mesh;
    nvis::bbox2 interior_bounds(bounds.min() + 0.01*diagonal,
                                bounds.min() + 0.99*diagonal);
                                
    std::cerr << "before refinement, mesh statistics are:\n";
    xavier::value_stats<double> area_stats, edge_stats;
    xavier::statistics(mesh, area_stats, edge_stats, interior_bounds);
    std::cerr << "area stats: mean=" << area_stats.mean << ", variance=" << area_stats.var
              << ", min=" << area_stats.min << ", max=" << area_stats.max << std::endl;
    std::cerr << "edge stats: mean=" << edge_stats.mean << ", variance=" << edge_stats.var
              << ", min=" << edge_stats.min << ", max=" << edge_stats.max << std::endl;
              
    // std::cerr << "refinement to achieve prescribed max area is... " << std::flush;
    // valid_mesh = xavier::refine(mesh, intg,
    //                             xavier::triangle_area_controller(max_area),
    //                             max_nb_triangles);
    std::cerr << "refinement to achieve prescribed max edge length is... " << std::flush;
    valid_mesh = xavier::refine(mesh, intg,
                                xavier::triangle_edge_length_controller(edge_stats.mean + sqrt(edge_stats.var)),
                                max_nb_triangles);
    std::cerr << (valid_mesh ? "successful" : "not successful") << " (" << mesh.get_nb_vertices()
              << " / " << mesh.get_nb_triangles() << ")\n";
              
    {
        std::ostringstream os_name, os_comment;
        os_name << outs << "-area-control-only.vtk";
        os_comment << "Sample points from " << file << " at t=" << ts << " with estimated safety factor";
        xavier::export_VTK(mesh, os_name.str(), os_comment.str(), false);
    }
    typedef boost::rational<int>                rational_type;
    typedef xavier::point_data                  data_type;
    typedef std::pair<nvis::vec2, data_type>    data_point_type;
    
    xavier::MarkedPos marked_positions(0.5);
    std::vector< std::vector< xavier::fixpoint > > chains;
    all_chains.clear();
    
    xavier::map_debug::verbose_level = 1;
    
    for (int p = 1 ; p <= max_period ; ++p) {
        std::cerr << "\nprocessing period " << p << "...\n";
        
        std::cerr << "there are currently " << xavier::__map_orbits.size() << " orbits on record\n";
        
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
            std::cerr << ref_values[i] << " = " << xavier::value(ref_values[i]) << '\n';
        }
        
        // refine triangulation around interesting rationals
        double _dq;
        for (std::set<rational_type>::const_iterator iter = rationals.begin() ; iter != rationals.end() ; ++iter) {
            double r = xavier::value(*iter);
            std::cerr << "refinement around q=" << r << " was... " << std::flush;
            _dq = (p == 1 ? 0.05 : dq);
            interval_type range(r - _dq, r + _dq);
            valid_mesh = xavier::refine(mesh, intg,
                                        xavier::triangle_interval_inclusion_controller(range, min_area),
                                        max_nb_triangles);
            std::cerr << (valid_mesh ? "successful" : "not successful") << '\n';
            // std::cout << mesh.get_nb_vertices() << " points sampled after value-driven 2D refinement\n";
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
            size_t sz = xavier::__map_orbits[orbit_id].size();
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
        for (int i = 0 ; i < xavier::__map_orbits.size() ; ++i) {
            const xavier::orbit& obt = xavier::__map_orbits[i];
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
            const xavier::orbit& orbit = xavier::__map_orbits[orbit_id];
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
        xavier::export_VTK(mesh, os_name.str(), os_comment.str(), true);
        
        std::cerr << "upon completing fixed points search there are " << xavier::__map_orbits.size()
                  << " orbits available\n";
                  
    } // for all periods
}
































