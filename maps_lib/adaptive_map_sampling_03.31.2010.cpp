#include "adaptive_map_sampling.hpp"

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

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>
#include <math/bounding_box.hpp>

// xavier's math stuff
#include <math/math.hpp>

// poincare lib
#include <poincare/metric.hpp>
#include <poincare/newton.hpp>
#include <poincare/topology.hpp>
#include <poincare/fixpoints.hpp>
#include <poincare/macros.hpp>
#include <poincare/map.hpp>
#include <poincare/newton.hpp>
#include <poincare/basic_math.hpp>

// map analysis
#include "sampling.hpp"
#include "experimental.hpp"
#include "orbits.hpp"
#include "period.hpp"
#include "regular_mesh.hpp"
#include "approximation.hpp"
#include "index.hpp"
#include "find_fixed_points.hpp"

// back to the deads
// #include "separatrix.hpp"

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
#include "definitions.hpp"
#include "triangulation.hpp"
#include "adaptive_triangulation.hpp"
#include "quality_control.hpp"
#include "IO.hpp"

using namespace map_analysis;

mesh_type* sampling_mesh;
std::vector< std::vector< xavier::fixpoint > > all_chains;

std::vector<rational_type> ref_values;

typedef xavier::point_locator<double, int, 2>       point_locator_type;

template<typename T, typename F>
struct Lt_point_norm {
    typedef std::pair<T, F> pair_type;
    
    bool operator()(const pair_type& p0, const pair_type& p1) {
        return p0.second < p1.second;
    }
    
};

void filter_orbits(std::vector<unsigned int>& ids,
                   const std::vector<orbit>& orbits,
                   double& min, double& max)
{
    ids.clear();
    if (!orbits.size()) {
        return;
    }
    double minq = min, maxq = max;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();
    for (int i = 0 ; i < orbits.size() ; ++i) {
        double q = orbits[i].period();
        if (q >= minq && q <= maxq) {
            ids.push_back(i);
            if (q < min) {
                min = q;
            }
            if (q > max) {
                max = q;
            }
        }
    }
}

inline void orbits_to_points(std::vector<nvis::vec2>& points,
                             std::vector<point_data>& data,
                             const std::vector<orbit>& orbits,
                             const std::vector<unsigned int>& valid_ids)
{
    points.clear();
    data.clear();
    for (unsigned int i = 0 ; i < valid_ids.size() ; ++i) {
        const orbit& orb = orbits[valid_ids[i]];
        for (unsigned int j = 0 ; j < orb.size() ; ++j) {
            points.push_back(orb[j]);
            data.push_back(point_data(valid_ids[i], j));
        }
    }
}

inline void orbits_to_points(std::vector<nvis::vec2>& points,
                             std::vector<point_data>& data,
                             const std::vector<orbit>& orbits)
{
    points.clear();
    data.clear();
    for (unsigned int i = 0 ; i < orbits.size() ; ++i) {
        const orbit& orb = orbits[i];
        for (unsigned int j = 0 ; j < orb.size() ; ++j) {
            points.push_back(orb[j]);
            data.push_back(point_data(i, j));
        }
    }
}

void map_analysis::
adaptive_map_sampling(adaptive_map_sampling_output& output,
                      const adaptive_map_sampling_params& params)
{
    typedef nvis::vec2                                  point_type;
    typedef point_data                                  data_type;
    typedef std::pair<nvis::vec2, data_type>            data_point_type;
    typedef poincare_index<mesh_type, poincare_map>     poincare_index_type;
    typedef mesh_type::triangle_type                    triangle_type;
    typedef xavier::Edge                                edge_index_type;
    
    output.p_orbits.clear();
    output.p_meshes.clear();
    output.p_cand_tris.clear();
    output.p_sing_tris.clear();
    output.p_prob_tris.clear();
    output.p_chains.clear();
    output.p_prob_edges.clear();
    output.p_degenerate_points.clear();
    output.p_separatrices.clear();
    output.p_rejected_separatrices.clear();
    //
    output.p_orbits.resize(params.mp);
    output.p_meshes.resize(params.mp);
    output.p_cand_tris.resize(params.mp);
    output.p_sing_tris.resize(params.mp);
    output.p_prob_tris.resize(params.mp);
    output.p_index_vectors.resize(params.mp);
    output.p_chains.resize(params.mp);
    output.p_prob_edges.resize(params.mp);
    output.p_degenerate_points.resize(params.mp);
    output.p_separatrices.resize(params.mp);
    output.p_rejected_separatrices.resize(params.mp);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(params.hdf_name),
                                          params.time_step);
    field->periodic_coordinates(false);
    
    bool per[2] = {true, false};
    metric_type metric(field->bounds(), per);
    
    static_data::metric = metric;
    
    double h = params.eps;
    
    poincare_map pmap(field);
    pmap.precision(params.eps);
    orbit_integrator<poincare_map> intg(pmap, params.m, metric);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double tiny_length = params.dx;
    double dtheta = params.dtheta;
    nvis::bbox2 inflated_bounds(bounds);
    inflated_bounds.min() -= 0.005 * diagonal;
    inflated_bounds.max() += 0.005 * diagonal;
    std::cerr << "inflated bounds are " << inflated_bounds << '\n';
    
    std::cerr << "approximation error bound is " << params.err << '\n';
    
    size_t nb_threads = 1;
#if _OPENMP
    nb_threads = omp_get_max_threads();
    std::cout << nb_threads << " threads available\n";
#endif
    
    srand48(time(0));
    std::vector< std::pair<point_type, point_data> > all_points;
    
    // reset central orbit repository
    static_data::central_map_orbits.clear();
    
    // initialize base triangulation
    std::pair<nvis::vec2, point_data> boundary[4];
    std::vector<nvis::vec2> corner(4);
    corner[0] = inflated_bounds.min();
    corner[2] = inflated_bounds.max();
    corner[1] = nvis::vec2(corner[2][0], corner[0][1]);
    corner[3] = nvis::vec2(corner[0][0], corner[2][1]);
    size_t frame_orbit_id = static_data::central_map_orbits.size();
    static_data::central_map_orbits.push_back(orbit(corner, 0));
    for (int i = 0 ; i < 4 ; ++i) {
        boundary[i].first = corner[i];
        boundary[i].second = point_data(frame_orbit_id, i);
    }
    output.base_mesh = mesh_type(boundary, point_locator_type());
    mesh_type& base_mesh = output.base_mesh;
    
    unsigned int res[2];
    double ratio = bounds.size()[1] / bounds.size()[0];
    res[0] = ceil(sqrt((double)params.n / ratio));
    res[1] = ceil(res[0] * ratio);
    nvis::timer _timer;
    experimental::sample_on_raster(static_data::central_map_orbits,
                                   pmap, metric, res, params.m);
    std::cerr << "initial sampling of the map at resolution "
              << res[0] << " x " << res[1] << " completed in "
              << _timer.elapsed() << " s.\n";
              
    // assign period to each orbit
    _timer.restart();
    for (int i = 0 ; i < static_data::central_map_orbits.size() ; ++i) {
        orbit& obt = static_data::central_map_orbits[i];
        // double q = (xavier::period_x_periodic(obt.points(), metric)).first;
        double q = dist_based_x_period(obt.points(), metric);
        for (int j = 0 ; j < obt.size() ; ++j) {
            obt[j] = metric.modulo(obt[j]);
        }
        obt.period() = q;
        if (q < 0) {
            std::cerr << "Returned a negative period!\n";
            std::cerr << "points were: \n";
            for (int j = 0 ; j < obt.size() ; ++j) {
                std::cerr << j << ": " << obt[j] << '\n';
            }
            assert(false);
        }
    }
    std::cerr << "period computation of initial orbits completed in " << _timer.elapsed() << " s.\n";
    
    {
        std::vector<nvis::vec2> points;
        std::vector<point_data> data;
        orbits_to_points(points, data, static_data::central_map_orbits);
        _timer.restart();
        base_mesh.insert_points(points, data);
        std::cerr << "meshing of initial orbits completed in " << _timer.elapsed() << " s.\n";
        
#ifdef EXPORT_TO_VTK
        std::ostringstream os;
        os << params.out_name << "-base_mesh.vtk";
        export_VTK(base_mesh, os.str(), "N / A", false);
#endif
    }
    base_mesh.set_tolerance(1.0e-7);
    
    xavier::map_debug::verbose_level = 1;
    
    typedef newton::cached_map<poincare_map>    cached_map_type;
    typedef newton::rhs_map<cached_map_type>    rhs_type;
    typedef newton::J_map<cached_map_type>      jacobian_type;
    
    /*************************************************************************
    *
    *           Main loop over all relevant periods of the map
    *
    *************************************************************************/
    for (int period = params.pm ; period <= params.mp ; ++period) {
        _timer.restart();
        
        std::cerr << "\nprocessing period " << period << "...\n";
        
        std::cerr << "there are currently " << static_data::central_map_orbits.size()
                  << " orbits on record and " << base_mesh.get_nb_triangles() << " points in base triangulation\n";
                  
        // this is a hack that I would only disclose under torture...
        double h_p = (period == 1 ? h : (0.25 * h / (double)(period - 1)));
        // double h_p = h;
        
        intg.precision(h_p);
        std::cerr << "integrator precision is currently " << h_p << '\n';
        
        output.p_meshes[period-1] = base_mesh; // copy mesh
        mesh_type& local_mesh = output.p_meshes[period-1];
        
        
        // determine what rationals should be considered
        
        std::set<rational_type> rationals;
        for (int den = 1 ; den <= period ; ++den) {
            rational_type q(period, den);
            if (q.numerator() != period) {
                continue;    // period and den are not mutually prime
            }
            rationals.insert(q);
        }
        
        std::vector<double> valid_qs;
        for (std::set<rational_type>::iterator it = rationals.begin() ; it != rationals.end() ; ++it) {
            valid_qs.push_back(xavier::value(*it));
        }
        
        std::cerr << "there are " << valid_qs.size() << " interesting rational periods:\n";
        for (int i = 0 ; i < valid_qs.size() ; ++i) {
            std::cerr << valid_qs[i] << "\t";
        }
        std::cerr << '\n';
        
        
        // refine mesh in vicinity of selected rational periods
        
        nvis::timer sub_timer;
        std::cerr << "refinement to achieve prescribed max angular variation\n";
        angular_variation_priority priority(period, params.err, metric, params.dq);
        bool ok = refine(local_mesh, intg, priority, params.mt);
        std::cerr << "refinement of " << period << "-mesh completed in " << sub_timer.elapsed() << " s.\n";
        
        std::cerr << "after refinement mesh contains " << local_mesh.get_nb_triangles()
                  << " triangles\n";
                  
        sub_timer.restart();
        triangle_filter filter(period, params.err, metric, params.dq);
        std::vector<unsigned int> included_triangles;
        for (int i = 0 ; i < local_mesh.get_nb_triangles() ; ++i) {
            mesh_type::point_type p[3];
            mesh_type::data_type d[3];
            local_mesh.get_triangle_info(p, d, i);
            if (filter.valid(p, d) > 0) {
                included_triangles.push_back(i);
            }
        }
        std::cerr << "filtering of newly created triangles completed in " << sub_timer.elapsed() << '\n';
        
        std::cerr << "of those refined triangles, " << included_triangles.size()
                  << "meet validity criterion\n";
                  
                  
        poincare_index_type pindex;
        std::vector<unsigned int> singular_triangles;
        
        poincare_map* precise_map = pmap.clone();
        precise_map->precision(0.1*h_p);
        
        
        // determine the set of edges that need to be inspected
        // NB: that could easily be done in parallel but should be fast anyway
        
        sub_timer.restart();
        std::set<edge_index_type> __unique_edges;
        for (int i = 0 ; i < included_triangles.size() ; ++i) {
            size_t id = included_triangles[i];
            mesh_type::point_type p[3];
            mesh_type::data_type d[3];
            local_mesh.get_triangle_info(p, d, id);
            
            // is this triangle remotely interesting?
            if (priority(p, d) < 0) {
                continue;
            }
            std::list<step_type> tmp_steps;
            std::pair<double, double> dpair = pindex.rotation_mag(d, period, metric, tmp_steps);
            if ((dpair.first > 0.9 && dpair.second > 0.45) || dpair.first > 0.98) {
                std::copy(tmp_steps.begin(), tmp_steps.end(), std::back_inserter(output.p_index_vectors[period-1]));
                output.p_cand_tris[period-1].push_back(id);
                const triangle_type& ids = local_mesh.get_triangle_vertices(id);
                __unique_edges.insert(edge_index_type(ids[0], ids[1]));
                __unique_edges.insert(edge_index_type(ids[1], ids[2]));
                __unique_edges.insert(edge_index_type(ids[2], ids[0]));
            }
        }
        std::cerr << "edge filtering completed in " << sub_timer.elapsed() << " s.\n";
        
        std::vector<edge_index_type> unique_edges(__unique_edges.begin(), __unique_edges.end());
        std::cerr << unique_edges.size() << " edges to inspect\n";
        std::cerr << output.p_cand_tris[period-1].size() << " candidate triangles\n";
        
        
        // determine unique edge vertices to avoid computing their p-map twice
        
        std::set<unsigned int> __unique_vertices;
        for (int i = 0 ; i < unique_edges.size() ; ++i) {
            const edge_index_type& eid = unique_edges[i];
            __unique_vertices.insert(eid.i0);
            __unique_vertices.insert(eid.i1);
        }
        std::vector<unsigned int> unique_vertices(__unique_vertices.begin(), __unique_vertices.end());
        std::cerr << "there are " << unique_vertices.size() << " unique vertices to sample\n";
        
        std::vector<std::map<unsigned int, nvis::vec2> > __unique_vertex_pmaps(nb_threads);
        std::vector<std::map<unsigned int, nvis::mat2> > __unique_vertex_jmaps(nb_threads);
        sub_timer.restart();
        
        // do that in parallel
        #pragma omp parallel
        {
            int thread_id = 0;
            #pragma omp for schedule(dynamic,1)
            for (int i = 0 ; i < unique_vertices.size() ; ++i) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                if (!thread_id) {
                    int pct = 100 * i / unique_vertices.size();
                    std::cerr << pct << "% of vertices completed...     \r";
                }
                
                unsigned int vert_id = unique_vertices[i];
                
                poincare_map* my_map = precise_map->clone();
                cached_map_type     base_map(*my_map, period);
                rhs_type            rhs(base_map);
                jacobian_type       jacobian(base_map);
                
                std::ostringstream os;
                const nvis::vec2& x = local_mesh.get_vertex(vert_id);
                __unique_vertex_pmaps[thread_id][vert_id] = rhs(x);
                __unique_vertex_jmaps[thread_id][vert_id] = jacobian(x);
            }
        }
        std::cerr << '\n';
        std::cerr << "unique vertex sampling completed in " << sub_timer.elapsed() << " s.\n";
        std::map<unsigned int, nvis::vec2> unique_vertex_pmaps;
        std::map<unsigned int, nvis::mat2> unique_vertex_jmaps;
        for (int i = 0 ; i < nb_threads ; ++i) {
            unique_vertex_pmaps.insert(__unique_vertex_pmaps[i].begin(),
                                       __unique_vertex_pmaps[i].end());
            unique_vertex_jmaps.insert(__unique_vertex_jmaps[i].begin(),
                                       __unique_vertex_jmaps[i].end());
        }
        
        
        // scrutinize interesting edges
        
        std::vector<std::list<step_type> > steps(nb_threads);
        std::vector<std::map<edge_index_type, double> > measured_angles(nb_threads);
        std::vector<std::list<edge_index_type> > __problematic_edges(nb_threads);
        
        sub_timer.restart();
        
        // do that in parallel
        #pragma omp parallel
        {
            int thread_id = 0;
            #pragma omp for schedule(dynamic,1)
            for (int i = 0 ; i < unique_edges.size() ; ++i) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                if (!thread_id) {
                    int pct = 100 * i / unique_edges.size();
                    std::cerr << pct << "% of edges completed...     \r";
                }
                
                edge_index_type edge_id = unique_edges[i];
                
                std::ostringstream os;
                const nvis::vec2& x0 = local_mesh.get_vertex(edge_id.i0);
                const nvis::vec2& x1 = local_mesh.get_vertex(edge_id.i1);
                const nvis::vec2& v0 = unique_vertex_pmaps[edge_id.i0];
                const nvis::vec2& v1 = unique_vertex_pmaps[edge_id.i1];
                const nvis::mat2& J0 = unique_vertex_jmaps[edge_id.i0];
                const nvis::mat2& J1 = unique_vertex_jmaps[edge_id.i1];
                
                poincare_map* my_map = precise_map->clone();
                cached_map_type     base_map(*my_map, period);
                rhs_type            rhs(base_map);
                jacobian_type       jacobian(base_map);
                
                try {
                    // measured_angles[thread_id][edge_id] = pindex.safe_rotation_angle(x0, x1,
                    //                                        v0, v1, true, true, // both values are known
                    //                                        *precise_map, period, dtheta,
                    //                                        tiny_length, metric, 0, steps[thread_id],
                    //                                        true); // yes, do warn me if things go South
                    
                    // measured_angles[thread_id][edge_id] =
                    //     pindex.rotation_angle_linear_predictor(x0, x1, v0, v1,
                    //                                            true, true, rhs,
                    //                                            dtheta, tiny_length, 0,
                    //                                            steps[thread_id], false);
                    
                    // measured_angles[thread_id][edge_id] =
                    //     pindex.rotation_angle_cubic_predictor(x0, x1, v0, v1, J0, J1,
                    //                                           true, true, rhs, jacobian,
                    //                                           dtheta, tiny_length, 0,
                    //                                           steps[thread_id], true);
                    
                    measured_angles[thread_id][edge_id] =
                        pindex.rotation_angle_linear_predictor_with_jacobian(
                            x0, x1, v0, v1, true, true, rhs, jacobian,
                            dtheta, 2, 0, steps[thread_id]);
                } catch (degenerate_point_exception& e) {
                    output.p_degenerate_points[period-1].push_back(e.pos());
                    __problematic_edges[thread_id].push_back(edge_id);
                    measured_angles[thread_id][edge_id] = 0; // ignore those edges for now
                } catch (std::runtime_error& e) {
                    __problematic_edges[thread_id].push_back(edge_id);
                    measured_angles[thread_id][edge_id] = 0; // ignore those edges for now
                }
            }
        }
        std::cerr << '\n';
        std::cerr << "edge angular measure completed in " << sub_timer.elapsed() << " s.\n";
        
        std::list<edge_index_type>& problematic_edges = output.p_prob_edges[period-1];
        for (int i = 0 ; i < nb_threads ; ++i) {
            std::copy(__problematic_edges[i].begin(), __problematic_edges[i].end(),
                      std::back_inserter(problematic_edges));
        }
        std::cerr << problematic_edges.size() << " problematic edges out of "
                  << unique_edges.size() << " edges have been found\n";
                  
        std::map<edge_index_type, double> known_angles;
        for (int i = 0 ; i < nb_threads ; ++i) {
            known_angles.insert(measured_angles[i].begin(), measured_angles[i].end());
        }
        
        
        // inspect triangles
        
        sub_timer.restart();
        
        int thread_id = 0;
        std::vector<std::list<size_t> >
        sing_tris(nb_threads),
                  prob_tris(nb_threads);
        std::vector<unsigned int> candidates(output.p_cand_tris[period-1].begin(),
                                             output.p_cand_tris[period-1].end());
                                             
        // do that in parallel
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (int i = 0 ; i < candidates.size() ; ++i) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                if (!thread_id) {
                    int pct = 100 * i / candidates.size();
                    std::cerr << pct << "% of triangles completed...     \r";
                }
                
                unsigned int id = candidates[i];
                
                // size_t init_size = steps[thread_id].size();
                int real_index = pindex.safe(local_mesh, id, *precise_map, period,
                                             0.75 * M_PI, tiny_length, metric,
                                             known_angles,
                                             steps[thread_id]);
                                             
                std::ostringstream os;
                
                if (real_index != 0) {
                    sing_tris[thread_id].push_back(id);
                }
            }
        }
        std::cerr << '\n';
        
        std::cerr << "triangle inspection completed in " << sub_timer.elapsed() << " s.\n";
        
        for (int n = 0 ; n < nb_threads ; ++n) {
            std::cerr << "adding primitives generated by thread #" << n << ": ";
            std::copy(sing_tris[n].begin(), sing_tris[n].end(), std::back_inserter(output.p_sing_tris[period-1]));
            std::cerr << '\t' << sing_tris[n].size() << " singular triangles, ";
            std::copy(sing_tris[n].begin(), sing_tris[n].end(), std::back_inserter(singular_triangles));
            std::copy(prob_tris[n].begin(), prob_tris[n].end(), std::back_inserter(output.p_prob_tris[period-1]));
            std::cerr << '\t' << prob_tris[n].size() << " problem triangles, ";
            std::copy(steps[n].begin(), steps[n].end(), std::back_inserter(output.p_index_vectors[period-1]));
            std::cerr << '\t' << steps[n].size() << " poincare index vectors\n";
        }
        
        
        
        /*********************************************************************
        *
        *                       Extract fixed points
        *
        *********************************************************************/
        
        // extract all chains for this period
        
        // start Newton search at found singular triangles while sorting them
        // by increasing norm of the p-map
        std::map<double, unsigned int> weighted_seeds;
        for (int i = 0 ; i < singular_triangles.size() ; ++i) {
            size_t id = singular_triangles[i];
            point_type p[3];
            data_type d[3];
            local_mesh.get_triangle_info(p, d, id);
            double norm[3];
            double per[3];
            for (int j = 0 ; j < 3 ; ++j) {
                norm[j] = nvis::norm(vector_value(d[j], period, metric));
            }
            double minnorm = *std::min_element(norm, &norm[3]);
            weighted_seeds[minnorm] = id;
        }
        
        std::map<double, unsigned int>::const_iterator map_iter;
        std::vector<unsigned int> sorted_seeds;
        for (map_iter = weighted_seeds.begin() ; map_iter != weighted_seeds.end() ; ++map_iter) {
            sorted_seeds.push_back(map_iter->second);
        }
        
        
        typedef xavier::fixpoint                    fixpoint_type;
        typedef std::list<fixpoint_type>            chain_type;
        
        xavier::map_debug::verbose_level = 1;
        xavier::__default_metric = metric;
        
        std::vector<std::list<chain_type> > __chains(nb_threads);
        
        sub_timer.restart();
        #pragma omp parallel
        {
            int thread_id = 0;
            
            #pragma omp for schedule(dynamic,1)
            for (int i = 0 ; i < sorted_seeds.size() ; ++i) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                
                if (!thread_id) {
                    int pct = 100 * i / sorted_seeds.size();
                    std::cerr << pct << "% seeds processed...  \r";
                }
                
                size_t id = sorted_seeds[i];
                // start in the middle of the triangle
                point_type p[3];
                data_type d[3];
                local_mesh.get_triangle_info(p, d, id);
                nvis::vec2 seed = (p[0] + p[1] + p[2]) / 3.;
                
                poincare_map* my_map = precise_map->clone();
                cached_map_type     base_map(*my_map, period);
                rhs_type            rhs(base_map);
                jacobian_type       jacobian(base_map);
                
                chain_type chain;
                if (find_fixed_points(chain, rhs, jacobian, seed, 1.0e-3)) {
                    __chains[thread_id].push_back(chain);
                }
            }
        }
        std::cerr << "fixed point extraction completed in " << sub_timer.elapsed() << " s.\n";
        
        // putting things together
        for (int i = 0 ; i < nb_threads ; ++i) {
            std::copy(__chains[i].begin(), __chains[i].end(),
                      std::back_inserter(output.p_chains[period-1]));
        }
        
        
        poincare_map* my_map = precise_map->clone();
        cached_map_type     base_map(*my_map, period);
        rhs_type            rhs(base_map);
        uniquify_chains(output.p_chains[period-1], rhs, params.close_d);
        
        std::cerr << "chains found for period " << period << ":\n";
        for (std::list<chain_type>::iterator it = output.p_chains[period-1].begin() ;
                it != output.p_chains[period-1].end(); ++it) {
            std::cerr << "new chain:\n";
            for (chain_type::iterator it2 = it->begin() ; it2 != it->end() ; ++it2) {
                std::cerr << *it2 << '\n';
            }
        }
        
#if 0
        /*********************************************************************
        *
        *                       Extract separatrices
        *
        *********************************************************************/
        
        std::vector<std::vector<xavier::fixpoint> > all_p_chains;
        for (std::list<chain_type>::iterator it = output.p_chains[period-1].begin() ;
                it != output.p_chains[period-1].end(); ++it) {
            all_p_chains.push_back(std::vector<xavier::fixpoint>(it->begin(), it->end()));
        }
        
        std::vector<bool> invalid(all_p_chains.size(), false);
        
        typedef std::vector<nvis::vec2>         separatrix_type;
        
        std::list<std::pair<int, std::vector<separatrix_type> > > all_separatrices;
        
        for (int i = 0 ; i < all_p_chains.size() ; ++i) {
            if (!all_p_chains[i].size() || !all_p_chains[i][0].saddle ||
                    invalid[i]) {
                continue;
            }
            
            std::vector<separatrix_type> separatrices;
            if (process_saddle_chain(all_p_chains, separatrices, invalid, i, *my_map, metric)) {
                all_separatrices.push_back(std::make_pair(i, separatrices));
                std::cerr << "we just received " << separatrices.size()
                          << " separatrices\n";
                for (int j = 0 ; j < separatrices.size() ; ++j) {
                    std::cerr << "separatrix #" << j << " contains " << separatrices[j].size()
                              << " points\n";
                }
            }
        }
        
        std::copy(all_separatrices.begin(), all_separatrices.end(),
                  std::back_inserter(output.p_separatrices[period-1]));
#endif
                  
        std::cerr << "entire processing of period " << period << " took " << _timer.elapsed() << " s.\n";
    } // for all periods
}
































