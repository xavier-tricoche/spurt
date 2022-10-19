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

// map analysis
#include <math/math.hpp>
#include <poincare/metric.hpp>
#include <poincare/newton.hpp>
#include <poincare/topology.hpp>
#include <poincare/fixpoints.hpp>
#include <poincare/macros.hpp>
#include "sampling.hpp"
#include "experimental.hpp"
#include "orbits.hpp"
#include "period.hpp"
#include "regular_mesh.hpp"
#include "approximation.hpp"
#include "index.hpp"

// christoph
#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <tokamak/poincare_map.hpp>

// teem
#include <teem/nrrd.h>

// OpenMP
// #if _OPENMP
// #include <omp.h>
// #endif
#include <boost/bind.hpp>
#include <boost/thread.hpp>


// meshing
#include "definitions.hpp"
#include "triangulation.hpp"
#include "adaptive_triangulation.hpp"
#include "quality_control.hpp"
#include "IO.hpp"

using namespace map_analysis;

mesh_type* sampling_mesh;
std::vector< std::vector< xavier::fixpoint > > all_chains;

double approx_error;
std::vector<rational_type> ref_values;

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

void parallel_check_triangle_index_threadmain(volatile int* count)
{
    while (true) {
        int current = __gnu_cxx::__exchange_and_add(count, -1);
        
        if (current < 0) {
            break;
        }
        
        nptr[current]->advance(t_max);
    }
}

void check_triangle_index(const mesh_type& local_mesh, unsigned int triangle_id,
                          const poincare_index<mesh_type, poincare_map>& pindex)
{
    size_t id = triangle_id;
    mesh_type::point_type p[3];
    mesh_type::data_type d[3];
    local_mesh.get_triangle_info(p, d, id);
    
    int index = pindex.direct(d, period, metric);
    std::pair<double, double> dpair = pindex.rotation_mag(d, period, metric);
#if TRUST_DIRECT_INDEX
    if (index != 0)
#else
    if (dpair.first > 0.8 && dpair.second > 0.45)
#endif
    {
        // std::cerr << "checking triangle (" << dpair.first << ", " << dpair.second
//  << ") for possible fixed point\n";

        cand_tris[thread_id].push_back(id);
        int real_index = pindex.safe(local_mesh, id, *precise_map, period,
                                     0.5 * M_PI, tiny_length, metric);
                                     
        if (real_index != 0) {
            singular_triangles.push_back(id);
            sing_tris[thread_id].push_back(id);
            
            std::cerr << "triangle " << p[0] << ", " << p[1] << ", " << p[2]
                      << " contains a " << period << "-"
                      << (real_index > 0 ? "center" : "saddle") << '\n';
                      
            nvis::vec2 v[3];
            for (int j = 0 ; j < 3 ; ++j) {
                v[j] = pindex.evaluate_step(p[j], *precise_map, period, metric).first;
            }
            nvis::vec3 b = zero_bary_coord(v);
            if (b[0] == -1) {
                std::cerr << "invalid barycentric coordinates for zero location in singular cell!\n";
            } else {
                nvis::vec2 approx_zero = b[0] * p[0] + b[1] * p[1] + b[2] * p[2];
                std::cerr << "approximate location of " << period << "-"
                          << (real_index > 0 ? "center" : "saddle") << " is " << approx_zero << '\n';
                double zero_norm = nvis::norm(pindex.evaluate_step(approx_zero, *precise_map, period, metric).first);
                std::cerr << "corresponding norm of " << period << "-map is " << zero_norm << '\n';
                if (zero_norm > 1) {
                
                    prob_tris[thread_id].push_back(id);
                    
                    std::cerr << "nonsensical result detected: vertex vectors were: " << v[0] << ", "
                              << v[1] << ", " << v[2] << '\n';
                              
                    std::vector<nvis::vec2> pos;
                    try {
                        precise_map->map(approx_zero, pos, 10*period);
                    } catch (...) {
                        // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
                    }
                    if (pos.size() >= 9*period) {
                        xavier::push_front(approx_zero, pos);
                        double q = period_x_periodic(pos, metric).first;
                        std::cerr << "rational period at this point is " << q << '\n';
                    } else {
                        std::cerr << "something went wrong with the integration\n";
                    }
                }
            }
        }
    }
}
for (int n = 0 ; n < nb_threads ; ++n)
{
    std::copy(cand_tris[n].begin(), cand_tris[n].end(), std::back_inserter(output.p_cand_tris[period-1]));
    std::copy(sing_tris[n].begin(), sing_tris[n].end(), std::back_inserter(output.p_sing_tris[period-1]));
    std::copy(prob_tris[n].begin(), prob_tris[n].end(), std::back_inserter(output.p_prob_tris[period-1]));
}
}

}


void map_analysis::
adaptive_map_sampling(adaptive_map_sampling_output& output,
                      const adaptive_map_sampling_params& params)
{
    typedef nvis::vec2                          point_type;
    typedef point_data                          data_type;
    typedef std::pair<nvis::vec2, data_type>    data_point_type;
    
    output.p_orbits.clear();
    output.p_meshes.clear();
    output.p_cand_tris.clear();
    output.p_sing_tris.clear();
    output.p_prob_tris.clear();
    output.p_orbits.resize(params.mp);
    output.p_meshes.resize(params.mp);
    output.p_cand_tris.resize(params.mp);
    output.p_sing_tris.resize(params.mp);
    output.p_prob_tris.resize(params.mp);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(params.hdf_name), params.time_step);
    field->periodic_coordinates(false);
    
#if 0
    {
        Nrrd* nrrd = field->to_nrrd();
        std::ostringstream os;
        os << outs << "-t=" << ts << ".nrrd";
        if (nrrdSave(os.str().c_str(), nrrd, NULL)) {
            std::cerr << biffGetDone(NRRD) << std::endl;
            exit(-1);
        }
    }
#endif
    
    bool per[2] = {true, false};
    metric_type metric(field->bounds(), per);
    
    static_data::metric = metric;
    
    double h = params.eps;
    
    poincare_map pmap(field);
    pmap.precision(params.eps);
    orbit_integrator<poincare_map> intg(pmap, params.m, metric);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double tiny_length = 1.0e-6 * nvis::norm(diagonal);
    nvis::bbox2 inflated_bounds(bounds);
    inflated_bounds.min() -= 0.005 * diagonal;
    inflated_bounds.max() += 0.005 * diagonal;
    std::cerr << "inflated bounds are " << inflated_bounds << '\n';
    
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
    output.base_mesh = mesh_type(boundary, xavier::point_locator());
    mesh_type& base_mesh = output.base_mesh;
    
    unsigned int res[2];
    double ratio = bounds.size()[1] / bounds.size()[0];
    res[0] = ceil(sqrt((double)params.n / ratio));
    res[1] = ceil(res[0] * ratio);
    experimental::sample_on_raster(static_data::central_map_orbits,
                                   pmap, metric, res, params.m);
                                   
    // assign period to each orbit
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
    
    {
        std::vector<nvis::vec2> points;
        std::vector<point_data> data;
        orbits_to_points(points, data, static_data::central_map_orbits);
        base_mesh.insert_points(points, data);
        
#ifdef EXPORT_TO_VTK
        std::ostringstream os;
        os << params.out_name << "-base_mesh.vtk";
        export_VTK(base_mesh, os.str(), "N / A", false);
#endif
    }
    base_mesh.set_tolerance(1.0e-7);
    
    xavier::map_debug::verbose_level = 1;
    
    for (int period = 1 ; period <= params.mp ; ++period) {
        std::cerr << "\nprocessing period " << period << "...\n";
        
        std::cerr << "there are currently " << static_data::central_map_orbits.size()
                  << " orbits on record and " << base_mesh.get_nb_triangles() << " points in base triangulation\n";
                  
        double h_p = (period == 1 ? h : 0.25 * h / (double)(period - 1));
        intg.precision(h_p);
        std::cerr << "integrator precision is currently " << 0.25*h / (double)period << '\n';
        
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
        
#ifdef EXPORT_TO_VTK
        if (true) {
            std::cerr << "export subset of base mesh satisfying refinement criterion\n";
#if 0
            conditional_angular_variation_priority priority(period, valid_qs, dq, approx_error, metric);
#else
            angular_variation_priority priority(period, params.err, metric, params.dq);
#endif
            std::vector<unsigned int> included_triangles;
            for (int i = 0 ; i < local_mesh.get_nb_triangles() ; ++i) {
                mesh_type::point_type p[3];
                mesh_type::data_type d[3];
                local_mesh.get_triangle_info(p, d, i);
                if (priority(p, d) > 0) {
                    included_triangles.push_back(i);
                }
            }
            
            std::ostringstream os_name, os_comment;
            os_name << params.out_name << "-pre-mesh_for_p=" << period << ".vtk";
            p_step_func functor_all(period, params.dq);
            export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_all,
                               included_triangles, true);
        }
#endif
        
        std::cerr << "building triangulation for period p = " << period << '\n';
        
        std::cerr << "refinement to achieve prescribed max angular variation\n";
        angular_variation_priority priority(period, approx_error, metric, params.dq);
        bool ok = refine(local_mesh, intg, priority, params.mt);
        
        std::cerr << "after refinement mesh contains " << local_mesh.get_nb_triangles()
                  << " triangles\n";
                  
        triangle_filter filter(period, approx_error, metric, params.dq);
        std::vector<unsigned int> included_triangles;
        for (int i = 0 ; i < local_mesh.get_nb_triangles() ; ++i) {
            mesh_type::point_type p[3];
            mesh_type::data_type d[3];
            local_mesh.get_triangle_info(p, d, i);
            if (filter.valid(p, d) > 0) {
                included_triangles.push_back(i);
            }
        }
        
#ifdef EXPORT_TO_VTK
        // export triangles matching period
        std::ostringstream os_name;
        os_name << params.out_name << "-" << period << "-specific_mesh.vtk";
        p_step_func functor_all(period, params.dq);
        export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_all, included_triangles, true);
        
        os_name.clear();
        os_name.str("");
        os_name << params.out_name << "-" << period << "-specific_mesh-error.vtk";
        error_func functor_err;
        export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_err, included_triangles, true);
#endif
        
        // export subset of triangles matching period with nonzero index
        poincare_index<mesh_type, poincare_map> pindex;
        std::vector<unsigned int> singular_triangles;
        
        poincare_map* precise_map = pmap.clone();
        precise_map->precision(0.1*h_p);
        
        #pragma omp parallel
        {
            int thread_id = 0;
            std::vector<std::list<size_t> >
            cand_tris(nb_threads),
                      sing_tris(nb_threads),
                      prob_tris(nb_threads);
                      
            #pragma omp for schedule(dynamic,1)
            for (int i = 0 ; i < included_triangles.size() ; ++i) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                size_t id = included_triangles[i];
                mesh_type::point_type p[3];
                mesh_type::data_type d[3];
                local_mesh.get_triangle_info(p, d, id);
                
                int index = pindex.direct(d, period, metric);
                std::pair<double, double> dpair = pindex.rotation_mag(d, period, metric);
#if TRUST_DIRECT_INDEX
                if (index != 0)
#else
                if (dpair.first > 0.8 && dpair.second > 0.45)
#endif
                {
                    // std::cerr << "checking triangle (" << dpair.first << ", " << dpair.second
//  << ") for possible fixed point\n";

                    cand_tris[thread_id].push_back(id);
                    int real_index = pindex.safe(local_mesh, id, *precise_map, period,
                                                 0.5 * M_PI, tiny_length, metric);
                                                 
                    if (real_index != 0) {
                        singular_triangles.push_back(id);
                        sing_tris[thread_id].push_back(id);
                        
                        std::cerr << "triangle " << p[0] << ", " << p[1] << ", " << p[2]
                                  << " contains a " << period << "-"
                                  << (real_index > 0 ? "center" : "saddle") << '\n';
                                  
                        nvis::vec2 v[3];
                        for (int j = 0 ; j < 3 ; ++j) {
                            v[j] = pindex.evaluate_step(p[j], *precise_map, period, metric).first;
                        }
                        nvis::vec3 b = zero_bary_coord(v);
                        if (b[0] == -1) {
                            std::cerr << "invalid barycentric coordinates for zero location in singular cell!\n";
                        } else {
                            nvis::vec2 approx_zero = b[0] * p[0] + b[1] * p[1] + b[2] * p[2];
                            std::cerr << "approximate location of " << period << "-"
                                      << (real_index > 0 ? "center" : "saddle") << " is " << approx_zero << '\n';
                            double zero_norm = nvis::norm(pindex.evaluate_step(approx_zero, *precise_map, period, metric).first);
                            std::cerr << "corresponding norm of " << period << "-map is " << zero_norm << '\n';
                            if (zero_norm > 1) {
                            
                                prob_tris[thread_id].push_back(id);
                                
                                std::cerr << "nonsensical result detected: vertex vectors were: " << v[0] << ", "
                                          << v[1] << ", " << v[2] << '\n';
                                          
                                std::vector<nvis::vec2> pos;
                                try {
                                    precise_map->map(approx_zero, pos, 10*period);
                                } catch (...) {
                                    // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
                                }
                                if (pos.size() >= 9*period) {
                                    xavier::push_front(approx_zero, pos);
                                    double q = period_x_periodic(pos, metric).first;
                                    std::cerr << "rational period at this point is " << q << '\n';
                                } else {
                                    std::cerr << "something went wrong with the integration\n";
                                }
                            }
                        }
                    }
                }
            }
            for (int n = 0 ; n < nb_threads ; ++n) {
                std::copy(cand_tris[n].begin(), cand_tris[n].end(), std::back_inserter(output.p_cand_tris[period-1]));
                std::copy(sing_tris[n].begin(), sing_tris[n].end(), std::back_inserter(output.p_sing_tris[period-1]));
                std::copy(prob_tris[n].begin(), prob_tris[n].end(), std::back_inserter(output.p_prob_tris[period-1]));
            }
        }
        
#ifdef EXPORT_TO_VTK
        os_name.clear();
        os_name.str("");
        os_name << params.out_name << "-" << period << "-specific_singular.vtk";
        export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_all, singular_triangles, true);
#endif
        
#ifdef EXPORT_TO_VTK
        if (false) {
            mesh_type tmp_mesh(boundary, xavier::point_locator());
            tmp_mesh.set_tolerance(1.0e-8);
            std::vector<point_type> points;
            std::vector<data_type> data;
            orbits_to_points(points, data, map_analysis::static_data::central_map_orbits);
            tmp_mesh.insert_points(points, data);
            
            std::ostringstream os_name, os_comment;
            os_name << params.out_name << "-complete_mesh_after_" << period << "-search.vtk";
            p_step_func functor_all(period, 100);
            export_VTK(tmp_mesh, os_name.str(), "N/A", functor_all, true);
        }
#endif
    } // for all periods
}































































