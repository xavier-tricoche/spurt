#ifndef __EXPERIMENTAL_HPP__
#define __EXPERIMENTAL_HPP__

#include <vector>
#include <list>
#include <map>
#include <boost/random.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <poincare/metric.hpp>
#include "maps_lib/definitions.hpp"
#include "maps_lib/triangulation.hpp"

#if _OPENMP
#include <omp.h>
#endif

namespace {
template<typename T>
struct offending_point {
    offending_point(const T& data, double err)
        : _data(data), _err(err) {}
        
    const T& data() const {
        return _data;
    }
    
    double error() const {
        return _err;
    }
    
    T _data;
    double _err;
};

template<typename T>
struct Gt_point {
    bool operator()(const T& p0, const T& p1) const {
        return p0.error() > p1.error();
    }
};

template<typename Mesh>
inline nvis::vec2 center(const Mesh& mesh, unsigned int tri)
{
    typedef typename Mesh::triangle_type        triangle_type;
    
    const triangle_type& verts = mesh.get_triangle_vertices(tri);
    nvis::vec2 p[3];
    for (int n = 0 ; n < 3 ; ++n) {
        p[n] = mesh.get_vertex(verts[n]);
    }
    return (p[0] + p[1] + p[2]) / 3.;
}
}

namespace map_analysis {
namespace experimental {
template<typename Map>
void sample_on_raster(std::vector<orbit>& orbits,
                      const Map& map, const metric_type& metric,
                      unsigned int res[2], unsigned int niter)
{
    orbits.clear();
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 d = bounds.size() / nvis::vec2(res[0], res[1]);
    srand48(time(0));
    
    unsigned int ncells = res[0] * res[1];
    std::vector<bool> sampled(ncells, false);
    
    // random list of starting cells obtained through index shuffling
    std::vector<unsigned int> cell_indices(ncells);
    for (unsigned int i = 0 ; i < ncells ; ++i) {
        cell_indices[i] = i;
    }
    std::random_shuffle(cell_indices.begin(), cell_indices.end());
    
    
    unsigned int nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    std::vector<std::vector<orbit> > __orbits(nbthreads);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int i = 0 ; i < ncells ; ++i) {
        
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            if (!thread_id) {
                std::cerr << i << " / " << ncells << "     \r";
            }
            unsigned int id = cell_indices[i];
            if (sampled[id]) {
                continue;
            }
            
            nvis::vec2 xi(id % res[0] + drand48(), id % res[1] + drand48());
            nvis::vec2 x = bounds.min() + xi * d;
            
            std::vector<nvis::vec2> steps;
            const Map* the_map = map.clone();
            try {
                the_map->map(x, steps, niter);
            } catch (...) {
                continue;
            }
            std::list<nvis::vec2> all_steps(steps.begin(), steps.end());
            all_steps.push_front(x);
            if (steps.size() < niter) {
                continue;
            }
            for (std::list<nvis::vec2>::iterator it = all_steps.begin() ; it != all_steps.end() ; ++it) {
                nvis::vec2 y = metric.modulo(*it) - bounds.min();
                nvis::vec2 lgc = y / d;
                unsigned int k = floor(lgc[0]) + floor(lgc[1]) * res[0];
                sampled[k] = true;
            }
            __orbits[thread_id].push_back(orbit(all_steps));
        }
    }
    
    for (int i = 0 ; i < nbthreads ; ++i) {
        std::copy(__orbits[i].begin(), __orbits[i].end(), std::back_inserter(orbits));
    }
    
    std::cerr << orbits.size() << " orbits were integrated\n";
}

template<typename Mesh, typename Error>
inline double
triangle_check(const typename Mesh::data_type& ground_truth, const Mesh& mesh, int& tri,
               const Error& error)
{
    typedef typename Mesh::triangle_type    triangle_type;
    typedef typename Mesh::data_type        data_type;
    
    nvis::vec2 x = ground_truth.pos();
    // if we ignore the location of this point, search for it first
    if (tri < 0) {
        tri = mesh.locate_point(x);
    }
    if (tri < 0) {
        return -1;
    }
    const triangle_type& verts = mesh.get_triangle_vertices(tri);
    nvis::vec2 p[3];
    data_type d[3];
    for (int n = 0 ; n < 3 ; ++n) {
        p[n] = mesh.get_vertex(verts[n]);
        d[n] = mesh.get_data(verts[n]);
    }
    nvis::vec3 beta = spurt::barycentric(x, p);
    return error(ground_truth, beta, d);
}

template <typename Mesh, typename Sampler, typename Error>
bool refine_map(Mesh& triangles, const Sampler& sampler,
                const Error& error, const size_t max_nb_triangles)
{
    // from triangulation.hpp
    typedef Mesh                                                triangulation_type;
    typedef typename triangulation_type::data_type              data_type;
    typedef typename triangulation_type::point_type             point_type;
    typedef typename triangulation_type::index_type             index_type;
    typedef typename triangulation_type::triangle_type          triangle_type;
    typedef typename triangulation_type::bounds_type            bounds_type;
    typedef offending_point<data_type>                          added_point;
    typedef Gt_point<added_point>                               Gt_added_point;
    typedef std::set<added_point, Gt_added_point>               priority_list_type;
    typedef typename priority_list_type::iterator               priority_list_iterator;
    
    size_t nb_triangles = triangles.get_nb_triangles();
    priority_list_type priority_list;
    
    // initially designate all triangles as potential seeds in no particular order
    std::vector<index_type> seed_triangles(nb_triangles);
    for (index_type i = 0 ; i < nb_triangles ; ++i) {
        seed_triangles[i] = i;
    }
    std::random_shuffle(seed_triangles.begin(), seed_triangles.end());
    
    try {
        std::cerr << "maximum number of triangles is " << max_nb_triangles << std::endl;
        
        while (nb_triangles < max_nb_triangles && !seed_triangles.empty()) {
            std::cerr << seed_triangles.size() << " seed triangles for " << nb_triangles << " total triangles\n";
            priority_list.clear();
            
            std::map<index_type, bool> is_touched;
            for (int i = 0 ; i < seed_triangles.size() ; ++i) {
                is_touched[seed_triangles[i]] = false;
            }
            
            unsigned int skipped = 0;
            
            // check approximation quality of all the seed triangles
            
            unsigned int nbthreads = 1;
// #if _OPENMP
//          nbthreads = omp_get_max_threads();
// #endif
//
//          std::vector<std::list<added_point> > __candidates(nbthreads);
//
// #pragma omp parallel
            {
// #pragma omp for schedule(dynamic,1)
                for (int i = 0 ; i < seed_triangles.size() ; ++i) {
                
//                  int thread_id = 0;
// #if _OPENMP
//                  thread_id = omp_get_thread_num();
// #endif
//                  if (!thread_id)
                    std::cerr << "seeding in triangle #" << i << ", " << skipped << " skipped       \r";
                    
                    unsigned int n = seed_triangles[i];
                    if (is_touched[n]) {
                        ++skipped;
                        continue;
                    }
                    
                    nvis::vec2 x = center(triangles, n);
                    
                    // sample the map at that location
                    std::vector<point_type> points;
                    std::vector<data_type> data;
                    sampler(x, points, data);
                    
                    // assess approximation error at each point along orbit
                    for (int j = 0 ; j < data.size() ; ++j) {
                        int tri = (!j ? n : -1);
                        double err = triangle_check(data[j], triangles, tri, error);
                        // std::cerr << "\tpos=" << data[j].pos() << ", tri=" << tri << ", err=" << err << '\n';
                        if (tri >= 0) {
                            typename std::map<index_type, bool>::iterator it = is_touched.find(tri);
                            if (it != is_touched.end()) {
                                it->second = true;
                            }
                        }
                        if (err <= 0) {
                            continue; // good fit or unable to find the point
                        } else {
                            // __candidates[thread_id].push_back(added_point(data[j], err));
                            priority_list.insert(added_point(data[j], err));
                        }
                    }
                }
            }
            std::cerr << '\n';
            // for (int i = 0 ; i < nbthreads ; ++i) {
            //  priority_list.insert(__candidates[i].begin(), __candidates[i].end());
            // }
            
            // insert points in order while not exceeding the maximum
            // allowed number of triangles.
            // each inserted point will create 2 more triangles on average (1->3)
            int nb_allowed = (max_nb_triangles - nb_triangles) / 2;
            int nb_inserted = std::min(nb_allowed, (int)priority_list.size());
            int counter = 0;
            std::vector<point_type> points(nb_inserted);
            std::vector<data_type> data(nb_inserted);
            for (priority_list_iterator it = priority_list.begin() ;
                    counter < nb_inserted ; ++it, ++counter) {
                points[counter] = it->data().pos();
                data[counter] = it->data();
            }
            std::vector<int> modified;
            triangles.insert_points(points, data, modified);
            
            std::cerr << "inserted " << points.size() << " in triangulation\n";
            
            seed_triangles.clear();
            std::copy(modified.begin(), modified.end(),
                      std::back_inserter(seed_triangles));
            nb_triangles = triangles.get_nb_triangles();
        }
    } catch (std::runtime_error& e) {
        std::cerr << "exception caught in refine(): " << e.what() << std::endl;
        throw;
    } catch (...) {
        std::cerr << "unknown exception caught in refine()" << std::endl;
        return true;
    }
    
    if (!seed_triangles.empty()) {
        std::cerr << seed_triangles.size() << " triangles do not meet the quality criterion\n";
    }
    return seed_triangles.empty(); // all triangles now meet the set criteria
}

} // experimental

} // map_analysis




#endif






























































