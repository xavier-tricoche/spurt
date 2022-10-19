#ifndef __ADAPTIVE_TRIANGULATION_HPP__
#define __ADAPTIVE_TRIANGULATION_HPP__

#include "maps_lib/triangulation.hpp"
#include <vector>
#include <set>
#include "maps_lib/priority_queue.hpp"

namespace map_analysis {
struct triangle_info {
    triangle_info(int idx, double p, size_t v = 0)
        : index(idx), priority(p), version(v) {}
        
    int index;
    double priority;
    size_t version;
};

struct Gt_triangle_priority {
    bool operator()(const triangle_info& i0, const triangle_info& i1) {
        // sort priorities in descending order for selection
        return (i0.priority > 0 && i0.priority > i1.priority);
    }
};

template<typename Mesh, typename TriEval>
inline double priority(size_t i, const Mesh& mesh, const TriEval& eval)
{
    typedef typename Mesh::data_type        data_type;
    typedef typename Mesh::point_type       point_type;
    
    point_type p[3];
    data_type v[3];
    mesh.get_triangle_info(p, v, i);
    return eval(p, v);
}

struct still_valid {
    still_valid(const std::vector<size_t>& current) : __current(current) {}
    
    bool operator()(const triangle_info& info) const {
        return (info.priority > 0 &&
                info.version == __current[info.index]);
    }
    
    const std::vector<size_t> __current;
};

template<typename T, typename PtInsert, typename TriEval>
inline bool refine(T& triangles, const PtInsert& inserter, const TriEval& evaluator,
                   const size_t max_nb_triangles)
{
    // from triangulation.hpp
    typedef T                                                   triangulation_type;
    typedef typename triangulation_type::data_type              data_type;
    typedef typename triangulation_type::point_type             point_type;
    typedef typename triangulation_type::index_type             index_type;
    typedef typename triangulation_type::triangle_type          triangle_type;
    typedef typename triangulation_type::bounds_type            bounds_type;
    typedef xavier::priority_queue<triangle_info, Gt_triangle_priority>     priority_queue_type;
    
    size_t nb_triangles = triangles.get_nb_triangles();
    std::vector<size_t> current_version(nb_triangles, 0);
    priority_queue_type priority_list;
    
    try {
        // evaluate each triangle and initialize priority list
        for (int i = 0 ; i < nb_triangles; ++i) {
            double p = priority(i, triangles, evaluator);
            if (p > 0) {
                // only consider triangles with positive (valid) priority value
                priority_list.add(triangle_info(i, p));
            }
        }
        
        std::cerr << priority_list.size() << " triangles out of "
                  << nb_triangles << " were added to priority list\n";
        std::cerr << "maximum number of triangles is " << max_nb_triangles << std::endl;
        
        while (nb_triangles < max_nb_triangles && !priority_list.empty()) {
        
            std::cerr << "there are " << nb_triangles << " (" << max_nb_triangles << ") so far and "
                      << "highest priority triangle is " << priority_list.first().index
                      << " with priority " << priority_list.first().priority<< '\n';
            std::cerr << "its version is " << priority_list.first().version
                      << " while the current version is " << current_version[priority_list.first().index] << '\n';
                      
            // select the element with highest priority
            const triangle_info& info = priority_list.remove();
            if (info.version < current_version[info.index]) {
                // we have already replaced that entry by a more recent one. remove it
                continue;
            } else if (info.priority <= 0) {
                // invalid priority value. remove it
                continue;
            } else {
                point_type p[3];
                const triangle_type& triangle = triangles.get_triangle_vertices(info.index);
                for (int n = 0 ; n < 3 ; ++n) {
                    p[n] = triangles.get_vertex(triangle[n]);
                }
                point_type seed = 1. / 3.*(p[0] + p[1] + p[2]);
                std::vector<point_type> points;
                std::vector<data_type> data;
                
                inserter(seed, points, data);
                
                if (!points.size()) {
                    // point insertion failed.
                    continue;
                }
                std::vector<int> modified;
                triangles.insert_points(points, data, modified);
                // extend version tracker and initialize new entries at 0
                current_version.resize(triangles.get_nb_triangles(), 0);
                // update record for modified triangles
                for (int n = 0 ; n < modified.size() ; ++n) {
                    int id = modified[n];
                    if (id < nb_triangles) {
                        ++current_version[id];
                    }
                    double p = priority(id, triangles, evaluator);
                    if (p <= 0) {
                        continue;
                    }
                    priority_list.add(triangle_info(id, p, current_version[id]));
                }
            }
            
            nb_triangles = triangles.get_nb_triangles();
        }
        
        // remove from list all remaining outdated triangles
        still_valid valid_checker(current_version);
        priority_list.remove(valid_checker);
    } catch (std::runtime_error& e) {
        std::cerr << "exception caught in refine(): " << e.what() << std::endl;
        throw;
    } catch (...) {
        std::cerr << "unknown exception caught in refine()" << std::endl;
        return true;
    }
    
    if (!priority_list.empty()) {
        std::cerr << priority_list.size() << " triangles do not meet the quality criterion\n";
        std::cerr << "worst offender has priority value: " << priority_list.first().priority << "\n";
    }
    return priority_list.empty(); // all triangles now meet the set criteria
}
}

#endif
