#ifndef __TRIANGULATION_HPP__
#define __TRIANGULATION_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>
#include <vector>
#include <list>
#include <set>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>
// #include <boost/random.hpp>

#include <kdtree++/kdtree.hpp>
#include <data/kdtree.hpp>

// for debugging only
#include "maps_lib/IO.hpp"

// #define DISPLAY_TIMINGS

namespace spurt {
inline nvis::vec3 barycentric(const nvis::vec2& x, const nvis::vec2 p[3])
{
    // x-p0 = b1*(p1-p0) + b2*(p2-p0)
    nvis::vec2 y = x - p[0];
    nvis::vec2 q[2] = { p[1] - p[0], p[2] - p[0] };
    double det = q[0][0] * q[1][1] - q[0][1] * q[1][0];
    nvis::vec3 b;
    b[1] = (y[0] * q[1][1] - y[1] * q[1][0]) / det;
    b[2] = (q[0][0] * y[1] - q[0][1] * y[0]) / det;
    b[0] = 1. - b[1] - b[2];
    return b;
}
}

namespace spurt {
struct triangulation_exception: public std::runtime_error {

    triangulation_exception() : std::runtime_error(""), offending_position(0, 0) {}
    triangulation_exception(const std::string& msg) : std::runtime_error(msg) {}
    virtual ~triangulation_exception() throw() {};

    std::vector<int> visited_triangles;
    nvis::vec2 offending_position;
};

}



// struct default_point_locator {
//  typedef int             index_type;
//  typedef nvis::vec2      pos_type;
//  index_type find_close_point(const pos_type& x, index_type id) const {
//      return 0;
//  }
// };
// }

namespace spurt {

template < typename T1, typename T2>
class triangulation {
public:
    typedef int                     index_type;
    typedef nvis::vec2              point_type;
    typedef T1                      data_type;
    typedef nvis::ivec3             triangle_type;
    typedef nvis::bbox2             bounds_type;
    typedef T2                      point_locator_type;
    typedef typename T2::data_point_type data_point_type;

    enum export_mode {
        NONE, EACH, ONE_PER_BATCH, SOME
    };

    triangulation() : __bounds_set(false) {}

    triangulation(const std::pair<point_type, data_type> boundary[4],
                  const point_locator_type& locator)
        : __nb_duplicate_pts(0), __nb_degeneracies(0), __bounds_set(false),
          __tolerance(0.00001), __step_counter(0), __locator(locator),
          __export_mode(NONE), __export_file_name("debug"), __export_frequency(0) {
        set_bounding_box(boundary);
    }

    void set_bounding_box(const std::pair<point_type, data_type> boundary[4]);
    const bounds_type& get_bounds() const {
        return __bounds;
    }

    void set_tolerance(double tol) {
        __tolerance = fabs(tol);
    }

    size_t get_nb_steps() const {
        return __step_counter;
    }

    const triangle_type& get_triangle_vertices(const index_type tri_id) const {
        assert(tri_id < __tri_to_verts.size());
        return __tri_to_verts[tri_id];
    }

    void get_triangle_info(point_type p[3], data_type v[3],
                           const index_type& tri_id) const {
        for (int i = 0 ; i < 3 ; ++i) {
            index_type id = __tri_to_verts[tri_id][i];
            p[i] = __verts[id];
            v[i] = __data[id];
        }
    }

    const point_type& get_vertex(const index_type& id) const {
        return __verts[id];
    }

    const data_type& get_data(const index_type& id) const {
        return __data[id];
    }

    void get_vertex_neighbors(std::vector<index_type>& ids,
                              const index_type id) const {
        assert(id < __verts.size());
        ids.clear();
        std::copy(__vert_to_tris[id].begin(), __vert_to_tris[id].end(),
                  std::back_inserter(ids));
    }

    index_type get_edge_neighbor(const index_type triangle_id,
                                 const index_type edge_id) const;

    index_type get_edge_neighbor(const index_type triangle_id,
                                 const index_type pos_id0,
                                 const index_type pos_id1) const;

    index_type get_nb_vertices() const {
        return __verts.size();
    }

    index_type get_nb_triangles() const {
        return __tri_to_verts.size();
    }

    void build_connectivity() {
        __vert_to_tris.resize(__verts.size());
        int nb_tris = __tri_to_verts.size();
        for (index_type i = 0 ; i < nb_tris ; ++i) {
            build_connectivity(i);
        }
    }

    void build_connectivity(const index_type& tri_id) {
        const triangle_type& ids = get_triangle_vertices(tri_id);
        for (int j = 0 ; j < 3 ; ++j) {
            __vert_to_tris[ids[j]].push_back(tri_id);
        }
    }

    void insert_points(const std::vector<point_type>& verts,
                       const std::vector<data_type>& data);

    void insert_points(const std::vector<point_type>& verts,
                       const std::vector<data_type>& data,
                       std::vector<index_type>& new_tris);

    void insert_points(const std::vector<std::pair<point_type, data_type> >& data_points);

    void insert_points(const std::vector<std::pair<point_type, data_type> >& data_points,
                       std::vector<index_type>& new_tris);

    void set_export_basename(const std::string& name) const {
        __export_file_name = name;
    }

    void set_export_mode(export_mode mode, size_t f) const;

    index_type locate_point(const nvis::vec2& x) const;

private:
    std::vector< point_type >                   __verts;
    std::vector< data_type >                    __data;
    std::vector< triangle_type >                __tri_to_verts;
    std::vector< std::list< index_type > >      __vert_to_tris;
    mutable int                                 __nb_duplicate_pts;
    mutable int                                 __nb_degeneracies;
    bounds_type                                 __bounds;
    bool                                        __bounds_set;
    double                                      __tolerance;
    std::set<index_type>                        __modified;
    // diagnostic info
    mutable size_t                              __step_counter;
    mutable export_mode                         __export_mode;
    mutable std::string                         __export_file_name;
    mutable size_t                              __export_frequency;
    // optimization
    point_locator_type                          __locator;

    void remove_tri_reference(index_type p_id, index_type tri_id) {
        std::list<index_type>::iterator it;
        it = std::find(__vert_to_tris[p_id].begin(), __vert_to_tris[p_id].end(), tri_id);
        assert(it != __vert_to_tris[p_id].end());
        __vert_to_tris[p_id].erase(it);
    }

    void set_triangle_ids(const triangle_type& tri_ids, index_type id) {
        assert(id < __tri_to_verts.size());
        __tri_to_verts[id] = tri_ids;
    }

    double outer(const point_type& x0, const point_type& x1) const {
        return x0[0]*x1[1] - x0[1]*x1[0];
    }

    double circumcircle(point_type& center, const point_type& x1,
                        const point_type& x2, const point_type& x3) const;

    bool in_circle(const point_type& x, const point_type x1,
                   const point_type& x2, const point_type& x3) const;

    void check_edge(index_type p_id,
                    index_type p1, index_type p2, index_type tri);

    void triangulate_new_points(index_type first_id);

    void shuffle(std::vector<index_type>& shuffled, size_t n) const;

    void barycentric(double b[3], const nvis::vec2 p[3], const nvis::vec2& x) const;

    void statistics(double& avg_area, double& min_area, double& max_area, double& var_area,
                    double& avg_edge_length, double& min_edge_length, double& max_edge_length,
                    double& var_edge_length) const;

    // returned values:
    // function: cell id if found, otherwise 0.
    // other: if negative, point was found inside the cell, otherwise on an edge
    int find_triangle(const point_type& x, index_type& tri,
                      double tol, index_type& other,
                      std::set<index_type>& been_there,
                      std::set<size_t>& visited) const;
};

template<typename T1, typename T2>
inline void triangulation<T1, T2>::
set_bounding_box(const std::pair<point_type, data_type> boundary[4])
{
    assert(!__bounds_set && !__verts.size());

    __verts.resize(4);
    __data.resize(4);
    __bounds.reset();

    for (int i = 0 ; i < 4 ; ++i) {
        __bounds.add(boundary[i].first);
        __verts[i] = boundary[i].first;
        __data[i] = boundary[i].second;
        __locator.insert(data_point_type(__verts[i], i));
    }
    std::cerr << "after initialization the bounding box of the triangulation is " << __bounds << std::endl;
    point_type diag = __bounds.size();
    assert(fabs(diag[0]*diag[1]) != 0.);

    // create initial triangulation
    __tri_to_verts.resize(2);
    __tri_to_verts[0] = triangle_type(0, 1, 2);
    __tri_to_verts[1] = triangle_type(0, 2, 3);

    build_connectivity();

    __bounds_set = true;
    __modified.insert(0);
    __modified.insert(1);

    export_VTK(*this, "very_first.vtk", "none", true, false);
}

template<typename T1, typename T2>
inline typename triangulation<T1, T2>::index_type
triangulation<T1, T2>::locate_point(const nvis::vec2& x) const
{
    index_type close_id = __locator.find_close_point(x, 0, false);
    index_type tri = 0;
    if (!__vert_to_tris[close_id].empty()) {
        tri = __vert_to_tris[close_id].front();
    }
    index_type other = -1;
    std::set<index_type> been_there;

    // for pathological cases...
    std::set<size_t> visited;
    return find_triangle(x, tri, __tolerance, other, been_there, visited);
}

template<typename T1, typename T2>
inline typename triangulation<T1, T2>::index_type
triangulation<T1, T2>::get_edge_neighbor(const index_type triangle_id,
        const index_type edge_id) const
{
    assert(triangle_id < get_nb_triangles());
    triangle_type ids = get_triangle_vertices(triangle_id);
    // intersect neighbors of both edge vertices
    index_type i1 = (edge_id + 1) % 3, i2 = (edge_id + 2) % 3;
    std::set<index_type> s0(__vert_to_tris[i1].begin(), __vert_to_tris[i1].end());
    std::set<index_type> s1(__vert_to_tris[i2].begin(), __vert_to_tris[i2].end());
    std::vector<index_type> two_cells;
    std::set_intersection(s0.begin(), s0.end(), s1.begin(), s1.end(),
                          std::back_inserter(two_cells));
    // determine which one of the two is the requested neighbor
    assert(two_cells.size());
    if (two_cells.size() == 1 && two_cells[0] == triangle_id) {
        return -1;
    } else if (two_cells[0] == triangle_id) {
        return two_cells[1];
    } else {
        return two_cells[0];
    }
}

template<typename T1, typename T2>
inline typename triangulation<T1, T2>::index_type
triangulation<T1, T2>::get_edge_neighbor(const index_type triangle_id,
        const index_type pos_id0,
        const index_type pos_id1) const
{
    assert(triangle_id < get_nb_triangles());
    triangle_type ids = get_triangle_vertices(triangle_id);
    // verify that the position ids do indeed correspond to vertices of the triangle
    std::set<index_type> tri(&ids[0], &ids[3]);
    assert(tri.find(pos_id0) != tri.end());
    assert(tri.find(pos_id1) != tri.end());
    // intersect neighbors of both edge vertices
    std::set<index_type> s0(__vert_to_tris[pos_id0].begin(), __vert_to_tris[pos_id0].end());
    std::set<index_type> s1(__vert_to_tris[pos_id1].begin(), __vert_to_tris[pos_id1].end());
    std::vector<index_type> two_cells;
    std::set_intersection(s0.begin(), s0.end(), s1.begin(), s1.end(),
                          std::back_inserter(two_cells));
    // determine which one of the two is the requested neighbor
    assert(two_cells.size());
    if (two_cells.size() == 1 && two_cells[0] == triangle_id) {
        return -1;
    } else if (two_cells[0] == triangle_id) {
        return two_cells[1];
    } else {
        return two_cells[0];
    }
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::insert_points(const std::vector<point_type>& verts,
        const std::vector<data_type>& data,
        std::vector<index_type>& new_tris)
{
#ifdef DISPLAY_TIMINGS
    nvis::timer _timer;
    std::cerr << "inserting " << verts.size() << " points in triangulation... " << std::flush;
#endif

    if (!verts.size()) {
        // std::cerr << "triangulation<>::insert_points(): empty point list\n";
        return;
    }

    __modified.clear();
    new_tris.clear();

    insert_points(verts, data);

    std::copy(__modified.begin(), __modified.end(), std::back_inserter(new_tris));

#ifdef DISPLAY_TIMINGS
    double dt = _timer.elapsed();
    std::cerr << "done after " << dt << " s. (" << dt / (double)verts.size() << " s. / point)\n";
#endif
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::insert_points(const std::vector<point_type>& verts,
        const std::vector<data_type>& data)
{
#ifdef DISPLAY_TIMINGS
    nvis::timer _timer;
    std::cerr << "inserting " << verts.size() << " points in triangulation... " << std::flush;
#endif

    if (!__bounds_set) {
        throw std::runtime_error("bounding box of triangulation not set");
    }
    assert(data.size() == verts.size());

    index_type first_id = __verts.size();
    std::copy(verts.begin(), verts.end(), std::back_inserter(__verts));
    std::copy(data.begin(), data.end(), std::back_inserter(__data));
    __vert_to_tris.resize(__verts.size());

    triangulate_new_points(first_id);

#ifdef DISPLAY_TIMINGS
    double dt = _timer.elapsed();
    std::cerr << "done after " << dt << " s. (" << dt / (double)verts.size() << " s. / point)\n";
#endif
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::insert_points(
    const std::vector<std::pair<point_type, data_type> >& data_points,
    std::vector<index_type>& new_tris)
{
#ifdef DISPLAY_TIMINGS
    nvis::timer _timer;
    std::cerr << "inserting " << data_points.size() << " points in triangulation... " << std::flush;
#endif

    __modified.clear();
    new_tris.clear();

    insert_points(data_points);

    std::copy(__modified.begin(), __modified.end(), std::back_inserter(new_tris));

#ifdef DISPLAY_TIMINGS
    double dt = _timer.elapsed();
    std::cerr << "done after " << dt << " s. (" << dt / (double)data_points.size() << " s. / point)\n";
#endif
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::insert_points(
    const std::vector<std::pair<point_type, data_type> >& data_points)
{
#ifdef DISPLAY_TIMINGS
    nvis::timer _timer;
    std::cerr << "inserting " << data_points.size() << " points in triangulation... " << std::flush;
#endif

    if (!__bounds_set) {
        throw std::runtime_error("bounding box of triangulation not set");
    }

    index_type first_id = __verts.size();
    __verts.resize(first_id + data_points.size());
    __data.resize(first_id + data_points.size());

    for (int i = 0 ; i < data_points.size() ; ++i) {
        __verts[first_id+i] = data_points[i].first;
        __data[first_id+i]  = data_points[i].second;
    }

    __vert_to_tris.resize(__verts.size());
    triangulate_new_points(first_id);

#ifdef DISPLAY_TIMINGS
    double dt = _timer.elapsed();
    std::cerr << "done after " << dt << " s. (" << dt / (double)data_points.size() << " s. / point)\n";
#endif
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::triangulate_new_points(index_type first_id)
{
    std::vector<index_type> shuffled;
    index_type nb_pts = __verts.size() - first_id;

#ifdef DISPLAY_TIMINGS
    nvis::timer _timer;
    std::cerr << "initializing triangulation with " << nb_pts << " points... " << std::flush;
#endif

    shuffle(shuffled, nb_pts);
    std::set<size_t> visited;
    for (index_type n = 0 ; n < nb_pts ; ++n) {

        // std::cerr << "inserting " << n << "th point from " << nb_pts << "    \r";

        index_type p_id = first_id + shuffled[n];
        const point_type& x = get_vertex(p_id);

        // std::cerr << "inserting " << n << "th position with index " << p_id << " and position "
        //           << x << std::endl;

        index_type close_id;
        try {
            close_id = __locator.find_close_point(x);
        } catch(...) {
            close_id = 0;
        }
        __locator.insert(data_point_type(x, p_id));
        index_type tri = 0;
        if (__vert_to_tris[close_id].empty()) {
            // std::cerr << "found neighbor (" << close_id << ") of point #"
            //           << p_id << " is disconnected!\n";
        } else {
            tri = __vert_to_tris[close_id].front();
        }
        index_type other = -1;

        // std::cerr << "starting search in triangle #" << tri
        // << " which contains following positions: ";
        // for (int i = 0 ; i < 3 ; ++i) {
        //  std::cerr << "(#" << __tri_to_verts[tri][i] << ", " << __verts[__tri_to_verts[tri][i]] << "), ";
        // }
        // std::cerr << '\n';

        std::set<index_type> been_there;

        if (__export_mode == EACH ||
                (__export_mode == SOME && __export_frequency > 0 && !(n % __export_frequency))) {
            std::ostringstream os;
            os << __export_file_name << "-" << n << ".vtk";
            export_VTK(*this, os.str(), "none", true, false);
        }

        visited.clear();

        if ((tri = find_triangle(x, tri, __tolerance, other, been_there, visited)) >= 0) {
            if (other < 0) { //in triangle

                //delete this triangle; create three new triangles
                //first triangle is replaced with one of the new ones
                triangle_type tri_ids = __tri_to_verts[tri];
                remove_tri_reference(tri_ids[2], tri);
                __tri_to_verts[tri] = triangle_type(p_id, tri_ids[0], tri_ids[1]);
                __vert_to_tris[p_id].push_back(tri);

                // {
                //  double b[3];
                //  nvis::vec2 p[3] = {__verts[tri_ids[0]], __verts[tri_ids[1]], __verts[tri_ids[2]]};
                //  barycentric(b, p, x);
                //  std::cerr << std::setprecision(12) << x << " is inside triangle "
                //            << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
                //  std::cerr << std::setprecision(12) << "barycentric coordinates: ("
                //            << b[0] << ", " << b[1] << ", " << b[2] << ")\n";
                // }

                //create two new triangles
                index_type tri1 = __tri_to_verts.size();
                index_type tri2 = tri1 + 1;
                __tri_to_verts.push_back(triangle_type(p_id, tri_ids[1], tri_ids[2]));
                __tri_to_verts.push_back(triangle_type(p_id, tri_ids[2], tri_ids[0]));
                build_connectivity(tri1);
                build_connectivity(tri2);

                // keep track of modified / created triangles
                __modified.insert(tri);
                __modified.insert(tri1);
                __modified.insert(tri2);

                // Check edge neighbors for Delaunay criterion. If not satisfied, flip
                // edge diagonal. (This is done recursively.)
                check_edge(p_id, tri_ids[0], tri_ids[1], tri);
                check_edge(p_id, tri_ids[1], tri_ids[2], tri1);
                check_edge(p_id, tri_ids[2], tri_ids[0], tri2);
            } else { // on triangle edge
                // identify edge

                const triangle_type& tri1_ids = __tri_to_verts[tri];
                const triangle_type& tri2_ids = __tri_to_verts[other];
                std::set<index_type> stri1(&tri1_ids[0], &tri1_ids[3]);
                std::set<index_type> stri2(&tri2_ids[0], &tri2_ids[3]);
                std::vector<index_type> edge_ids;
                std::set_intersection(stri1.begin(), stri1.end(),
                                      stri2.begin(), stri2.end(),
                                      std::back_inserter(edge_ids));
                assert(edge_ids.size() == 2);
                index_type id_1, id_2;
                for (int i = 0 ; i < 3 ; ++i) {
                    index_type id = tri1_ids[i];
                    if (id != edge_ids[0] && id != edge_ids[1]) {
                        id_1 = id;
                    }
                    id = tri2_ids[i];
                    if (id != edge_ids[0] && id != edge_ids[1]) {
                        id_2 = id;
                    }
                }

                // std::cerr << x << "is on edge: "
                //           << __verts[id_1] << " - " << __verts[id_2] << std::endl;

                //replace two triangles
                remove_tri_reference(edge_ids[1], tri);
                remove_tri_reference(edge_ids[1], other);
                __tri_to_verts[tri]  = triangle_type(p_id, id_1, edge_ids[0]);
                __tri_to_verts[other] = triangle_type(p_id, edge_ids[0], id_2);
                __vert_to_tris[p_id].push_back(tri);
                __vert_to_tris[p_id].push_back(other);

                //create two new triangles
                index_type tri1 = __tri_to_verts.size();
                index_type tri2 = tri1 + 1;
                __tri_to_verts.push_back(triangle_type(p_id, edge_ids[1], id_1));
                __tri_to_verts.push_back(triangle_type(p_id, id_2, edge_ids[1]));
                build_connectivity(tri1);
                build_connectivity(tri2);

                // keep track of modified / created triangles
                __modified.insert(tri);
                __modified.insert(other);
                __modified.insert(tri1);
                __modified.insert(tri2);

                // Check edge neighbors for Delaunay criterion.
                check_edge(p_id, id_1, edge_ids[0], tri);
                check_edge(p_id, id_2, edge_ids[0], other);
                check_edge(p_id, id_1, edge_ids[1], tri1);
                check_edge(p_id, id_2, edge_ids[1], tri2);
            }
        } //if triangle found
        else {
            if (tri == -4) {
                // std::cerr << "position " << n << "/" << nb_pts << " is known already\n";
            } else if (tri == -5) {
                // std::cerr << "position " << n << "/" << nb_pts << " led to degenerate case\n";
            } else if (tri == -6 || tri == -10) {
                std::ostringstream os;
                os << "position " << x << " led to circular search through " << been_there.size() << " triangles";
                std::cerr << os.str() << "\n";
                continue; // let's not get dramatic about it

                spurt::triangulation_exception e(os.str());
                std::copy(been_there.begin(), been_there.end(), std::back_inserter(e.visited_triangles));
                e.offending_position = x;
                throw e;

                std::cerr << "position " << n << "/" << nb_pts << " led to endless search\n";
            }
            tri = 0; //no triangle found
        }
    }//for all points

    if (__export_mode == ONE_PER_BATCH) {
        std::ostringstream os;
        os << __export_file_name << ".vtk";
        export_VTK(*this, os.str(), "none", true, false);
    }

#ifdef DISPLAY_TIMINGS
    double dt = _timer.elapsed();
    std::cerr << "done after " << dt << " s. (" << dt / (double)nb_pts << " s. / point)\n";
#endif
}

#define DEL2D_TOLERANCE 1.0e-9

// from VTK::vtkDelaunay2D::FindTriangles()
template<typename T1, typename T2>
inline int triangulation<T1, T2>::find_triangle(const point_type& x, index_type& tri,
        double tol, index_type& other, std::set<index_type>& been_there,
        std::set<size_t>& visited) const
{
    assert(tri < get_nb_triangles());

    if (been_there.find(tri) != been_there.end()) {
        // std::cerr << "turning in circle to find triangle for " << x << ". Giving up.\n";
        // throw std::runtime_error("turning in circles");
        return -10;
    }
    been_there.insert(tri);

    // std::cerr << "visiting triangle #" << tri << std::endl;

    int i, j, ir, ic, inside, i2, i3;
    index_type npts, newNei, nei[3];
    point_type p[3], n, vp, vx;
    double dp, minProj;

    while (true) {

        // check that triangle index is valid.
        // otherwise we have left the triangulated domain
        if (tri < 0) {
            // std::cerr << "invalid triangle index\n";
            return -2;
        } else if (visited.find(tri) != visited.end()) {
            // std::cerr << "turning in circle to insert " << x << ". giving up\n";
            return -6;
        }
        visited.insert(tri);

        // get local triangle info
        const triangle_type& tri_ids = get_triangle_vertices(tri);
        for (int k = 0 ; k < 3 ; ++k) {
            p[k] = get_vertex(tri_ids[k]);
        }

        // double b[3];
        // barycentric(b, p, x);
        // std::cerr << "triangle #" << tri << ", x=" << x
        //           << ", (" << b[0] << ", " << b[1] << ", " << b[2] << ")\n";

        // Randomization (of find edge neighbors) avoids walking in
        // circles in certain weird cases
        srand(tri);
        ir = rand() % 3;
        // evaluate in/out of each edge

        nvis::vec3 __dp;
        for (inside = 1, minProj = 0.0, ic = 0; ic < 3; ic++) {
            i  = (ir + ic) % 3;
            i2 = (i + 1) % 3;
            i3 = (i + 2) % 3;

            // create a 2D edge normal to define a "half-space"; evaluate points (i.e.,
            // candidate point and other triangle vertex not on this edge).
            n[0] = -(p[i2][1] - p[i][1]);
            n[1] = p[i2][0] - p[i][0];
            n /= nvis::norm(n);

            // compute local vectors
            vp = p[i3] - p[i];
            vp /= nvis::norm(vp);
            vx = x - p[i];

            //check for duplicate point
            double vxnorm = nvis::norm(vx);
            vx /= vxnorm;
            if (vxnorm <= tol) {
                // std::cerr << "duplicate point!\n";
                // std::cerr << vxnorm << '\n';
                __nb_duplicate_pts++;
                return -4;
            }

            // see if two points are in opposite half spaces
            dp = nvis::inner(n, vx) * (nvis::inner(n, vp) < 0 ? -1.0 : 1.0);
            // __dp[ic] = dp;
            if (dp < DEL2D_TOLERANCE) {
                if (dp < minProj) { //track edge most orthogonal to point direction
                    inside = 0;
                    nei[1] = tri_ids[i];
                    nei[2] = tri_ids[i2];
                    minProj = dp;
                }
            }//outside this edge
        }//for each edge

        if (inside) { // all edges have tested positive
            other = -1;
            // std::cerr << "triangle " << tri << " contains " << x << ", dps = " << __dp << std::endl;
            return tri;
        } else if (!inside && (fabs(minProj) < DEL2D_TOLERANCE)) { // on edge
            other = get_edge_neighbor(tri, nei[1], nei[2]);
            return tri;
        } else { //walk towards point
            tri = get_edge_neighbor(tri, nei[1], nei[2]);
            if (tri == other) {
                // std::cerr << "degeneracy!\n";
                __nb_degeneracies++;
                return -5;
            } else {
                other = tri;
                // return find_triangle(x, tri, tol, other, been_there);
            }
        }
    }
}

template<typename T1, typename T2>
inline double triangulation<T1, T2>::circumcircle(point_type& center, const point_type& x1,
        const point_type& x2, const point_type& x3) const
{
    // cf. http://en.wikipedia.org/wiki/Circumscribed_circle
    // use local coordinates
    point_type y2 = x2 - x1;
    point_type y3 = x3 - x1;
    double D = 2. * outer(y2, y3);
    if (fabs(D) < DEL2D_TOLERANCE) {
        return -1;    // degenerate case: triangle is flat
    }
    point_type y;
    y[0] = (y3[1] * (y2[0] * y2[0] + y2[1] * y2[1]) - y2[1] * (y3[0] * y3[0] + y3[1] * y3[1])) / D;
    y[1] = (y2[0] * (y3[0] * y3[0] + y3[1] * y3[1]) - y3[0] * (y2[0] * y2[0] + y2[1] * y2[1])) / D;
    center = x1 + y;

    //determine average value of radius squared
    y3 -= y;
    y2 -= y;
    double sqr_radius = nvis::inner(y, y) + nvis::inner(y2, y2) + nvis::inner(y3, y3);

    return sqr_radius / 3.;
}

template<typename T1, typename T2>
inline bool triangulation<T1, T2>::in_circle(const point_type& x, const point_type x1,
        const point_type& x2, const point_type& x3) const
{
    double radius2, dist2;
    point_type center;
    radius2 = circumcircle(center, x1, x2, x3);

    // check if inside/outside circumcircle
    dist2 = nvis::inner(x - center, x - center);
    return (dist2 < (0.999999999999*radius2));
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::check_edge(index_type p_id,
        index_type p1, index_type p2, index_type tri)
{
    index_type p3;
    index_type next;
    triangle_type tri_ids;

    const point_type& x  = get_vertex(p_id);
    const point_type& x1 = get_vertex(p1);
    const point_type& x2 = get_vertex(p2);
    next = get_edge_neighbor(tri, p1, p2);

    if (next >= 0) { //i.e., not a boundary edge
        // get neighbor info including opposite point
        tri_ids = get_triangle_vertices(next);
        for (int i = 0; i < 3; ++i) {
            if (tri_ids[i] != p1 && tri_ids[i] != p2) {
                p3 = tri_ids[i];
                break;
            }
        }
        const point_type& x3 = get_vertex(p3);

        // see whether point is in circumcircle
        if (in_circle(x3, x, x1, x2)) {// swap diagonal
            remove_tri_reference(p1, tri);
            remove_tri_reference(p2, next);
            __vert_to_tris[p_id].push_back(next);
            __vert_to_tris[p3].push_back(tri);
            __tri_to_verts[tri]  = triangle_type(p_id, p3, p2);
            __tri_to_verts[next] = triangle_type(p_id, p1, p3);

            // keep track of modified triangles
            __modified.insert(tri);
            __modified.insert(next);

            // two new edges become suspect
            check_edge(p_id, p3, p2, tri);
            check_edge(p_id, p1, p3, next);
        } //in circle
    } //interior edge
}

#undef DEL2D_TOLERANCE

template<typename T1, typename T2>
inline void triangulation<T1, T2>::shuffle(std::vector<index_type>& shuffled, size_t n) const
{
    assert(n > 0);
    shuffled.resize(n);
    for (int i = 0 ; i < n ; ++i) {
        shuffled[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(shuffled.begin(), shuffled.end(), g);
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::
barycentric(double b[3], const nvis::vec2 p[3], const nvis::vec2& x) const
{
    nvis::vec2 e0 = p[1] - p[0], e1 = p[2] - p[0];
    nvis::vec2 y = x - p[0];
    // b0*p0 + b1*p1 + b2*p2 = p
    // (1-b1-b2)*p0 + b1*p1 + b2*p2 = p
    // b1*(p1-p0) + b2*(p2-p0) = p-p0
    // b1*e0 + b2*e1 = y
    double delta = e0[0] * e1[1] - e0[1] * e1[0];
    b[1] = (y[0] * e1[1] - y[1] * e1[0]) / delta;
    b[2] = (e0[0] * y[1] - e0[1] * y[0]) / delta;
    b[0] = 1. - b[1] - b[2];
}

template<typename T1, typename T2>
inline void triangulation<T1, T2>::
set_export_mode(export_mode mode, size_t f) const
{
    __export_mode = (mode < 4 ? mode : NONE);
    __export_frequency = f; // controls export frequency when "SOME" is selected
    std::cerr << "export mode ";
    switch (__export_mode) {
        case NONE:
            std::cerr << "none";
            break;
        case EACH:
            std::cerr << "each frame";
            break;
        case ONE_PER_BATCH:
            std::cerr << "one frame per batch";
            break;
        case SOME:
            std::cerr << "every " << f << "th frame";
            break;
        default:
            std::cerr << "ERROR in triangulation<T1,T2>::set_export_mode()\n";
            __export_mode = NONE;
            return;
    }
    std::cerr << " has been selected\n";
}

}

#endif
