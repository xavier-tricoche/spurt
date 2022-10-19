#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>


// voxel stuff

typedef nvis::fixed_vector<int, 8>  voxel;
typedef nvis::lexicographical_order Lt_voxel;
Lt_voxel __Lt_voxel;
inline voxel create_voxel(int i, const nvis::ivec3& size)
{
    voxel v;
    v[0] = i;
    v[1] = i + 1;
    v[2] = v[1] + size[0];
    v[3] = v[0] + size[0];
    v[4] = v[0] + size[0] * size[1];
    v[5] = v[1] + size[0] * size[1];
    v[6] = v[2] + size[0] * size[1];
    v[7] = v[3] + size[0] * size[1];
    return v;
}

// face stuff

typedef nvis::fixed_vector<int, 4>  face;
typedef nvis::lexicographical_order Lt_face;
Lt_face __Lt_face;
const int voxel_faces[][4] = {
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {3, 0, 4, 7}
};
inline face create_face(voxel v, int i)
{
    face f;
    for (int n = 0 ; n < 4 ; ++n) {
        f[n] = v[voxel_faces[i][n]];
    }
    return f;
}

// edge stuff

typedef std::pair<int, int> edge;
inline edge create_edge(int a, int b)
{
    if (a <= b) {
        return edge(a, b);
    } else {
        return edge(b, a);
    }
}
struct Lt_edge {
    bool operator()(const edge& e0, const edge& e1) const {
        if (e0.first < e1.first) {
            return true;
        } else if (e0.first > e1.first) {
            return false;
        } else {
            return (e0.second < e1.second);
        }
    }
} __Lt_edge;
const int voxel_edges[][2] = {
    {0, 1}, {1, 2}, {2, 3}, {0, 3},
    {4, 5}, {5, 6}, {6, 7}, {4, 7},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};
inline edge create_edge(const voxel& v, int i)
{
    return create_edge(v[voxel_edges[i][0]], v[voxel_edges[i][1]]);
}
inline edge create_edge(const face& f, int i)
{
    return create_edge(f[i], f[(i+1)%4]);
}

// interface stuff

struct vertex {
    vertex() : type(-1), label("") {}
    vertex(const edge& _e, const std::vector<float>& vals, int edgeid = -1)
        : e(_e), type(0) {
        val.insert(vals[e.first]);
        val.insert(vals[e.second]);
        std::ostringstream os;
        os << "E" << edgeid;
        label = os.str();
    }
    vertex(const face& _f, const std::vector<float>& vals, int faceid = -1)
        : f(_f), type(1) {
        for (int i = 0 ; i < 4 ; ++i) {
            val.insert(vals[f[i]]);
        }
        std::ostringstream os;
        os << "F" << faceid;
        label = os.str();
    }
    vertex(const voxel& _v, const std::vector<float>& vals)
        : vx(_v), type(2) {
        for (int i = 0 ; i < 8 ; ++i) {
            val.insert(vals[vx[i]]);
        }
        label = "V";
    }
    vertex(const vertex& v)
        : e(v.e), f(v.f), vx(v.vx),
          val(v.val.begin(), v.val.end()), type(v.type), label(v.label) {}
          
    bool operator<(const vertex& v) const {
        if (v.type < type) {
            return false;
        }
        if (type < v.type) {
            return true;
        }
        if (type == 0) {
            return __Lt_edge(e, v.e);
        } else if (type == 1) {
            return __Lt_face(f, v.f);
        } else {
            return __Lt_voxel(vx, v.vx);
        }
    }
    
    bool operator==(const vertex& v) const {
        return (!((*this) < v) && !(v < (*this)));
    }
    
    bool operator!=(const vertex& v) const {
        return ((*this) < v || v < (*this));
    }
    
    bool match(float f) const {
        return val.find(f) != val.end();
    }
    
    bool match(float f0, float f1) const {
        return match(f0) && match(f1);
    }
    
    bool all_values_match(const vertex& v) const {
        std::set<float>::const_iterator it;
        for (it = val.begin() ; it != val.end() ; ++it) {
            if (!v.match(*it)) {
                return false;
            }
        }
        return true;
    }
    
    edge                e;
    face                f;
    voxel               vx;
    std::set<float>     val;
    int                 type;
    std::string         label;
};

std::ostream& operator<<(std::ostream& out, const vertex& v)
{
    out << v.label;
    return out;
}

typedef std::pair<vertex, vertex>   link;
inline link make_link(const vertex& v0, const vertex& v1)
{
    if (v1 < v0) {
        return link(v1, v0);
    } else {
        return link(v0, v1);
    }
}
struct Lt_link {
    bool operator()(const link& c0, const link& c1) {
        if (c0.first < c1.first) {
            return true;
        } else if (c1.first < c0.first) {
            return false;
        }
        return (c0.second < c1.second);
    }
};

template<typename T>
struct Lt_set {
    bool operator()(const std::set<T>& s0, const std::set<T>& s1) const {
        if (s0.size() < s1.size()) {
            return true;
        } else if (s1.size() < s0.size()) {
            return false;
        }
        std::set<T>::const_iterator it0, it1;
        for (it0 = s0.begin(), it1 = s1.begin() ; it0 != s0.end() && it1 != s1.end() ; ++it0, ++it1) {
            if (*it0 < *it1) {
                return true;
            } else if (*it1 < *it0) {
                return false;
            }
        }
        return false;
    }
};

template<typename T>
struct Eq_set {
    bool operator()(const std::set<T>& s0, const std::set<T>& s1) const {
        if (s0.size() != s1.size()) {
            return false;
        }
        std::set<T>::const_iterator it0, it1;
        for (it0 = s0.begin(), it1 = s1.begin() ; it0 != s0.end() && it1 != s1.end() ; ++it0, ++it1) {
            if (*it0 != *it1) {
                return false;
            }
        }
        return true;
    }
};

std::ostream& operator<<(std::ostream& out, const link& c)
{
    out << "( " << c.first << " - " << c.second << " )";
}

void process_face(std::vector<link>& links,
                  std::vector<vertex>& points,
                  const vertex& center)
{
    if (points.size() == 2) {
        links.push_back(make_link(points[0], points[1]));
    } else if (points.size() == 3) {
        for (int i = 0 ; i < points.size() ; ++i) {
            links.push_back(make_link(points[i], center));
        }
        points.push_back(center);
    } else if (points.size() == 4) {
        typedef std::set<float> set_type;
        std::set<set_type, Lt_set<float> > unique_sets;
        for (int i = 0 ; i < point.size() ; ++i) {
            unique_set.insert(points[i].val);
        }
        if (unique_set.size() == 2) {
            int minid;
            Eq_set<float> equal_set;
            Lt_set<float> less_set;
            for (minid = 0 ; minid < 4 ; ++minid) {
                if (equal_set(points[minid].val, unique_set.front())) {
                    break;
                }
            }
            if (less_set(points[(minid+1)%4], points[(minid+3)%4])) {
                links.push_back(make_link(points[minid], points[(minid+1)%3]));
                links.push_back(make_link(points[(minid+2)%3], points[(minid+3)%3]));
            } else {
                links.push_back(make_link(points[minid], points[(minid+3)%3]));
                links.push_back(make_link(points[(minid+2)%3], points[(minid+1)%3]));
            }
        } else if (unique_set.size() == 3) {
            // identify vertices sharing same value set
            if (equal_set(points[0].val, points[2].val)) {
                links.push_back(make_link(points[0], points[1]));
                links.push_back(make_link(points[2], points[3]));
            } else {
                links.push_back(make_link(points[0], points[3]));
                links.push_back(make_link(points[1], points[2]));
            }
        } else if (unique_set.size() == 4) {
            for (int i = 0 ; i < points.size() ; ++i) {
                links.push_back(make_link(points[i], center));
            }
            points.push_back(center);
        }
    }
}

struct triangle {
    triangle(const vertex& v0, const vertex& v1, const vertex& v2) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
    }
    
    vertex& operator[](unsigned int i) {
        return v[i];
    }
    
    const vertex& operator[](unsigned int i) const {
        return v[i];
    }
    
    vertex v[3];
};

void triangulate_loop(std::vector<triangle>& triangles, const std::vector<vertex>& loop)
{
    if (loop.size() == 3) {
        triangles.push_back(triangle(loop[0], loop[1], loop[2]));
    } else if (loop.size() == 4) {
        triangles.push_back(triangle(loop[0], loop[1], loop[2]));
        triangles.push_back(triangle(loop[0], loop[2], loop[3]));
    } else if (loop.size() >= 5) {
        std::vector<vertex> subloop0, subloop1;
        int cut = loop.size() / 2;
        std::copy(loop.begin(), loop.begin() + cut + 1, std::back_inserter(subloop0));
        triangulate_loop(triangles, subloop0);
        subloop1.push_back(loop[0]);
        std::copy(loop.begin() + cut, loop.end(), std::back_inserter(subloop1));
        triangulate_loop(triangles, subloop1);
    }
}

bool extend_loop(std::vector<vertex>& loop, const std::vector<link>& links,
                 std::vector<bool>& included)
{
    while (loop.front() != loop.back()) {
        std::cerr << "loop currently contains: ";
        std::copy(loop.begin(), loop.end(),
                  std::ostream_iterator<vertex>(std::cerr, " "));
        std::cerr << '\n';
        std::cerr << "links are: ";
        std::copy(links.begin(), links.end(),
                  std::ostream_iterator<link>(std::cerr, " "));
        std::cerr << '\n';
        std::cerr << "included flags are: ";
        std::copy(included.begin(), included.end(),
                  std::ostream_iterator<bool>(std::cerr, " "));
        std::cerr << '\n';
        bool found = false;
        for (int i = 0 ; !found && i < links.size() ; ++i) {
            const link& c = links[i];
            if (included[i]) {
                continue;
            }
            if (c.first == loop.back()) {
                loop.push_back(c.second);
                included[i] = true;
                found = true;
            } else if (c.second == loop.back()) {
                loop.push_back(c.first);
                included[i] = true;
                found = true;
            }
        }
        
        if (!found) {
            break;
        }
    }
    
    return (loop.front() == loop.back());
}

template<typename T>
struct Lt_pair {
    bool operator()(const std::pair<T, T>& p0, const std::pair<T, T>& p1) {
        if (p0.first < p1.first) {
            return true;
        }
        if (p1.first < p0.first) {
            return false;
        }
        return p0.second < p1.second;
    }
};

bool triangulate(std::vector<triangle>& triangles,
                 const std::vector<link>& links)
{
    typedef std::pair<float, float>                     Key;
    typedef std::pair<Key, int>                         value_type;
    typedef std::multimap<Key, int, Lt_pair<float> >    multimap;
    typedef multimap::iterator                          iterator;
    typedef std::pair<iterator, iterator>               range;
    
    multimap all_edge_values;
    
    // create list of unique matching between link and pair of values
    for (int i = 0 ; i < links.size() ; ++i) {
        const link& c = links[i];
        std::vector<float> overlap;
        std::set_intersection(c.first.val.begin(), c.first.val.end(),
                              c.second.val.begin(), c.second.val.end(),
                              std::back_inserter(overlap));
        for (int j = 0 ; j < overlap.size() - 1 ; ++j) {
            for (int k = j + 1 ; k < overlap.size() ; ++k) {
                Key p(overlap[j], overlap[k]);
                all_edge_values.insert(value_type(p, i));
            }
        }
    }
    
    // loop over found value pairs
    for (iterator it = all_edge_values.begin() ; it != all_edge_values.end() ;
            it = all_edge_values.upper_bound(it->first)) {
        std::vector<link> __links;
        range r = all_edge_values.equal_range(it->first);
        for (iterator __it = r.first ; __it != r.second ; ++__it) {
            __links.push_back(links[__it->second]);
        }
        std::cerr << "links for Key (" << it->first.first << ", " << it->first.second
                  << ") are:\n";
        for (int i = 0 ; i < __links.size() ; ++i) {
            std::cerr << __links[i] << '\n';
        }
        if (__links.size() == 1) {
            continue;
        } else if (__links.size() < 3) {
            std::cerr << "invalid list of links: only " << __links.size()
                      << " links found for Key (" << it->first.first << ", " << it->first.second
                      << ")\n";
            throw std::runtime_error("giving up");
        };
        
        // sanity check
        std::map<vertex, int> vertex_counter;
        for (int i = 0 ; i < __links.size() ; ++i) {
            std::map<vertex, int>::iterator it;
            it = vertex_counter.find(__links[i].first);
            if (it != vertex_counter.end()) {
                ++(it->second);
            } else {
                vertex_counter[__links[i].first] = 1;
            }
            it = vertex_counter.find(__links[i].second);
            if (it != vertex_counter.end()) {
                ++(it->second);
            } else {
                vertex_counter[__links[i].second] = 1;
            }
        }
        for (std::map<vertex, int>::iterator it = vertex_counter.begin() ;
                it != vertex_counter.end() ; ++it) {
            if (it->second > 2) {
                std::cerr << "link list contains redundant vertex\n";
                return false;
            }
        }
        
        std::vector<bool> included(__links.size(), false);
        for (int i = 0 ; i < __links.size() ; ++i) {
            if (included[i]) {
                continue;
            }
            std::vector<vertex> loop;
            loop.push_back(__links[i].first);
            loop.push_back(__links[i].second);
            included[i] = true;
            if (!extend_loop(loop, __links, included)) {
                return false;
            }
            
            // triangulate loop
            triangulate_loop(triangles, loop);
        }
    }
    
    return true;
}

inline nvis::fvec3 coord(int i, const nvis::ivec3& size, const nvis::vec3& step)
{
    int u, v, w, j = i / size[0];
    u = i % size[0];
    v = j % size[1];
    w = j / size[1];
    return nvis::fvec3((float)u*step[0], (float)v*step[1], (float)w*step[2]);
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "Multi-material marching cubes (Wu and Sullivan, 2003)\n";
        std::cerr << "USAGE: " << argv[0] << " <input.nrrd> <output.vtk>\n";
        return -1;
    }
    
    Nrrd* nin = xavier::readNrrd(argv[1]);
    std::vector<float> values;
    xavier::to_vector(values, nin);
    nvis::ivec3 size;
    nvis::vec3 step;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i] = nin->axis[i].size;
        step[i] = nin->axis[i].spacing;
    }
    
    int nbvoxels = (size[0] - 1) * (size[1] - 1) * (size[2] - 1);
    
    std::set<vertex>        vertices;
    std::vector<triangle>   triangles;
    
    // find all intersected edges
    for (int n = 0 ; n < nbvoxels ; ++n) {
        std::cerr << n*100 / nbvoxels << "% completed, "
                  "so far " << vertices.size() << " vertices and " << triangles.size() << " triangles  \r";
        voxel v = create_voxel(n, size);
        
        std::cerr << "processing voxel: " << v << '\n';
        std::cerr << "corresponding labels are ";
        for (int i = 0 ; i < 8 ; ++i) {
            std::cerr << values[v[i]] << " ";
        }
        std::cerr << '\n';
        
        std::vector<link>   links;
        
        // loop over edges, add edge vertices
        
        for (int i = 0 ; i < 12 ; ++i) {
            edge e = create_edge(v, i);
            float v0 = values[e.first];
            float v1 = values[e.second];
            if (v0 != v1) {
                std::cerr << "E" << i << " (" << v0 << ", " << v1 << ")\n";
                vertex v(e, values, i);
                vertices.insert(v);
            }
        }
        
        // loop over faces, add face vertices
        
        std::vector<vertex> face_points;
        for (int i = 0 ; i < 6 ; ++i) {
            std::set<float> local_values;
            std::vector<vertex> face_edge_points;
            face f = create_face(v, i);
            std::set<vertex>::const_iterator it;
            
            // loop over face edges
            for (int j = 0 ; j < 4 ; ++j) {
                vertex tmp(create_edge(f, j), values);
                it = vertices.find(tmp);
                if (it != vertices.end()) {
                    local_values.insert(it->val.begin(), it->val.end());
                    face_edge_points.push_back(*it);
                }
            }
            
            // check the results
            if (local_values.size() == 2) {
                if (face_edge_points.size() == 2) {
                    links.push_back(make_link(face_edge_points[0],
                                              face_edge_points[1]));
                } else {
                    // arbitrary but globally consistent connectivity
                    int minpos = std::distance(&f[0], std::min_element(&f[0], &f[4]));
                    links.push_back(make_link(face_edge_points[minpos],
                                              face_edge_points[(minpos+1)%4]));
                    links.push_back(make_link(face_edge_points[(minpos+2)%4],
                                              face_edge_points[(minpos+3)%4]));
                }
            } else if (local_values.size() >= 3) {
                vertex fp(f, values, i);
                // 4 values -> face point
                // 3 values, 3 points -> face point
                // 3 values, 4 points -> no face point
                if (!(fp.val.size() == 3 && local_values.size() == 4)) {
                    vertices.insert(fp);
                    face_points.push_back(fp);
                    std::cerr << "F" << i << " (";
                    std::copy(fp.val.begin(), fp.val.end(), std::ostream_iterator<float>(std::cerr, " "));
                    std::cerr << ")\n";
                    for (int j = 0 ; j < face_edge_points.size() ; ++j) {
                        links.push_back(make_link(face_edge_points[j], fp));
                    }
                } else {
                    if (face_edge_points[0].all_values_match(face_edge_points[1])) {
                        links.push_back(make_link(face_edge_points[0], face_edge_points[1]));
                        links.push_back(make_link(face_edge_points[2], face_edge_points[3]));
                    } else {
                        links.push_back(make_link(face_edge_points[0], face_edge_points[3]));
                        links.push_back(make_link(face_edge_points[1], face_edge_points[2]));
                    }
                }
            }
        }
        
        std::set<link, Lt_link> unique_links(links.begin(), links.end());
        links.clear();
        std::copy(unique_links.begin(), unique_links.end(),
                  std::back_inserter(links));
                  
        std::cerr << "available links for this voxel:\n";
        for (int k = 0 ; k < links.size() ; ++k) {
            std::cerr << links[k] << '\n';
        }
        
        // form closed loops
        if (face_points.size() == 2) {
            links.push_back(make_link(face_points[0], face_points[1]));
            if (!triangulate(triangles, links)) {
                std::cerr << "triangulation failed with 2 face vertices\n";
                throw std::runtime_error("giving up");
            }
        } else if (face_points.size() >= 3) {
            // first connect all face points that match exactly
            for (int j = 0 ; j + 1 < face_points.size() ; ++j) {
                const vertex& fp0 = face_points[j];
                for (int k = j + 1 ; k < face_points.size() ; ++k) {
                    const vertex& fp1 = face_points[k];
                    if (fp0.all_values_match(fp1)) {
                        links.push_back(make_link(fp0, fp1));
                    }
                }
            }
            
            std::vector<triangle> tmp_triangles;
            bool done = triangulate(tmp_triangles, links);
            if (!done) {
                tmp_triangles.clear();
                vertex vp(v, values);
                std::cerr << "V (";
                std::copy(vp.val.begin(), vp.val.end(), std::ostream_iterator<float>(std::cerr, " "));
                std::cerr << ")\n";
                vertices.insert(vp);
                for (int i = 0 ; i < face_points.size() ; ++i) {
                    links.push_back(make_link(vp, face_points[i]));
                }
                if (!triangulate(tmp_triangles, links)) {
                    std::cerr << "triangulation failed with " << face_points.size() << " face vertices after adding a voxel vertex\n";
                    throw std::runtime_error("giving up");
                } else {
                    std::cerr << "triangulation fixed with " << face_points.size() << " face vertices after adding a voxel vertex\n";
                }
            }
            std::copy(tmp_triangles.begin(), tmp_triangles.end(), std::back_inserter(triangles));
        } else {
            if (!triangulate(triangles, links)) {
                std::cerr << "triangulation failed without face vertex\n";
                throw std::runtime_error("giving up");
            }
        }
    }
    
    std::cerr << '\n';
    std::cerr << vertices.size() << " vertices in triangulation\n";
    std::cerr << triangles.size() << " triangles in triangulation\n";
    
    std::map<vertex, int>       pos_ids;
    std::vector<nvis::fvec3>    pos;
    for (std::set<vertex>::const_iterator it = vertices.begin() ; it != vertices.end() ; ++it) {
        const vertex& v = *it;
        nvis::fvec3 x;
        if (v.type == 0) {
            x = 0.5 * (coord(v.e.first, size, step) + coord(v.e.second, size, step));
        } else if (v.type == 1) {
            x = 0.5 * (coord(v.f[0], size, step) + coord(v.f[2], size, step));
        } else {
            x = 0.5 * (coord(v.vx[0], size, step) + coord(v.vx[6], size, step));
        }
        pos.push_back(x);
        pos_ids[v] = pos.size() - 1;
    }
    
    std::fstream vtk(argv[2], std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << "Generated through multi-material marching cubes from file " << argv[1] << '\n'
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << pos.size() << " float\n";
    for (int i = 0 ; i < pos.size() ; ++i) {
        vtk << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
    }
    vtk << "CELLS " << triangles.size() << " " << 4*triangles.size() << '\n';
    for (int i = 0 ; i < triangles.size() ; ++i) {
        vtk << "3";
        for (int j = 0 ; j < 3 ; ++j) {
            vtk << " " << pos_ids[triangles[i][j]];
        }
        vtk << '\n';
    }
    vtk << "CELL_TYPES " << triangles.size() << '\n';
    for (int i = 0 ; i < triangles.size() ; ++i) {
        vtk << "5\n";
    }
    vtk << "POINT_DATA " << pos.size() << '\n';
    vtk << "SCALARS nb_materials float 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (std::set<vertex>::const_iterator it = vertices.begin() ; it != vertices.end() ; ++it) {
        const vertex& v = *it;
        vtk << it->val.size() << '\n';
    }
    vtk.close();
    
    return 0;
}

























































































































