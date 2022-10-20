#include <vector>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>
#include <teem/nrrd.h>

#define QUIET

char* outs, *file;
float alpha;
int niter;
bool verbose;
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
    hestOptAdd(&hopt, "i",      "file",             airTypeString,  1, 1, &file,    NULL,       "input NRRD file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output file name");
    hestOptAdd(&hopt, "n",      "niter",            airTypeInt,     0, 1, &niter,   "10",       "number of iterations of Laplacian smoothing");
    hestOptAdd(&hopt, "a",      "alpha",            airTypeFloat,   0, 1, &alpha,   "0.5",      "weight coefficient of Laplacian smoothing");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &verbose, "0",        "verbose mode (debugging)");
    
    _hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                    me, "Extract multi-material boundaries from raster volume and export smooth geometry",
                    AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline nvis::ivec3 int_to_ivec(int i, const nvis::ivec3& s)
{
    int x = i % s[0];
    int aux = i / s[0];
    int y = aux % s[1];
    int z = aux / s[1];
    return nvis::ivec3(x, y, z);
}

inline int ivec_to_int(const nvis::ivec3& c, const nvis::ivec3& s)
{
    return c[0] + s[0]*(c[1] + s[1]*c[2]);
}

inline int vert_to_vox(int i, const nvis::ivec3& s)
{
    nvis::ivec3 svx = s - nvis::ivec3(1, 1, 1);
    return ivec_to_int(int_to_ivec(i, s), svx);
}

inline int vox_to_vert(int i, const nvis::ivec3& s)
{
    nvis::ivec3 svx = s - nvis::ivec3(1, 1, 1);
    return ivec_to_int(int_to_ivec(i, svx), s);
}

// voxel stuff

typedef nvis::fixed_vector<int, 8>  voxel;
typedef nvis::lexicographical_order Lt_voxel;
Lt_voxel __Lt_voxel;
inline voxel create_voxel(int __i, const nvis::ivec3& size)
{
    int i = vox_to_vert(__i, size);
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
inline void voxel_data(std::set<float>& r, const voxel& v, const std::vector<float>& data)
{
    for (int i = 0 ; i < 8 ; ++i) {
        r.insert(data[v[i]]);
    }
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
    face f, sf;
    for (int n = 0 ; n < 4 ; ++n) {
        f[n] = v[voxel_faces[i][n]];
    }
    
    // re-order coefficients
    int min = std::distance(&f[0], std::min_element(&f[0], &f[4]));
    if (f[(min+3)%4] < f[(min+1)%4])
        // backward loop
        for (int i = 0 ; i < 4 ; ++i) {
            sf[i] = f[(min+4-i)%4];
        }
    else
        // forward loop
        for (int i = 0 ; i < 4 ; ++i) {
            sf[i] = f[(min+i)%4];
        }
        
    return sf;
}

// edge stuff

typedef std::pair<int, int> edge;
template<typename T>
inline std::pair<T, T> create_pair(T a, T b)
{
    if (a <= b) {
        return std::make_pair(a, b);
    } else {
        return std::make_pair(b, a);
    }
}
template<typename T>
struct Lt_pair {
    bool operator()(const std::pair<T, T>& e0, const std::pair<T, T>& e1) const {
        if (e0.first < e1.first) {
            return true;
        } else if (e0.first > e1.first) {
            return false;
        } else {
            return (e0.second < e1.second);
        }
    }
};
Lt_pair<int> __Lt_edge;
template<typename T>
bool operator==(const std::pair<T, T>& e0, const std::pair<T, T>& e1)
{
    return e0.first == e1.first &&
           e0.second == e1.second;
}
template<typename >T
std::ostream& operator<<(std::ostream& out, const std::pair<T, T>& e)
{
    out << "(" << e.first << ", " << e.second << ")";
    return out;
}

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

inline std::string to_str(const vertex& v)
{
    std::ostringstream os;
    switch (v.type) {
        case 0:
            os << "edge vertex " << v.e;
            break;
        case 1:
            os << "face vertex " << v.f;
            break;
        case 2:
            os << "voxel vertex " << v.vx;
            break;
    }
    return os.str();
}

typedef std::pair<vertex, vertex>   t_edge;

template<typename T>
struct Lt_set {
    bool operator()(const std::set<T>& s0, const std::set<T>& s1) const {
        if (s0.size() < s1.size()) {
            return true;
        } else if (s1.size() < s0.size()) {
            return false;
        }
        typename std::set<T>::const_iterator it0, it1;
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
        typename std::set<T>::const_iterator it0, it1;
        for (it0 = s0.begin(), it1 = s1.begin() ; it0 != s0.end() && it1 != s1.end() ; ++it0, ++it1) {
            if (*it0 != *it1) {
                return false;
            }
        }
        return true;
    }
};

std::ostream& operator<<(std::ostream& out, const t_edge& c)
{
    out << "( " << c.first << " - " << c.second << " )";
}

inline std::string to_str(const t_edge& c)
{
    std::ostringstream os;
    os << "(" << to_str(c.first) << " - " << to_str(c.second) << ")";
    return os.str();
}

bool process_face(std::vector<t_edge>& t_edges,
                  std::vector<vertex>& points,
                  const vertex& center)
{
    Eq_set<float> equal_set;
    Lt_set<float> less_set;
    if (!points.size()) {
        return false;
    }
    
    bool included_center = false;
    
    std::cerr << "process face: " << points.size() << " vertices in input\n";
    for (int i = 0 ; i < points.size() ; ++i) {
        std::cerr << i << ": " << points[i] << '\n';
    }
    if (points.size() == 2) {
        std::cerr << "2 points on face - 1 edge\n";
        t_edges.push_back(make_t_edge(points[0], points[1]));
    } else if (points.size() == 3) {
        std::cerr << "3 points on face - 3 edges\n";
        for (int i = 0 ; i < points.size() ; ++i) {
            t_edges.push_back(make_t_edge(points[i], center));
        }
        points.push_back(center);
        included_center = true;
    } else if (points.size() == 4) {
        typedef std::set<float> set_type;
        std::set<set_type, Lt_set<float> > unique_sets;
        for (int i = 0 ; i < points.size() ; ++i) {
            unique_sets.insert(points[i].val);
        }
        std::cerr << unique_sets.size() << " unique sets of values\n";
        if (unique_sets.size() == 1) {
            // simply connect the min vertex to its smallest neighbor
            int minid = std::distance(points.begin(), std::min_element(points.begin(), points.end()));
            if (points[(minid+1)%4] < points[(minid+3)%4]) {
                t_edges.push_back(make_t_edge(points[minid], points[(minid+1)%4]));
                t_edges.push_back(make_t_edge(points[(minid+2)%4], points[(minid+3)%4]));
            } else {
                t_edges.push_back(make_t_edge(points[minid], points[(minid+3)%4]));
                t_edges.push_back(make_t_edge(points[(minid+2)%4], points[(minid+1)%4]));
            }
            std::cerr << "4 points on face, 1 unique set of values - 2 edges\n";
        } else if (unique_sets.size() == 2) {
            // identify vertices sharing same value set and connect them
            if (equal_set(points[0].val, points[2].val)) {
                t_edges.push_back(make_t_edge(points[0], points[2]));
                t_edges.push_back(make_t_edge(points[1], points[3]));
            } else {
                t_edges.push_back(make_t_edge(points[0], points[1]));
                t_edges.push_back(make_t_edge(points[2], points[3]));
            }
            std::cerr << "4 points on face, 2 unique sets of values - 2 edges\n";
        } else if (unique_sets.size() >= 3) {
            for (int i = 0 ; i < points.size() ; ++i) {
                t_edges.push_back(make_t_edge(points[i], center));
            }
            points.push_back(center);
            included_center = true;
            std::cerr << "4 points on face, " << unique_sets.size() << " unique sets of values - 4 edges\n";
        }
    }
    
    return included_center;
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

namespace face_query {
const int faces[][4] = {
    {0, 1, 3, 5},
    {0, 1, 2, 4},
    {2, 3, 4, 5}
};

void filter(std::vector<vertex>& selected,
            const std::vector<vertex>& vertices,
            const voxel& v,
            int faceid)
{
    if (faceid < 6) {
        face f = create_face(v, faceid);
        for (int i = 0 ; i < vertices.size() ; ++i) {
            if (vertices[i].type != 0) {
                continue;
            }
            for (int j = 0 ; j < 4 ; ++j) {
                edge e = create_edge(f[j], (f[(j+1)%4]));
                if (vertices[i].e == e) {
                    selected.push_back(vertices[i]);
                    break;
                }
            }
        }
    } else {
        for (int i = 0 ; i < vertices.size() ; ++i) {
            if (vertices[i].type != 1) {
                continue;
            }
            for (int j = 0 ; j < 4 ; ++j) {
                face f = create_face(v, faces[faceid-6][j]);
                if (nvis::all(vertices[i].f == f)) {
                    selected.push_back(vertices[i]);
                    break;
                }
            }
        }
    }
}
}

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

bool extend_loop(std::vector<vertex>& loop, const std::vector<t_edge>& t_edges,
                 std::vector<bool>& included)
{
    while (loop.front() != loop.back()) {
        std::cerr << "loop currently contains: ";
        std::copy(loop.begin(), loop.end(),
                  std::ostream_iterator<vertex>(std::cerr, " "));
        std::cerr << '\n';
        std::cerr << "t_edges are: ";
        std::copy(t_edges.begin(), t_edges.end(),
                  std::ostream_iterator<t_edge>(std::cerr, " "));
        std::cerr << '\n';
        std::cerr << "included flags are: ";
        std::copy(included.begin(), included.end(),
                  std::ostream_iterator<bool>(std::cerr, " "));
        std::cerr << '\n';
        bool found = false;
        for (int i = 0 ; !found && i < t_edges.size() ; ++i) {
            const t_edge& c = t_edges[i];
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

bool triangulate(std::vector<triangle>& triangles,
                 const std::vector<t_edge>& t_edges,
                 const std::vector<float>& corner_values)
{
    typedef std::pair<float, float>                     Key;
    typedef std::set<Key, Lt_pair<float> >              key_set;
    typedef std::pair<Key, int>                         value_type;
    typedef std::multimap<Key, int, Lt_pair<float> >    multimap;
    typedef multimap::iterator                          iterator;
    typedef std::pair<iterator, iterator>               range;
    
    multimap all_edge_values;
    
    key_set all_keys;
    for (int i = 0 ; i < 12 ; ++i) {
        float v0 = corner_values[voxel_edges[i][0]];
        float v1 = corner_values[voxel_edges[i][1]];
        if (v0 == v1) {
            continue;
        }
        all_keys.insert(Key(std::min(v0, v1), std::max(v0, v1)));
    }
    std::cerr << "all valid keys for this voxel are ";
    for (key_set::iterator it = all_keys.begin() ; it != all_keys.end() ; ++it) {
        std::cerr << "(" << it->first << ", " << it->second << ") ";
    }
    std::cerr << '\n';
    
    // identify matching between edges and value pairs
    for (key_set::const_iterator it = all_keys.begin() ; it != all_keys.end() ; ++it) {
        float v0 = it->first;
        float v1 = it->second;
        for (int i = 0 ; i < t_edges.size() ; ++i) {
            const t_edge& e = t_edges[i];
            if (e.first.match(v0, v1) && e.second.match(v0, v1)) {
                all_edge_values.insert(value_type(*it, i));
                std::cerr << "edge " << e.first << " - " << e.second << " matches key (" << v0 << ", " << v1 << ")\n";
            }
        }
    }
    
    // loop over found value pairs
    for (iterator it = all_edge_values.begin() ; it != all_edge_values.end() ;
            it = all_edge_values.upper_bound(it->first)) {
        std::vector<t_edge> __t_edges;
        range r = all_edge_values.equal_range(it->first);
        for (iterator __it = r.first ; __it != r.second ; ++__it) {
            __t_edges.push_back(t_edges[__it->second]);
        }
        std::cerr << "t_edges for Key (" << it->first.first << ", " << it->first.second
                  << ") are:\n";
        for (int i = 0 ; i < __t_edges.size() ; ++i) {
            std::cerr << __t_edges[i] << '\n';
        }
        if (__t_edges.size() == 1) {
            continue;
        } else if (__t_edges.size() < 3) {
            std::cerr << "invalid list of t_edges: only " << __t_edges.size()
                      << " t_edges found for Key (" << it->first.first << ", " << it->first.second
                      << ")\n";
            throw std::runtime_error("giving up");
        };
        
        // sanity check
        std::map<vertex, int> vertex_counter;
        for (int i = 0 ; i < __t_edges.size() ; ++i) {
            std::map<vertex, int>::iterator it;
            it = vertex_counter.find(__t_edges[i].first);
            if (it != vertex_counter.end()) {
                ++(it->second);
            } else {
                vertex_counter[__t_edges[i].first] = 1;
            }
            it = vertex_counter.find(__t_edges[i].second);
            if (it != vertex_counter.end()) {
                ++(it->second);
            } else {
                vertex_counter[__t_edges[i].second] = 1;
            }
        }
        for (std::map<vertex, int>::iterator it = vertex_counter.begin() ;
                it != vertex_counter.end() ; ++it) {
            if (it->second > 2) {
                std::cerr << "t_edge list contains redundant vertex\n";
                return false;
            }
        }
        
        std::vector<bool> included(__t_edges.size(), false);
        for (int i = 0 ; i < __t_edges.size() ; ++i) {
            if (included[i]) {
                continue;
            }
            std::vector<vertex> loop;
            loop.push_back(__t_edges[i].first);
            loop.push_back(__t_edges[i].second);
            included[i] = true;
            if (!extend_loop(loop, __t_edges, included)) {
                std::cerr << "unable to complete the loop\n";
                throw std::runtime_error("giving up");
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
    nvis::vec3 c = int_to_ivec(i, size);
    return c*step;
}

inline nvis::vec3 coord(const vertex& v, const nvis::ivec3& size, const nvis::vec3& step)
{
    nvis::vec3 r;
    if (v.type == 0) {
        r = 0.5 * (coord(v.e.first, size, step) + coord(v.e.second, size, step));
    } else if (v.type == 1) {
        r = 0.5 * (coord(v.f[0], size, step) + coord(v.f[2], size, step));
    } else {
        r = 0.5 * (coord(v.vx[0], size, step) + coord(v.vx[6], size, step));
    }
    return r;
}

void classify(std::vector<int>& simple, std::vector<int>& edge,
              std::vector<int>& corner, const std::vector<float>& values,
              const nvis::ivec3& size)
{
    nvis::ivec3 vsize = size - nvis::ivec3(1, 1, 1);
    int nbvoxels = vsize[0] * vsize[1] * vsize[2];
    std::vector<int> tmp;
    for (int i = 0 ; i < nbvoxels ; ++i) {
        std::set<float> data;
        voxel v = create_voxel(i, size);
        voxel_data(data, v, values);
        if (data.size() == 2) {
            simple.push_back(i);
        } else if (data.size() >= 3) {
            tmp.push_back(i);
        }
    }
    
    // distinguish between edge and corner voxels
    for (int n = 0 ; n < tmp.size() ; ++n) {
        std::set<float> d;
        voxel v = create_voxel(tmp[n], size);
        voxel_data(d, v, values);
        nvis::ivec3 vc = int_to_ivec(tmp[n], vsize);
        int i = vc[0], j = vc[1], k = vc[2];
        int counter = 0;
        for (int ii = i - 1 ; counter < 3 && ii <= i + 1 ; ++ii) {
            if (ii < 0 || ii >= vsize[0]) {
                continue;
            }
            for (int jj = j - 1 ; counter < 3 && jj <= j + 1 ; ++jj) {
                if (jj < 0 || jj >= vsize[1]) {
                    continue;
                }
                for (int kk = k - 1 ; counter < 3 && kk <= k + 1 ; ++kk) {
                    if (ii == i && jj == j && kk == k) {
                        continue;
                    }
                    if (kk < 0 || kk >= vsize[2]) {
                        continue;
                    }
                    int id = ivec_to_int(nvis::ivec3(ii, jj, kk), vsize);
                    voxel vv = create_voxel(id, size);
                    std::set<float> dd;
                    voxel_data(dd, vv, values);
                    std::vector<float> all;
                    std::set_union(d.begin(), d.end(), dd.begin(), dd.end(),
                                   std::back_inserter(all));
                    if (all.size() == dd.size()) {
                        ++counter;
                    }
                }
            }
        }
        
        if (counter == 2) {
            edge.push_back(tmp[n]);
        } else {
            corner.push_back(tmp[n]);
        }
    }
}

inline void vertex_to_set(std::set<int>& out, const vertex& in)
{
    out.clear();
    switch (in.type) {
        case 0:
            out.insert(in.e.first);
            out.insert(in.e.second);
            break;
        case 1: {
            std::vector<int> tmp(&in.f[0], &in.f[4]);
            std::sort(tmp.begin(), tmp.end());
            out.insert(tmp.begin(), tmp.end());
            break;
        }
        case 2: {
            std::vector<int> tmp(&in.vx[0], &in.vx[8]);
            std::sort(tmp.begin(), tmp.end());
            out.insert(tmp.begin(), tmp.end());
            break;
        }
    }
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
    os << "(";
    for (typename std::set<T>::const_iterator it = s.begin() ; it != s.end() ; ++it) {
        os << *it << ", ";
    }
    os << ")";
    return os;
}

void topology(std::vector<vertex>& corners, std::vector<t_edge>& edges,
              const std::vector<t_edge>& all_edges)
{
    std::map<vertex, int> counter;
    std::map<vertex, int>::iterator it;
    for (int i = 0 ; i < all_edges.size() ; ++i) {
        const vertex& v0 = all_edges[i].first;
        const vertex& v1 = all_edges[i].second;
        if (v0.val.size() > 2 && v1.val.size() > 2) {
            edges.push_back(all_edges[i]);
            it = counter.find(v0);
            if (it != counter.end()) {
                ++(it->second);
            } else {
                counter[v0] = 1;
            }
            it = counter.find(v1);
            if (it != counter.end()) {
                ++(it->second);
            } else {
                counter[v1] = 1;
            }
        }
    }
    for (it = counter.begin() ; it != counter.end() ; ++it) {
        if (it->second > 2) {
            corners.push_back(it->first);
        }
    }
    
    // sanity check - verify that no corner vertex has exactly 2 edge neighbors
    if (false) {
        typedef std::set<int> alt_face;
        std::set<alt_face, Lt_set<int> > neigh;
        
        std::map<vertex, std::list<int> > incident;
        std::map<vertex, std::list<int> >::iterator it;
        for (int i = 0 ; i < corners.size() ; ++i) {
            incident[corners[i]] = std::list<int>();
        }
        for (int i = 0 ; i < edges.size() ; ++i) {
            it = incident.find(edges[i].first);
            if (it != incident.end()) {
                it->second.push_back(i);
            }
            it = incident.find(edges[i].second);
            if (it != incident.end()) {
                it->second.push_back(i);
            }
        }
        for (it = incident.begin() ; it != incident.end() ; ++it) {
            neigh.clear();
            std::cout << "corner vertex " << to_str(it->first) << " is neighbored by vertices:";
            for (std::list<int>::iterator j = it->second.begin() ; j != it->second.end() ; ++j) {
                const vertex& v0 = edges[*j].first;
                const vertex& v1 = edges[*j].second;
                std::set<int> out;
                if (it->first == v0) {
                    std::cout << " " << to_str(v1);
                    vertex_to_set(out, v1);
                    neigh.insert(out);
                } else if (it->first == v1) {
                    vertex_to_set(out, v0);
                    std::cout << " " << to_str(v0);
                    neigh.insert(out);
                }
            }
            std::cout << '\n';
            std::cout << "of those neighbors, following distinct sets were identified:";
            for (std::set<alt_face, Lt_set<int> >::iterator it = neigh.begin() ; it != neigh.end() ; ++it) {
                std::cout << " " << *it;
            }
            std::cout << "\n\n";
        }
    }
}

typedef nvis::ivec3 tri_face;

void smooth_mesh(std::vector<nvis::vec3>& pos, const std::vector<tri_face>& faces,
                 const std::vector<int>& corners, double alpha, int niter)
{
    // identify constraints
    std::vector<bool> fixed(pos.size(), false);
    for (int i = 0 ; i < corners.size() ; ++i) {
        fixed[corners[i]] = true;
    }
    
    // determine connectivity
    typedef std::set<int>       neighbor;
    typedef neighbor::iterator  iterator;
    std::vector<neighbor> neigh(pos.size());
    for (int i = 0 ; i < faces.size() ; ++i) {
        const tri_face& tri = faces[i];
        neigh[tri[0]].insert(tri[1]);
        neigh[tri[0]].insert(tri[2]);
        neigh[tri[1]].insert(tri[2]);
        neigh[tri[1]].insert(tri[0]);
        neigh[tri[2]].insert(tri[0]);
        neigh[tri[2]].insert(tri[1]);
    }
    
    // loop over laplacian smoothing steps
    std::vector<nvis::vec3> next(pos.size());
    for (int n = 0 ; n < niter ; ++n) {
        for (int i = 0 ; i < pos.size() ; ++i) {
            if (fixed[i] || !neigh[i].size()) {
                next[i] = pos[i];
                continue;
            }
            nvis::vec3 avg(0, 0, 0);
            for (iterator it = neigh[i].begin() ; it != neigh[i].end() ; ++it) {
                avg += pos[*it];
            }
            avg /= (double)neigh[i].size();
            next[i] = (1. - alpha) * pos[i] + alpha * avg;
        }
        std::swap(next, pos);
    }
}


struct Lt_close_enough {
    bool operator()(const nvis::vec3& x, const nvis::vec3& y) const {
        double d = nvis::norm(x - y);
        if (d < 1.0e-8) {
            return false;
        } else {
            return order(x, y);
        }
    }
    
    nvis::lexicographical_order order;
};

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    std::ofstream blah;
    std::streambuf* init_buf = std::cerr.rdbuf();
    std::streambuf* alt_buf;
    if (!verbose) {
        blah.open("/tmp/alt-cerr.txt");
        alt_buf = blah.rdbuf();
        std::cerr.rdbuf(alt_buf);
    }
    
    Nrrd* nin = spurt::readNrrd(file);
    std::vector<float> values;
    spurt::to_vector(values, nin);
    nvis::ivec3 size;
    nvis::vec3 step;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i] = nin->axis[i].size;
        step[i] = nin->axis[i].spacing;
        if (!step[i] || std::isnan(step[i]) || std::isinf(step[i])) {
            step[i] = 1;
        }
    }
    
    int nbvoxels = (size[0] - 1) * (size[1] - 1) * (size[2] - 1);
    
    std::set<vertex>                    vertices;
    std::set<t_edge, Lt_pair<vertex> >  edges;
    std::vector<triangle>               triangles;
    
    int last_pct = -1;
    for (int n = 0 ; n < nbvoxels ; ++n) {
        int pct = n * 100 / nbvoxels;
        if (pct > last_pct) {
            std::cout << pct << "% completed, "
                      << "so far " << vertices.size() << " vertices and " << triangles.size() << " triangles  \r"
                      << std::flush;
            last_pct = pct;
        }
        voxel v = create_voxel(n, size);
        std::vector<float> voxel_values(8);
        for (int i = 0 ; i < 8 ; ++i) {
            voxel_values[i] = values[v[i]];
        }
        
        std::set <
        
        if (v[0] == 1017) {
            std::cerr.rdbuf(init_buf);
        } else {
            std::cerr.rdbuf(alt_buf);
        }
        
        std::cerr << "processing voxel: " << v << '\n';
        std::cerr << "corresponding labels are ";
        for (int i = 0 ; i < 8 ; ++i) {
            std::cerr << voxel_values[i] << " ";
        }
        std::cerr << '\n';
        
        std::set<t_edge, Lt_pair<vertex> >  local_edges;
        std::set<vertex>                    local_vertices;
        std::set<vertex>                    face_vertices;
        
        // loop over edges, add edge vertices
        
        for (int i = 0 ; i < 12 ; ++i) {
            edge e = create_edge(v, i);
            float v0 = values[e.first];
            float v1 = values[e.second];
            if (v0 != v1) {
                std::cerr << "E" << i << " (" << v0 << ", " << v1 << ")\n";
                vertex v(e, values, i);
                local_vertices.insert(v);
            }
        }
        if (!local_vertices.size()) {
            continue;
        }
        
        // loop over faces, create edges
        
        for (int i = 0 ; i < 6 ; ++i) {
            std::vector<vertex> selected, __vert(local_vertices.begin(), local_vertices.end());
            std::vector<t_edge> __edges;
            vertex center = vertex(create_face(v, i), values, i);
            face_query::filter(selected, __vert, v, i);
            bool included_center = process_face(__edges, selected, center);
            if (included_center) {
                face_vertices.insert(center);
            }
            
            std::sort(selected.begin(), selected.end());
            local_vertices.insert(selected.begin(), selected.end());
            std::sort(__edges.begin(), __edges.end(), Lt_t_edge());
            local_edges.insert(__edges.begin(), __edges.end());
        }
        if (face_vertices.size() == 2) {
            local_edges.insert(make_t_edge(face_vertices.front(), face_vertices.back()));
        } else if (face_vertices.size() >= 3) {
            // determine all boundary value pairs
            typedef std::pair<float, float>         value_pair;
            std::set<value_pair, Lt_pair<float> >   all_value_pairs;
            typedef std::set<value_pair, Lt_pair<float> >::const_iterator val_iterator;
            for (int i = 0 ; i < 12 ; ++i) {
                int j0 = voxel_edges[i][0];
                int j1 = voxel_edges[i][1];
                all_value_pairs.insert(create_pair(voxel_values[j0], voxel_values[j1]));
            }
            std::vector<vertex> __face_vertices(face_vertices.begin(), face_vertices.end());
            for (int i = 0 ; i < __face_vertices.size() - 1 ; ++i) {
                const vertex& v0 = __face_vertices[i];
                for (int j = i + 1 ; j < __face_vertices.size() ; ++j) {
                    const vertex& v1 = __face_vertices[j];
                    std
                    for (val_iterator it = all_value_pairs.begin() ;
                            it != all_value_pairs.end() ; ++it) {
                    }
                }
            }
        }
        
        edges.insert(local_edges.begin(), local_edges.end());
        
        std::cerr << "available t_edges for this voxel:\n";
        for (std::set<t_edge>::const_iterator it = local_edges.begin() ; it != local_edges.end() ; ++it) {
            std::cerr << *it << '\n';
        }
        
        // triangulate
        std::vector<t_edge> __edges(local_edges.begin(), local_edges.end());
        triangulate(triangles, __edges, voxel_values);
        vertices.insert(local_vertices.begin(), local_vertices.end());
    }
    
    std::cout << '\n';
    std::cout << vertices.size() << " vertices in triangulation\n";
    std::cout << triangles.size() << " triangles in triangulation\n";
    std::cout << "step = " << step << "\n";
    
    std::map<vertex, int>       pos_ids;
    std::vector<nvis::vec3>     pos;
    std::vector<int>            voxel_centers;
    std::vector<int>            face_centers;
    for (std::set<vertex>::const_iterator it = vertices.begin() ; it != vertices.end() ; ++it) {
        const vertex& v = *it;
        if (pos_ids.find(v) != pos_ids.end()) {
            // std::cout << "vertex " << v << " is already known (" << pos_ids[v] << ")\n";
            continue;
        }
        nvis::vec3 x = coord(v, size, step);
        pos.push_back(x);
        pos_ids[v] = pos.size() - 1;
        if (v.type == 2) {
            voxel_centers.push_back(pos.size() - 1);
        } else if (v.type == 1) {
            face_centers.push_back(pos.size() - 1);
        }
        
        if (pos.size() - 1 == 871) {
            std::cout << "vertex #" << pos.size() - 1 << " corresponds to " << to_str(v)
                      << " and has coordinates " << x << '\n';
        }
    }
    
    std::ostringstream os;
    std::string base(outs);
    size_t where = base.rfind('.');
    if (where != std::string::npos) {
        base.resize(where);
    }
    
    if (false) {
        std::vector<int> corner_voxels, edge_voxels, simple_voxels;
        classify(simple_voxels, edge_voxels, corner_voxels, values, size);
        
        std::cout << simple_voxels.size() << " simple voxels, "
                  << edge_voxels.size() << " edge voxels, "
                  << corner_voxels.size() << " corner voxels\n";
                  
        std::string name(base);
        name.append("-classified.nrrd");
        Nrrd* nout = nrrdNew();
        int* tags = (int*)calloc(nbvoxels, sizeof(int));
        for (int i = 0 ; i < simple_voxels.size() ; ++i) {
            tags[simple_voxels[i]] = 1;
        }
        for (int i = 0 ; i < edge_voxels.size() ; ++i) {
            tags[edge_voxels[i]] = 2;
        }
        for (int i = 0 ; i < corner_voxels.size() ; ++i) {
            tags[corner_voxels[i]] = 3;
        }
        size_t dims[3] = { size[0] - 1, size[1] - 1, size[2] - 1 };
        if (nrrdWrap_nva(nout, tags, nrrdTypeInt, 3, dims)) {
            std::cerr << biffGetDone(NRRD) << '\n';
        }
        if (nrrdSave(name.c_str(), nout, NULL)) {
            std::cerr << biffGetDone(NRRD) << '\n';
        }
    }
    
    std::vector<vertex> corner_vertex;
    std::vector<t_edge> real_edges, tmp_edges(edges.begin(), edges.end());
    topology(corner_vertex, real_edges, tmp_edges);
    
    // smooth the mesh while preserving corners
    if (niter > 0) {
        std::vector<tri_face> all_faces(triangles.size());
        std::vector<int> fixed(corner_vertex.size());
        for (int i = 0 ; i < triangles.size() ; ++i) {
            tri_face& f = all_faces[i];
            for (int j = 0 ; j < 3 ; ++j) {
                f[j] = pos_ids[triangles[i][j]];
            }
        }
        for (int i = 0 ; i < fixed.size() ; ++i) {
            fixed[i] = pos_ids[corner_vertex[i]];
        }
        smooth_mesh(pos, all_faces, fixed, alpha, niter);
    }
    
    {
        std::string name(base);
        name.append("-edges.vtk");
        std::fstream vtk(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << pos.size() << " float\n";
        for (int i = 0 ; i < pos.size() ; ++i) {
            vtk << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
        }
        vtk << "VERTICES 1 " << corner_vertex.size() + 1 << '\n' << corner_vertex.size();
        for (int i = 0 ; i < corner_vertex.size() ; ++i) {
            vtk << " " << pos_ids[corner_vertex[i]];
        }
        vtk << "\nLINES " << real_edges.size() << " " << real_edges.size()*3 << '\n';
        for (int i = 0 ; i < real_edges.size() ; ++i) {
            vtk << "2 " << pos_ids[real_edges[i].first] << " " << pos_ids[real_edges[i].second] << '\n';
        }
        vtk << '\n';
        vtk.close();
        
        name = base;
        name.append("-corners.vtk");
        vtk.open(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << corner_vertex.size() << " float\n";
        for (int i = 0 ; i < corner_vertex.size() ; ++i) {
            int j = pos_ids[corner_vertex[i]];
            vtk << pos[j][0] << " " << pos[j][1] << " " << pos[j][2] << '\n';
        }
        vtk.close();
        
        // check what we just exported
        if (false) {
            // initialize counter
            typedef std::set<nvis::vec3, Lt_close_enough>           pos_set;
            typedef std::map<nvis::vec3, pos_set, Lt_close_enough>  pos_map;
            pos_map index;
            for (int i = 0 ; i < corner_vertex.size() ; ++i) {
                nvis::vec3 x = pos[pos_ids[corner_vertex[i]]];
                index[x] = pos_set();
            }
            Lt_close_enough order;
            std::vector<double> dist;
            for (int i = 0 ; i < real_edges.size() ; ++i) {
                const t_edge& e = real_edges[i];
                pos_map::iterator it;
                nvis::vec3 x0 = pos[pos_ids[e.first]];
                nvis::vec3 x1 = pos[pos_ids[e.second]];
                it = index.find(x0);
                if (it != index.end()) {
                    it->second.insert(x1);
                }
                it = index.find(x1);
                if (it != index.end()) {
                    it->second.insert(x0);
                }
            }
            std::cout << "connectivity of exported corner vertices\n";
            for (pos_map::iterator it = index.begin() ;
                    it != index.end() ; ++it) {
                std::cout << "vertex at " << it->first << " has " << it->second.size() << " neighbors\n";
                
                if (it->second.size() == 2) {
                    std::cout << "WARNING\n\n\n\n";
                }
            }
        }
    }
    
    {
        std::string name(base);
        name.append("-all.vtk");
        std::fstream vtk(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << pos.size() << " float\n";
        for (int i = 0 ; i < pos.size() ; ++i) {
            vtk << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
        }
        vtk << "VERTICES 1 " << corner_vertex.size() + 1 << '\n' << corner_vertex.size();
        for (int i = 0 ; i < corner_vertex.size() ; ++i) {
            vtk << " " << pos_ids[corner_vertex[i]];
        }
        vtk << "\nLINES " << real_edges.size() << " " << real_edges.size()*3 << '\n';
        for (int i = 0 ; i < real_edges.size() ; ++i) {
            vtk << "2 " << pos_ids[real_edges[i].first] << " " << pos_ids[real_edges[i].second] << '\n';
        }
        vtk << "\nPOLYGONS " << triangles.size() << " " << 4*triangles.size() << '\n';
        for (int i = 0 ; i < triangles.size() ; ++i) {
            vtk << "3";
            for (int j = 0 ; j < 3 ; ++j) {
                vtk << " " << pos_ids[triangles[i][j]];
            }
            vtk << '\n';
        }
        vtk << '\n';
        vtk.close();
        
        name = base;
        name.append("-corners.vtk");
        vtk.open(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << corner_vertex.size() << " float\n";
        for (int i = 0 ; i < corner_vertex.size() ; ++i) {
            int j = pos_ids[corner_vertex[i]];
            vtk << pos[j][0] << " " << pos[j][1] << " " << pos[j][2] << '\n';
        }
        vtk.close();
    }
    
    {
        std::string name(base);
        name.append("-voxel-centers.vtk");
        std::fstream vtk(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << pos.size() << " float\n";
        for (int i = 0 ; i < pos.size() ; ++i) {
            vtk << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
        }
        vtk << "VERTICES 1 " << voxel_centers.size() + 1 << '\n' << voxel_centers.size();
        for (int i = 0 ; i < voxel_centers.size() ; ++i) {
            vtk << " " << voxel_centers[i];
        }
        vtk << '\n';
    }
    
    {
        std::string name(base);
        name.append("-face-centers.vtk");
        std::fstream vtk(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << pos.size() << " float\n";
        for (int i = 0 ; i < pos.size() ; ++i) {
            vtk << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
        }
        vtk << "VERTICES 1 " << face_centers.size() + 1 << '\n' << face_centers.size();
        for (int i = 0 ; i < face_centers.size() ; ++i) {
            vtk << " " << face_centers[i];
        }
        vtk << '\n';
    }
    
    {
        std::string name(base);
        name.append("-mesh.vtk");
        std::fstream vtk(name.c_str(), std::ios::out);
        vtk << "# vtk DataFile Version 2.0\n"
            << "Generated through multi-material marching cubes from file " << file << '\n'
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
    }
    
    if (!verbose) {
        std::cerr.rdbuf(init_buf);
    }
    
    return 0;
}









