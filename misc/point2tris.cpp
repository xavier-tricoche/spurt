#include <iostream>
#include <fstream>
#include <math/types.hpp>
#include <vector>
#include <list>
#include <teem/nrrd.h>
#include <misc/progress.hpp>
#include <limits>


std::vector<spurt::fvec3> points;
std::vector<spurt::ivec2> coords;
std::vector<bool> tagged;
std::vector<float> data;
std::vector<spurt::ivec3> tris;

typedef std::list<unsigned int>     list_type;
typedef list_type::const_iterator   const_list_iterator_type;
typedef list_type::iterator         list_iterator_type;

char* input, *output;
int size[2];
double bmin[3], bmax[3];
double h;
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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &input,       NULL,       "points text file");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &output,      NULL,       "output file name");
    hestOptAdd(&hopt, "r",      "resolution",       airTypeInt,     2, 2, size,         NULL,       "image raster resolution");
    hestOptAdd(&hopt, "min",    "min",              airTypeDouble,  3, 3, bmin,         NULL,       "bounding box minimum");
    hestOptAdd(&hopt, "max",    "max",              airTypeDouble,  3, 3, bmax,         NULL,       "bounding box maximum");
    hestOptAdd(&hopt, "h",      "step",             airTypeDouble,  1, 1, &h,           NULL,       "depth discontinuity threshold");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert image-based point cloud to triangulated surface",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}


struct raster {

    unsigned int idx(unsigned int u, unsigned int v) const {
        return u + v*__width;
    }
    
    raster(size_t width, size_t height)
        : __width(width), __height(height),
          __raster(width* height) {
    }
    
    void add(unsigned int u, unsigned int v, unsigned int i) {
        unsigned int id = idx(u, v);
        __raster[id].push_back(i);
    }
    
    bool empty(unsigned int u, unsigned int v) const {
        unsigned int id = idx(u, v);
        return __raster[id].empty();
    }
    
    const list_type& contents(unsigned int u, unsigned int v) const {
        return __raster[idx(u,v)];
    }
    
    size_t __width, __height;
    std::vector<list_type> __raster;
};

inline int get_neighbor(unsigned int u, unsigned int v, const raster& r,
                        const spurt::fvec3& x)
{
    if (u >= r.__width || v >= r.__height) {
        return -1;
    }
    
    const list_type& pixel = r.contents(u, v);
    const_list_iterator_type it;
    
    float mind = std::numeric_limits<float>::max();
    int mini = -1;
    
    for (it = pixel.begin() ; it != pixel.end() ; ++it) {
        if (tagged[*it]) {
            continue;
        }
        float d = spurt::norm(x - points[*it]);
        if (d < mind) {
            mind = d;
            mini = *it;
        }
    }
    
    return mini;
}

inline float diameter(const spurt::ivec3& tri)
{
    std::vector<float> d(3);
    for (int i = 0 ; i < 3 ; ++i) {
        d[i] = spurt::norm(points[tri[i]] - points[tri[(i+1)%3]]);
    }
    return *std::max_element(d.begin(), d.end());
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    std::fstream in(input, std::ios::in);
    
    size_t width = size[0];
    size_t height = size[1];
    raster image(width, height);
    
    spurt::fvec3 min(bmin[0], bmin[1], bmin[2]), max(bmax[0], bmax[1], bmax[2]);
    spurt::fvec3 diam(256., 256., 128.);
    
    spurt::timer timer;
    
    unsigned int npts;
    in >> npts;
    float x, y, z, u, v, lambda;
    for (int i = 0 ; i < npts ; ++i) {
        in >> u >> v >> x >> y >> z >> lambda;
        points.push_back(spurt::fvec3(x, y, z));
        coords.push_back(spurt::ivec2(u, v));
        data.push_back(lambda);
        image.add(u, v, points.size() - 1);
    }
    
    std::cerr << "input file read in " << timer.elapsed() << " seconds\n";
    
    tagged.resize(points.size());
    std::fill(tagged.begin(), tagged.end(), false);
    
    timer.start();
    for (unsigned int i = 0 ; i < points.size() ; ++i) {
        spurt::ivec2 uv = coords[i];
        int p0 = get_neighbor(uv[0] + 1, uv[1] - 1, image, points[i]);
        int p1 = get_neighbor(uv[0] + 1, uv[1], image, points[i]);
        int p2 = get_neighbor(uv[0], uv[1] + 1, image, points[i]);
        
        if (p0 >= 0 && p1 >= 0) {
            spurt::ivec3 tri(i, p0, p1);
            if (diameter(tri) < h) {
                tris.push_back(tri);
            }
        }
        if (p1 >= 0 && p2 >= 0) {
            spurt::ivec3 tri(i, p1, p2);
            if (diameter(tri) < h) {
                tris.push_back(tri);
            }
        }
    }
    std::cerr << points.size() << " processed in " << timer.elapsed() << " seconds\n";
    
    std::fstream out(output, std::ios::out);
    out << "# vtk DataFile Version 2.0\n"
        << "Converted from point cloud file " << argv[1] << '\n'
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << points.size() << " float\n";
    for (int i = 0 ; i < points.size() ; ++i) {
        spurt::fvec3 p = min + (points[i] / diam) * (max - min);
        out << p[0] << " " << p[1] << " " << p[2] << '\n';
    }
    out << "CELLS " << tris.size() << " " << 4*tris.size() << '\n';
    for (int i = 0 ; i < tris.size() ; ++i) {
        out << "3 " << tris[i][0] << " " << tris[i][1] << " " << tris[i][2] << '\n';
    }
    out << "CELL_TYPES " << tris.size() << '\n';
    for (int i = 0 ; i < tris.size() ; ++i) {
        out << "5\n";
    }
    out << "POINT_DATA " << points.size() << '\n';
    out << "SCALARS value float\n";
    out << "LOOKUP_TABLE default\n";
    for (int i = 0 ; i < data.size() ; ++i) {
        out << data[i] << '\n';
    }
    out.close();
    
    return 0;
}































