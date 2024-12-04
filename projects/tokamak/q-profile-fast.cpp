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

// spurt
#include <math/math.hpp>
#include <poincare/metric.hpp>

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
#include <maps/triangulation.hpp>
#include <maps/adaptive_triangulation.hpp>

#include <maps/period.hpp>
#include <data/quadtree.hpp>

int res, maxp, maxit, n1, n2, it, m, maxround, max_nb_triangles, max_period;
char* outs, *file, *ts;
double h, max_diff, min_area, max_area, max_ratio;

struct point_counter {
    typedef spurt::quadtree<int, nvis::vec2>   quadtree_type;
    
    point_counter() : n(0) {}
    
    void operator()(const quadtree_type& q) {
        n += q.data().size();
    }
    
    int n;
};

struct leaf_display {
    typedef spurt::quadtree<int, nvis::vec2>   quadtree_type;
    
    leaf_display() : os() {}
    void operator()(const quadtree_type& q) {
        os << q << '\n';
    }
    
    void reset() {
        os.str("");
        os.clear();
    }
    
    std::ostringstream os;
};

struct point_locator {

    typedef spurt::quadtree<int, nvis::vec2>   quadtree_type;
    typedef quadtree_type::data_type            data_type;
    
    point_locator(const nvis::bbox2& bounds, size_t max_depth, size_t max_nb_pts)
        : qt(bounds, max_depth, max_nb_pts) {}
        
    point_locator(const point_locator& pl) : qt(pl.qt) {}
    
    int find_close_point(const nvis::vec2& x, int id) {
        quadtree_type& leaf = qt.find(x);
        int res = (leaf.data().size() ? leaf.data()[0].second : 0);
        leaf.insert(x, id);
        return res;
    }
    
    void display() const {
        std::string str("");
        display(str);
    }
    
    void display(const std::string& preamble) const {
        leaf_display display;
        qt.for_each_leaf(display);
        std::cerr << preamble << "quadtree contains: " << display.os.str() << std::endl;
    }
    
    quadtree_type qt;
};

template<typename IntType>
double value(const boost::rational<IntType>& r)
{
    return (double)r.numerator() / (double)r.denominator();
}

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
    hestOptAdd(&hopt, "f",      "file",                 airTypeString,  1, 1, &file,                NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",                 airTypeString,  1, 1, &ts,                  NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &outs,                NULL,       "output nrrd file name");
    hestOptAdd(&hopt, "eps",    "eps",                  airTypeDouble,  0, 1, &h,                   "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "n1",     "# 1D samples",         airTypeInt,     0, 1, &n1,                  "100",      "number of 1D samples");
    hestOptAdd(&hopt, "n2",     "# 2D samples",         airTypeInt,     0, 1, &n2,                  "2000",     "number of 2D samples");
    hestOptAdd(&hopt, "m",      "# rotations",          airTypeInt,     0, 1, &m,                   "50",       "number of rotations");
    hestOptAdd(&hopt, "d",      "safety factor step",   airTypeDouble,  0, 1, &max_diff,            "0.01",     "maximum discrepancy between consecutive safety factors");
    hestOptAdd(&hopt, "mr",     "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,           "3.",       "maximum triangle aspect ratio");
    // hestOptAdd(&hopt, "nr",  "max # rounds",         airTypeInt,     0, 1, &maxround,            "10",       "max number of 1D approximation iterations");
    hestOptAdd(&hopt, "mt",     "max # triangles",      airTypeInt,     0, 1, &max_nb_triangles,    "500000",   "max number of triangles in adaptive sampling");
    hestOptAdd(&hopt, "ma",     "min triangle area",    airTypeDouble,  1, 1, &min_area,            NULL,       "min triangle area in adaptive sampling");
    hestOptAdd(&hopt, "Ma",     "max triangle area",    airTypeDouble,  1, 1, &max_area,            NULL,       "max triangle area in adaptive sampling");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0, 1, &max_period,          "25",       "max considered period in adaptive refinement");
    
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Compute approximative q-profile using magnetic field integration",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct Chain {
    std::vector< nvis::vec2 > x;
    double q;
    double varq;
};

typedef std::pair<double, double>                                       pair_type;
typedef std::pair<nvis::vec2, double>                                   value_type;
typedef std::vector<nvis::vec2>::const_iterator                         const_iter_type;
typedef spurt::triangulation<double, point_locator>                    mesh_type;
typedef mesh_type::index_type                                           index_type;
typedef mesh_type::triangle_type                                        triangle_type;
typedef mesh_type::point_type                                           point_type;

std::vector<boost::rational<int> > ref_values;

struct triangle_priority {
    index_type id;
    double span;
    double area;
    int version;
    double ratio;
    bool interesting;
};

struct Gt_triangle_priority {
    int operator()(const triangle_priority& p0,
                   const triangle_priority& p1) const {
                   
        // priority rules:
        // 1) area must be less than max area
        if (p0.area > p1.area && p0.area > max_area) {
            return true;
        } else if (p0.area <= min_area && p1.area > max_area) {
            return false;
        }
        
        // 2) value span must be smaller than max span
        if (p0.span > p1.span && p0.span > max_diff) {
            return true;
        } else if (p0.span < p1.span && p1.span > max_diff) {
            return false;
        }
        
        // 3) triangle must be "interesting"
        if (p0.interesting && !p1.interesting) {
            return true;
        } else if (!p0.interesting && p1.interesting) {
            return false;
        }
        
        // 4) aspect ratio must be under max aspect ratio
        if (p0.ratio > p1.ratio && p0.ratio > max_ratio) {
            return true;
        } else if (p0.ratio < p1.ratio && p1.ratio > max_ratio) {
            return true;
        }
        
        // if none of the above criteria can be used to sort the triangles, use
        // their value span
        return (p0.span > p1.span);
    }
};

inline double area(const nvis::vec2 p[3])
{
    nvis::vec2 e0 = p[1] - p[0];
    nvis::vec2 e1 = p[2] - p[0];
    double det = e0[0] * e1[1] - e0[1] * e1[0];
    return 0.5*fabs(det);
}

inline double span(const double v[3])
{
    return *std::max_element(&v[0], &v[3]) - *std::min_element(&v[0], &v[3]);
}

inline double aspect_ratio(const nvis::vec2 p[3])
{
    double l[3];
    for (int i = 0 ; i < 3 ; ++i) {
        l[i] = nvis::norm(p[i] - p[(i+1)%3]);
    }
    return *std::max_element(&l[0], &l[3]) / *std::min_element(&l[0], &l[3]);
}

template<typename T>
struct interval {
    interval(const T& min, const T& max) : __min(min), __max(max) {
        if (__max < __min) {
            T tmp = __min;
            __min = __max;
            __max = tmp;
        }
    }
    
    bool inside(const T& val) const {
        return (!(val < __min) && !(__max < val));
    }
    
    T __min, __max;
};
typedef interval<double>    interval_type;

bool interesting(const double v[3])
{
    interval_type span(*std::min_element(&v[0], &v[3]), *std::max_element(&v[0], &v[3]));
    for (int i = 0 ; i < ref_values.size() ; ++i) {
        double v = value(ref_values[i]);
        if (span.inside(v)) {
            return true;
        }
    }
    return false;
}

triangle_priority priority(const mesh_type& mesh, index_type id, int version = 0)
{
    nvis::vec2 p[3];
    double v[3];
    const triangle_type& triangle = mesh.get_triangle_vertices(id);
    for (int n = 0 ; n < 3 ; ++n) {
        p[n] = mesh.get_vertex(triangle[n]);
        v[n] = mesh.get_data(triangle[n]);
    }
    
    triangle_priority q;
    q.id = id;
    q.area = area(p);
    q.span = span(v);
    q.version = version;
    q.ratio = aspect_ratio(p);
    q.interesting = interesting(v);
    return q;
}

double period(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
              const poincare_map& __pmap, const spurt::map_metric& metric)
{
    const poincare_map* pmap = __pmap.clone();
    try {
        pmap->map(x0, steps, m);
    } catch (...) {
        return -1;
    }
    
    pair_type sf = spurt::period_x_periodic(steps, x0, metric);
    
    // double check = spurt::period_x_periodic_fourier(steps, x0, metric).first;
    // std::cerr << "standard method found " << sf.first << " while fourier found " << check
    // << std::endl;
    
    return sf.first;
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(file), std::string(ts));
    field->periodic_coordinates(false);
    
    bool per[2] = {true, false};
    spurt::map_metric metric(field->bounds(), per);
    
    poincare_map pmap(field);
    pmap.precision(h);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double width = diagonal[0];
    double height = diagonal[1];
    double dh = height / (double)(n1 - 1);
    double x_center = bounds.center()[0];
    
    typedef boost::rational<int>    rational_type;
    std::set<rational_type> rationals;
    for (int num = 1 ; num <= max_period ; ++num) {
        for (int den = 1 ; den <= num ; ++den) {
            rationals.insert(rational_type(num, den));
        }
    }
    ref_values.clear();
    std::copy(rationals.begin(), rationals.end(),
              std::back_inserter(ref_values));
    std::cerr << "there are " << ref_values.size() << " interesting values\n";
    for (int i = 0 ; i < ref_values.size() ; ++i) {
        std::cerr << ref_values[i] << " = " << value(ref_values[i]) << '\n';
    }
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    std::vector< Chain > chains;
    
    srand48(time(0));
    std::vector< value_type > all_points;
    
#pragma openmp parallel
    {
        const poincare_map* clone = pmap.clone();
        std::vector< nvis::vec2 > steps;
        
        #pragma omp for schedule(dynamic,1)
        for (int i = 0 ; i < n1 ; ++i) {
#if _OPENMP
            const int thread = omp_get_thread_num();
            if (!thread) {
                std::cout << i << " from " << n1 << " values computed\n";
            }
#endif
            nvis::vec2 x0 = bounds.min() + nvis::vec2(0.5 * diagonal[0], i * dh);
            double q = period(x0, steps, pmap, metric);
            if (q <= 0) {
                continue;
            }
            
            for (const_iter_type it = steps.begin() ; it != steps.end() ; ++it) {
                all_points.push_back(std::make_pair(metric.modulo(*it), q));
            }
        }
    }
    std::cout << all_points.size() << " points sampled after 1D seeding\n";
    
#pragma openmp parallel
    {
        const poincare_map* clone = pmap.clone();
        std::vector< nvis::vec2 > steps;
        
        #pragma omp for schedule(dynamic,1)
        for (int i = 0 ; i < n2 ; ++i) {
#if _OPENMP
            const int thread = omp_get_thread_num();
            if (!thread) {
                std::cout << i << " from " << n2 << " values computed\n";
            }
#endif
            nvis::vec2 x0 = bounds.min() + nvis::vec2(drand48() * width, drand48() * height);
            double q = period(x0, steps, pmap, metric);
            if (q <= 0) {
                continue;
            }
            
            for (const_iter_type it = steps.begin() ; it != steps.end() ; ++it) {
                all_points.push_back(std::make_pair(metric.modulo(*it), q));
            }
        }
    }
    std::cout << all_points.size() << " points sampled after random 2D seeding\n";
    
    // 2D quality control
    
    // create initial Delaunay triangulation first
    std::cerr << "building initial triangulation\n";
    mesh_type mesh(bounds, -1, point_locator(bounds, 10, 10));
    mesh.insert_points(all_points);
    
    std::cerr << "after initial triangulation, mesh contains " << mesh.get_nb_triangles()
              << " triangles\n";
              
    int nb_triangles = mesh.get_nb_triangles();
    std::set<triangle_priority, Gt_triangle_priority> priority_list;
    for (int i = 0 ; i < nb_triangles ; ++i) {
        priority_list.insert(priority(mesh, i));
    }
    
    int last_pct = 100 * nb_triangles / max_nb_triangles;
    
    std::map<index_type, int> current_version;
    
    while (nb_triangles < max_nb_triangles && priority_list.begin() != priority_list.end()) {
        triangle_priority q = *priority_list.begin();
        std::map<int, int>::iterator iter = current_version.find(q.id);
        if (iter != current_version.end() && q.version < iter->second) {
            // we have already replaced that entry by a more recent one. remove it
            priority_list.erase(priority_list.begin());
            continue;
        } else if (q.area < min_area) {
            priority_list.erase(priority_list.begin());
            continue;
        } else if (q.span < max_diff && q.area < max_area && !q.interesting) {
            priority_list.erase(priority_list.begin());
            continue;
        } else {
            nvis::vec2 p[3];
            const triangle_type& triangle = mesh.get_triangle_vertices(q.id);
            for (int n = 0 ; n < 3 ; ++n) {
                p[n] = mesh.get_vertex(triangle[n]);
            }
            nvis::vec2 seed = 1. / 3.*(p[0] + p[1] + p[2]);
            const poincare_map* clone = pmap.clone();
            std::vector< nvis::vec2 > steps;
            try {
                clone->map(seed, steps, m);
            } catch (...) {
                priority_list.erase(priority_list.begin());
                continue;
            }
            
            pair_type sf = spurt::period_x_periodic(steps, seed, metric);
            
            double val = sf.first;
            if (val <= 0) {
                priority_list.erase(priority_list.begin());
                continue;
            }
            
            for (int i = 0 ; i < steps.size() ; ++i) {
                steps[i] = metric.modulo(steps[i]);
            }
            
            steps.push_back(seed);
            std::vector<double> data(steps.size());
            std::fill(data.begin(), data.end(), val);
            std::vector<int> modified;
            mesh.insert_points(steps, data, modified);
            for (int n = 0 ; n < modified.size() ; ++n) {
                int id = modified[n];
                int version = 0;
                if (id < nb_triangles) {
                    iter = current_version.find(id);
                    if (iter != current_version.end()) {
                        version = iter->second + 1;
                    } else {
                        version = 1;
                    }
                    current_version[id] = version;
                }
                priority_list.insert(priority(mesh, id, version));
            }
        }
        
        nb_triangles = mesh.get_nb_triangles();
        int pct = 100 * nb_triangles / max_nb_triangles;
        if (pct > last_pct) {
            std::cerr << nb_triangles << " triangles for " << mesh.get_nb_vertices() << " vertices (" <<
                      pct << "%)  \r";
        }
        last_pct = pct;
    }
    std::cerr << '\n';
    
    std::ostringstream os;
    os << outs << "-n1=" << n1 << "-n2=" << n2 << "-m=" << m << ".vtk";
    std::fstream vtk(os.str().c_str(), std::ios::out);
    vtk << "# vtk DataFile Version 2.0\n"
        << "Sample points from " << file << " at t=" << ts << " with estimated safety factor\n"
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";
        
    int nb_pts = mesh.get_nb_vertices();
    vtk << "POINTS " << nb_pts << " float\n";
    for (int i = 0 ; i < nb_pts ; ++i) {
        const nvis::vec2& p = mesh.get_vertex(i);
        vtk << p[0] << " " << p[1] << " 0.\n";
    }
    
    int ncells = mesh.get_nb_triangles();
    vtk << "CELLS " << ncells << " " << 4*ncells << '\n';
    for (int i = 0 ; i < ncells ; ++i) {
        const triangle_type& tri = mesh.get_triangle_vertices(i);
        vtk << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }
    
    vtk << "CELL_TYPES " << ncells << '\n';
    for (int i = 0 ; i < ncells ; ++i) {
        vtk << "5\n";
    }
    
    vtk << "POINT_DATA " << nb_pts << '\n'
        << "SCALARS q float\n"
        << "LOOKUP_TABLE default\n";
    for (int i = 0 ; i < nb_pts ; ++i) {
        double v = mesh.get_data(i);
        vtk << v << '\n';
    }
    vtk.close();
}







































































































































































