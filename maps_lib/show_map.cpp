#include <iostream>
#include <list>
#include <map>
#include <iomanip>

#include <math/fixed_vector.hpp>
#include <vector>
#include "display.hpp"
#include "definitions.hpp"
#include "xmt_poincare_map.hpp"
#include <data/grid.hpp>
#include <data/raster.hpp>
#include <maps/period.hpp>
#include "map_field_wrapper.hpp"

#include <util/wall_timer.hpp>

#include <graphics/GLUT_helper.hpp>

#include <image/nrrd_wrapper.hpp>

#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>

#if _OPENMP
#include <omp.h>
#endif

using namespace xavier;
using namespace map_analysis;
using namespace map_display;
// using namespace div_cleaning;

namespace {
double __mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}
}

std::vector<std::pair<nvis::vec2, nvis::fvec3> > seeds;

typedef xavier::grid<double, 3>                 grid_type;
typedef raster_data<nvis::vec3, 3, double>      trilinear_field_type;
typedef divfree_field<trilinear_field_type>     divfree_field_type;
typedef nvis::fixed_vector<int, 3>              ivec_type;

xavier::map_metric  metric2d;

typedef xmt_poincare_map<xavier::map::wrapper<trilinear_field_type> >   trilinear_map_type;
typedef xmt_poincare_map<xavier::map::wrapper<divfree_field_type> >     divfree_map_type;

nvis::bbox2 _bounds;

char*    in, *seed_l, *seed_s, *out;
int     npts, niter, rx, dolink;
double  eps;
bool    div_free;

nvis::vec3 hsv2rgb(double hue, double sat, double val)
{
    double chroma = val * sat;
    double h_ = hue / 60.;
    double x = chroma * (1. - fabs(fmod(h_, 2.) - 1));
    
    nvis::vec3 rgb_;
    if(val == 0) {
        rgb_ = nvis::vec3(0, 0, 0);
    } else if(0 <= h_ && h_ < 1) {
        rgb_ = nvis::vec3(chroma, x, 0);
    } else if(1 <= h_ && h_ < 2) {
        rgb_ = nvis::vec3(x, chroma, 0);
    } else if(2 <= h_ && h_ < 3) {
        rgb_ = nvis::vec3(0, chroma, x);
    } else if(3 <= h_ && h_ < 4) {
        rgb_ = nvis::vec3(0, x, chroma);
    } else if(4 <= h_ && h_ < 5) {
        rgb_ = nvis::vec3(x, 0, chroma);
    } else if(5 <= h_ && h_ < 6) {
        rgb_ = nvis::vec3(chroma, 0, x);
    }
    
    double m = val - chroma;
    return rgb_ + nvis::vec3(m, m, m);
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
    hestOptAdd(&hopt, "i",  "input",                airTypeString,  1,  1,  &in,        NULL,           "input file name (NRRD)");
    hestOptAdd(&hopt, "o",  "output",               airTypeString,  0,  1,  &out,       "none",         "output file name (NRRD)");
    hestOptAdd(&hopt, "sl", "seed load",            airTypeString,  0,  1,  &seed_l,    "none",         "file containing seed locations (txt)");
    hestOptAdd(&hopt, "ss", "seed save",            airTypeString,  0,  1,  &seed_s,    "none",         "file containing seed locations (txt)");
    hestOptAdd(&hopt, "n",  "# pts",                airTypeInt,     0,  1,  &npts,      "500",          "number of seeds");
    hestOptAdd(&hopt, "m",  "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "link", NULL,                 airTypeInt,     0,  0,  &dolink,    NULL,           "connect iterates to form curves");
    hestOptAdd(&hopt, "e",  "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "df", "div free",             airTypeBool,    0,  1,  &div_free,  "0",            "divergence free interpolation?");
    hestOptAdd(&hopt, "rx", "resolution factor",    airTypeInt,     0,  1,  &rx,        "10",           "output resolution upsample");
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Poincare map visualization",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef std::list<nvis::vec2>                           orbit_type;
typedef std::list<std::pair<nvis::vec2, nvis::vec2> >   orbit_edges_type;
typedef std::pair<orbit_type, nvis::fvec3>              colored_orbit_type;
typedef std::pair<orbit_edges_type, nvis::fvec3>        colored_orbit_edges_type;

std::list<colored_orbit_type>               pplot;      // poincare plot
std::list<colored_orbit_edges_type>         pplot_links;

template<typename T1, typename T2>
void compare_fields(const grid_type& domain, const T1& field1, const T2& field2)
{
    nvis::bbox3 bounds = domain.bounds();
    srand48(time(0));
    
    try {
        for(int i = 0 ; i < 100 ; ++i) {
            nvis::vec3 u(drand48(), drand48(), drand48());
            nvis::vec3 x = bounds.min() + u * bounds.size();
            
            nvis::vec3 f1 = field1(x);
            nvis::vec3 f2 = field2(x);
            if(nvis::norm(f1 - f2) / nvis::norm(f1) > 0.1) {
                std::cerr << "WARNING: f1 = " << f1 << " != f2 = " << f2 << '\n';
                field1.verbose(true);
                field2.verbose(true);
                f1 = field1(x);
                f2 = field2(x);
                field1.verbose(false);
                field2.verbose(false);
            }
        }
    } catch(...) {
        field1.verbose(false);
        field2.verbose(false);
    }
}

// --------------------------------------------------------------------------

inline double eval(const rational_type& q)
{
    return (double)q.numerator() / (double)q.denominator();
}

inline double mindist(double q, int n)
{
    double _min = std::numeric_limits<double>::max();
    for(int i = 1 ; i <= n ; ++i) {
        rational_type r(n, i);
        double d = fabs(eval(r) - q);
        if(d < _min) {
            _min = d;
        }
    }
    
    return _min;
}

template<typename Field>
double cost(const std::vector<nvis::vec2>& steps, unsigned int period,
            const Field& field)
{
    if(period >= steps.size()) {
        return std::numeric_limits<double>::max();
    }
    double c = 0;
    for(int i = 0 ; i < steps.size() - period ; ++i) {
        // std::cerr << "distance between " << steps[i] << " and " << steps[i+period]
        //  << " is " << field.distance(steps[i], steps[i+period]) << std::endl;
        c += field.distance(steps[i], steps[i+period]);
    }
    return c / (double)(steps.size() - period);
}

inline double __angle(const nvis::vec2& x, const nvis::vec2& y)
{
    double cosine = nvis::inner(x, y) / (nvis::norm(x) * nvis::norm(y));
    return acos(cosine) / M_PI * 180.;
}

template<typename M>
void connect(const std::vector<nvis::vec2>& steps, const M& pmap)
{
    typedef std::pair<int, int>     edge_type;
    std::map<int, int>              edges;
    std::multimap<int, int>         p2p;
    
    const bool periodic[] = {true, false};
    std::ostringstream os;
    
    std::map<double, int> scores;
    std::map<int, double> id2score;
    for(unsigned int i = 1 ; i < steps.size() / 3 ; ++i) {
        double c = cost(steps, i, pmap.field());
        scores.insert(std::pair<double, int>(c, i));
        id2score.insert(std::pair<int, double>(i, c));
        if(c == 0) {
            os << "WARNING: zero cost for " << steps.size() << " steps" << std::endl;
            std::cerr << os.str();
        }
    }
    
    if(scores.empty()) {
        return;
    }
    
    int bestp = scores.begin()->second;
    
    os.str("");
    os.clear();
    os << "best period found is " << bestp << ", for average distance of " << std::setprecision(20) << scores.begin()->first << std::endl;
    std::cerr << os.str();
    
    for(int i = 0 ; i < steps.size() - bestp ; ++i) {
        edges.insert(std::pair<int, int>(i, i + bestp));
    }
    
    pplot_links.push_back(colored_orbit_edges_type());
    colored_orbit_edges_type& coe = pplot_links.back();
    orbit_edges_type& orb = coe.first;
    for(std::map<int, int>::const_iterator it = edges.begin() ; it != edges.end() ; ++it) {
        nvis::vec2 x = pmap.field().modulo(steps[it->first]);
        nvis::vec2 y = pmap.field().modulo(steps[it->second]);
        
#if 0
        std::list<std::pair<nvis::vec2, nvis::vec2> > tmp;
        xavier::clip(tmp, x, y, metric2d);
        for(std::list<std::pair<nvis::vec2, nvis::vec2> >::const_iterator jit = tmp.begin();
                jit != tmp.end() ; ++jit) {
            orb.push_back(*jit);
        }
#else
        std::vector<std::pair<nvis::vec2, nvis::vec2> > segments;
        metric2d.clip_segment(segments, x, y);
        for(int j = 0 ; j < segments.size() ; ++j) {
            orb.push_back(segments[j]);
        }       //
        // nvis::vec2 dx = metric2d.displacement(x, y, periodic);
        // if (dx[0] > 50) {
        //  os.str(""); os.clear();
        //  os << "ERROR: displacement between " << x << " and " << y
        //      << " is " << dx << std::endl;
        //  std::cerr << os.str();
        // }
        // if (metric2d.bounds().inside(x+dx))
        //  orb.push_back(std::pair<nvis::vec2, nvis::vec2>(x, y));
        // else {
        //  orb.push_back(std::pair<nvis::vec2, nvis::vec2>(x, x+dx));
        //  orb.push_back(std::pair<nvis::vec2, nvis::vec2>(y, y-dx));
        // }
#endif
    }
    
#if 1
    
    for(std::map<int, int>::const_iterator it = edges.begin() ; it != edges.end() ; ++it) {
        int a = it->first;
        int b = it->second;
        p2p.insert(edge_type(a, b));
        p2p.insert(edge_type(b, a));
    }
    
    // identify loose ends
    std::map<int, int> loose;
    for(int i = 0 ; i < steps.size() ; ++i) {
        if(p2p.count(i) == 1) {
            loose[i] = edges[i];
        }
    }
    os.clear();
    os.str("");
    os << "there are " << loose.size() << " loose ends in current curve" << std::endl;
    std::cerr << os.str();
    
    std::set<int> available;
    for(std::map<int, int>::const_iterator it = loose.begin() ; it != loose.end() ; ++it) {
        available.insert(it->first);
    }
    
    for(std::set<int>::iterator it = available.begin() ; it != available.end() ; ++it) {
        std::map<double, int>::const_iterator jt = scores.begin();
        for(++jt; jt != scores.end() ; ++jt) {
            if(available.find(jt->second) != available.end()) { //&& jt->first < 0.005*metric2d.diameter()) {
                available.erase(jt->second);
                std::vector<std::pair<nvis::vec2, nvis::vec2> > segments;
                nvis::vec2 x = pmap.field().modulo(steps[*it]);
                nvis::vec2 y = pmap.field().modulo(steps[jt->second]);
                metric2d.clip_segment(segments, x, y);
                for(int j = 0 ; j < segments.size() ; ++j) {
                    orb.push_back(segments[j]);
                }
            }
        }
    }
#else
    for(std::set<int>::iterator it = available.begin() ; it != available.end() ; ++it) {
        std::set<int>::iterator jt = it, best_it;
        double mind = std::numeric_limits<double>::max();
        nvis::vec2 prev = metric2d.displacement(steps[edges[*it]], steps[*it]);
        for(++jt ; jt != available.end() ; ++jt) {
            double d = pmap.field().distance(steps[*jt], steps[*it]);
            if(d >= mind || d >= 0.05 * metric2d.diameter()) {
                continue;
            }
            nvis::vec2 cur = metric2d.displacement(steps[*it], steps[*jt]);
            nvis::vec2 next = metric2d.displacement(steps[*jt], steps[edges[*jt]]);
            if(__angle(prev, cur) < 10. && __angle(cur, next) < 10.) {
                mind = d;
                best_it = jt;
            }
        }
        if(mind < std::numeric_limits<double>::max()) {
            std::vector<std::pair<nvis::vec2, nvis::vec2> > segments;
            nvis::vec2 x = pmap.field().modulo(steps[*it]);
            nvis::vec2 y = pmap.field().modulo(steps[*best_it]);
            metric2d.clip_segment(segments, x, y);
            for(int j = 0 ; j < segments.size() ; ++j) {
                orb.push_back(segments[j]);
            }
            available.erase(best_it);
        }
    }
#endif
    
    coe.second = hsv2rgb(drand48() * 360., 1, 0.5);
}


void display();

template<typename M>
void poincare_plot(std::list<colored_orbit_type>& points, const M& pmap,
                   size_t nseeds, size_t niter)
{
    typedef M       map_type;
    
    points.clear();
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<std::list<colored_orbit_type> > orbits(nbthreads);
    
    srand48(time(0));
    
    std::cerr << "bounds are " << _bounds << '\n';
    
    if(strcmp(seed_l, "none")) {
        // load seeds from file
        std::fstream input(seed_l, std::ios::in);
        input >> nseeds;
        std::cerr << nseeds << " seeds in input\n";
        seeds.resize(nseeds);
        float x, y, r, g, b;
        for(int i = 0 ; i < nseeds ; ++i) {
            input >> x >> y >> r >> g >> b;
            seeds[i].first = nvis::vec2(x, y);
            seeds[i].second = nvis::fvec3(r, g, b);
        }
        input.close();
        std::cerr << "seeds loaded\n";
    } else {
        seeds.resize(nseeds);
        for(int i = 0 ; i < nseeds ; ++i) {
            nvis::vec2 loc(drand48(), drand48());
            seeds[i].first = _bounds.min() + loc * _bounds.size();
            seeds[i].second = nvis::fvec3(drand48(), drand48(), drand48());
        }
        if(strcmp(seed_s, "none")) {
            std::cerr << "saving seeds to file\n";
            std::fstream output(seed_s, std::ios::out);
            output << nseeds << '\n';
            for(int i = 0 ; i < nseeds ; ++i) {
                const nvis::vec2& xy = seeds[i].first;
                const nvis::fvec3& rgb = seeds[i].second;
                output << xy[0] << " " << xy[1] << " " << rgb[0] << " " << rgb[1] << " " << rgb[2] << '\n';
            }
            output.clear();
            std::cerr << nseeds << " seeds exported\n";
        }
    }
    
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int i = 0 ; i < nseeds ; ++i) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            if(!thread_id) {
                std::cerr << i << "/" << nseeds << "          \r" << std::flush;
            }
            
            nvis::vec2 x0 = seeds[i].first;
            
            std::vector<nvis::vec2> steps;
            map_type* my_map = pmap.clone();
            try {
                my_map->map(x0, steps, niter);
            } catch(...) {
                // continue;
            }
            
            if(dolink) {
                connect(steps, pmap);
            }
            
            orbits[thread_id].push_back(colored_orbit_type());
            colored_orbit_type& orb = orbits[thread_id].back();
            for(unsigned int n = 0 ; n < steps.size() ; ++n) {
                orb.first.push_back(pmap.field().modulo(steps[n]));
            }
            orb.second = seeds[i].second;
        }
    }
    for(int i = 0 ; i < orbits.size() ; ++i) {
        std::copy(orbits[i].begin(), orbits[i].end(), std::back_inserter(points));
        std::cerr << orbits[i].size() << " orbits computed by thread #" << i << '\n';
    }
    
    std::cerr << '\n';
    std::cerr << "integration done\n" << points.size() << " orbits computed\n";
}

// --------------------------------------------------------------------------------

static void init()
{
    Nrrd* nin = nrrdNew();
    nin = xavier::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> __array;
    xavier::to_vector(__array, nin);
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    std::vector<nvis::vec3> vectors(__array.size() / 3);
    for(int i = 0 ; i < __array.size() / 3 ; ++i) {
        vectors[i][0] = __array[3*i  ];
        vectors[i][1] = __array[3*i+1];
        vectors[i][2] = __array[3*i+2];
    }
    trilinear_field_type basic_field(domain, vectors);
    
    // // begin debug
    // {
    //  int int_count = 0;
    //  int float_count = 0;
    //  int local_count = 0;
    //  srand48(time(0));
    //
    //  for (int i = 0 ; i < 2000 ; ++i) {
    //      // random cell coordinates
    //      int i = floor((dims[0]-1) * drand48());
    //      int j = floor((dims[1]-1) * drand48());
    //      int k = floor((dims[2]-1) * drand48());
    //      ivec_type c(i, j, k);
    //
    //      // random offset
    //      int N = floor(-10 + 20 * drand48());
    //      int M = floor(-10 + 20 * drand48());
    //      ivec_type d(i, j + N*(dims[1]-1), k + M*(dims[2]-1));
    //      try {
    //          if (nvis::any(c != domain.imodulo(d))) {
    //              std::cerr << c << " is different from modulo(" << d << ") with dims = " << dims << " which was found to be "
    //                        << domain.imodulo(d) << '\n';
    //              exit(-1);
    //          };
    //      }
    //      catch (std::runtime_error& e) {
    //          std::cerr << "exception caught: " << e.what() << " with d = " << d << '\n';
    //          exit(-1);
    //      }
    //
    //      // random local coordinates
    //      double x = drand48();
    //      double y = drand48();
    //      double z = drand48();
    //      vec_type loc(x, y, z);
    //      vec_type p = domain.bounds().min() + (vec_type(d) + loc) * domain.spacing();
    //      vec_type q = domain.bounds().min() + (vec_type(c) + loc) * domain.spacing();
    //
    //      std::pair<ivec_type, vec_type> tmp0 = domain.local_coordinates(p);
    //      std::pair<ivec_type, vec_type> tmp1 = domain.local_coordinates(q);
    //      if (nvis::any(tmp0.first != c) || nvis::any(tmp0.second != loc) || nvis::any(tmp1.first != c) || nvis::any(tmp1.second != loc)) ++local_count;
    //      if (nvis::any(q != domain.dmodulo(p))) ++float_count;
    //  }
    //
    //  std::cerr << int_count << "\% int modulo were wrong\n"
    //            << float_count << "\% float modulo were wrong\n"
    //            << local_count << "\% local coordinates were wrong\n";
    // }
    // // end debug
    
    
    double h = eps;
    
    trilinear_map_type basic_map(basic_field);
    basic_map.precision(eps);
    
    divfree_field<trilinear_field_type> divfree_field(domain, basic_field);
    divfree_field.verbose(false);
    
    // compare_fields(domain, basic_field, divfree_field);
    
    divfree_map_type divfree_map(divfree_field);
    divfree_map.precision(eps);
    
    grid_type::bounds_type bbox(domain.bounds());
    _bounds = nvis::bounding_box<nvis::vec2> (nvis::vec2(bbox.min()[1], bbox.min()[0]),
              nvis::vec2(bbox.max()[1], bbox.max()[0]));
              
    metric2d.bounds() = _bounds;
    metric2d.periodic(0) = true;
    metric2d.periodic(1) = false;
    
    std::cerr << "distance between (119, 1) and (1, 1) is " << metric2d.distance(nvis::vec2(119, 1), nvis::vec2(1, 1)) << std::endl;
    
    if(!div_free) {
        poincare_plot(pplot, basic_map, npts, niter);
    } else {
        poincare_plot(pplot, divfree_map, npts, niter);
    }
    
    // initialize OpenGL
    
    if(strcmp(out, "none")) {
        int nx = _bounds.size()[0] * rx;
        int ny = _bounds.size()[1] * rx;
        unsigned char* raster = (unsigned char*)calloc(3 * nx * ny, sizeof(unsigned char));
        std::list<colored_orbit_type>::iterator it;
        for(it = pplot.begin() ; it != pplot.end() ; ++it) {
            nvis::fvec3 c = it->second;
            c *= 256;
            unsigned char col[3] = { static_cast<unsigned char>(floor(c[0])), 
                                     static_cast<unsigned char>(floor(c[1])), 
                                     static_cast<unsigned char>(floor(c[2])) };
            std::list<nvis::vec2>::const_iterator it2;
            for(it2 = it->first.begin() ; it2 != it->first.end() ; ++it2) {
                nvis::vec2 x =  rx * *it2;
                int id = floor(x[0]) + nx * floor(x[1]);
                
                raster[3*id] = col[0];
                raster[3*id+1] = col[1];
                raster[3*id+2] = col[2];
            }
        }
        
        Nrrd* nout = nrrdNew();
        size_t __size[] = {3, static_cast<size_t>(nx), static_cast<size_t>(ny)};
        nrrdWrap_nva(nout, raster, nrrdTypeUChar, 3, __size);
        nrrdSave(out, nout, NULL);
        std::cerr << "poincare plot exported\n";
    }
    
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    
    glShadeModel(GL_FLAT);
    glClearColor(1.0, 1.0, 1.0, 0.0);
}

// --------------------------------------------------------------------------------

void display(void)
{
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    
#if 0
    // draw mesh
    glEnable(GL_BLEND);
    glLineWidth(1.0);
    glColor3f(.5, .5, .5);
    glEnable(GL_LINE_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glVertexPointer(2, GL_DOUBLE, 0, &grid.front());
    glDrawArrays(GL_QUADS, 0, grid.size());
    
    // draw cursor
    glEnable(GL_BLEND);
    glPointSize(10.0);
    glColor3f(1.0, 0.7, 0);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    glVertex2dv((GLdouble*)&cursor);
    glEnd();
    
    // draw seed candidates
    glEnable(GL_BLEND);
    glPointSize(3.0);
    glColor3f(0.0, 0.0, 0.0);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    for(std::vector<seedpoint>::iterator si = spts.begin(); si != spts.end(); ++si) {
        glVertex2dv((GLdouble*)&si->pos);
    }
    glEnd();
#endif
    
    if(dolink) {
        glEnable(GL_BLEND);
        // glEnable(GL_LINE_SMOOTH);
        glLineWidth(1.0);
        for(std::list<colored_orbit_edges_type>::const_iterator it = pplot_links.begin() ;
                it != pplot_links.end() ; ++it) {
            const orbit_edges_type& edges = it->first;
            glColor3f(it->second[0], it->second[1], it->second[2]);
            glBegin(GL_LINES);
            for(orbit_edges_type::const_iterator jt = edges.begin() ; jt != edges.end() ; ++jt) {
                const nvis::vec2& x = jt->first;
                const nvis::vec2& y = jt->second;
                glVertex2f(x[0], x[1]);
                glVertex2f(y[0], y[1]);
            }
            glEnd();
        }
    }
    // else
    {
        // draw poincare plot - each orbit in a different color
        std::list<colored_orbit_type>::iterator it;
        unsigned int count = 0;
        glDisable(GL_POINT_SMOOTH);
        for(it = pplot.begin() ; it != pplot.end() ; ++it, ++count) {
            glDisable(GL_BLEND);
            glPointSize(0.5);
            const nvis::fvec3& c = it->second;
            glColor3f(c[0], c[1], c[2]);
            glVertexPointer(2, GL_DOUBLE, 0, &it->first.front());
            glDrawArrays(GL_POINTS, 0, it->first.size());
        }
    }
    
#if 0
    glDisable(GL_BLEND);
    glPointSize(1.0);
    glColor3f(0.0, 0.0, 0.5);
    glDisable(GL_POINT_SMOOTH);
    glVertexPointer(2, GL_DOUBLE, 0, &dpts.front());
    glDrawArrays(GL_POINTS, 0, dpts.size());
    
    glDisable(GL_BLEND);
    glLineWidth(1.0);
    glColor3f(0.4, 0.2, 0.2);
    glDisable(GL_LINE_SMOOTH);
    glVertexPointer(2, GL_DOUBLE, 0, &lpts.front());
    glDrawArrays(GL_LINES, 0, lpts.size());
    
    // draw newton iteration points
    glEnable(GL_BLEND);
    glLineWidth(2.0);
    glColor3f(0.5, 0.5, 0.5);
    glEnable(GL_LINE_SMOOTH);
    glVertexPointer(2, GL_DOUBLE, 0, &npts.front());
    glDrawArrays(GL_LINE_STRIP, 0, npts.size());
    
    glDisableClientState(GL_VERTEX_ARRAY);
    
    // draw fixpoints
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glPointSize(8.0);
    glLineWidth(1.0);
    for(std::vector<fixpoint>::iterator fi = fpts.begin(); fi != fpts.end(); ++fi) {
        if(fi->saddle) {
            glColor3f(1, 0, 0);
            
            vec2 p00 = fi->pos + 0.01 * fi->evec[0] / norm(fi->evec[0]);
            vec2 p01 = fi->pos - 0.01 * fi->evec[0] / norm(fi->evec[0]);
            vec2 p10 = fi->pos + 0.01 * fi->evec[1] / norm(fi->evec[1]);
            vec2 p11 = fi->pos - 0.01 * fi->evec[1] / norm(fi->evec[1]);
            
            glBegin(GL_LINES);
            glVertex2dv((GLdouble*)&p00);
            glVertex2dv((GLdouble*)&p01);
            glVertex2dv((GLdouble*)&p10);
            glVertex2dv((GLdouble*)&p11);
            glEnd();
        } else {
            glColor3f(0, 0, 1);
        }
        
        glBegin(GL_POINTS);
        glVertex2dv((GLdouble*)&fi->pos);
        glEnd();
    }
#endif
    
    glutSwapBuffers();
}

// --------------------------------------------------------------------------------

void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(_bounds.min()[0], _bounds.max()[0],
            _bounds.min()[1], _bounds.max()[1],
            0, 1);
            
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// --------------------------------------------------------------------------------

void motion(int button, int state, int x, int y)
{
    // if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    //  GLint vp[4];
    //  glGetIntegerv(GL_VIEWPORT, vp);
    //
    //  GLdouble mv[16], pr[16];
    //
    //  glGetDoublev(GL_MODELVIEW_MATRIX, mv);
    //  glGetDoublev(GL_PROJECTION_MATRIX, pr);
    //
    //  double tmp;
    //  gluUnProject(x, vp[3] - (y - vp[1]), 0.0, mv, pr, vp, &cursor[0], &cursor[1], &tmp);
    //
    //  glutPostRedisplay();
    // }
}

// --------------------------------------------------------------------------------

void keyboard(unsigned char key, int, int)
{
    std::cout << key << '\n';
    
#if 0
    switch(key) {
        case 'p':
            poincare_plot(cursor);
            break;
        case 'f':
            find_fixpoint(cursor);
            break;
        case 's':
            find_seedpoints();
            break;
        case 'i':
            for(int i = 0; i < spts.size(); ++i) {
                cursor = spts[i].pos;
                find_fixpoint(cursor);
            }
            break;
        case 'C':
            lpts.clear();
            fpts.clear();
        case 'c':
            ppts.clear();
            npts.clear();
            spts.clear();
            break;
        case 'k':
            ++K;
            std::cout << "K=" << K << std::endl;
            break;
        case 'j':
            if(K > 1) {
                --K;
                std::cout << "K=" << K << std::endl;
            }
            break;
    }
#endif
    
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    initialize(argc, argv);
    
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(1500, 1500);
    glutInitWindowPosition(1400 - 768, 20);
    glutCreateWindow(argv[0]);
    
    init();
    
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(motion);
    glutMainLoop();
    return 0;
}


