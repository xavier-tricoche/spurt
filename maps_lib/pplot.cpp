#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>

#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>
#include <kdtree++/kdtree.hpp>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

#include <teem/nrrd.h>

#include <data/grid.hpp>
#include <data/raster.hpp>
#include <data/kdtree.hpp>
#include <graphics/colors.hpp>
#include <image/nrrd_wrapper.hpp>
#include <maps_lib/display.hpp>
#include <maps_lib/period.hpp>
#include <math/fixed_vector.hpp>
#include "definitions.hpp"
#include "xmt_poincare_map.hpp"
#include "map_field_wrapper.hpp"
#include <math/divergence_cleaning.hpp>

// GLUT business
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>

#if _OPENMP
#include <omp.h>
#endif

using namespace spurt;
using namespace map_analysis;
using namespace map_display;
// using namespace div_cleaning;

using namespace GLUT_helper;
int main_window;
nvis::vec2 wc;

const double max_angle = 90;
const int max_count = 20;

namespace {
double __mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}
}

std::vector<std::pair<nvis::vec2, nvis::fvec3> > seeds;

typedef spurt::grid::uniform_grid<double, 3>  grid_type;

template<typename Val_> 
using image_type = image<Val_, 3, double, size_t>;

template<typename Val_>
class legacy_wrapper : public image_type<Val_> {
public:
    typedef image_type<Val_> base_type;
    typedef Val_             value_type;
    
    legacy_wrapper(base_type& image, const default_metric_type& metric) 
        : base_type(image), _metric(metric), _verbose(false) {}
        
    legacy_wrapper(const grid_type& grid, const std::vector<value_type>& val) 
        : base_type(grid, val), _verbose(false) {}
        
    grid_type get_grid() const {
        return grid_type(this->grid());
    }
    
    value_type interpolate(const typename base_type::point_type& p) const {
        return this->value(p);
    }
    
    void verbose(bool v) const {
        _verbose = v;
    }
    
private:
    spurt::default_metric_type    _metric;
    mutable bool          _verbose;
};

typedef legacy_wrapper<nvis::vec3>              trilinear_field_type;
typedef divfree_field<trilinear_field_type>     divfree_field_type;
typedef nvis::ivec3                             ivec_type;

typedef xmt_poincare_map<spurt::map::wrapper<trilinear_field_type> > trilinear_map_type;
typedef xmt_poincare_map<spurt::map::wrapper<divfree_field_type> >   divfree_map_type;

spurt::default_metric_type  metric2d;
nvis::bbox2 _bounds;

char*    in, *seed_l, *seed_s, *out;
int     npts, niter, rx, dolink;
double  eps;
bool    div_free;
float   pt_sz, ln_w;

struct vertex_data {
    vertex_data() : chain_id(-1), id(-1), next(-1), prev(-1) {}
    vertex_data(int ci, int i, int n=-1, int p=-1) : chain_id(ci), id(i), next(n), prev(p) {}
    vertex_data(const vertex_data& vd) : chain_id(vd.chain_id), id(vd.id), next(vd.next), prev(vd.prev) {}
    
    int chain_id, id, next, prev;
};

typedef spurt::point_locator<double, vertex_data, 2>   point_locator_type;
typedef spurt::data_point<double, vertex_data, 2>      data_point_type;
point_locator_type locator;

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
    hestOptAdd(&hopt, "edge", "create edges",       airTypeInt,     0,  1,  &dolink,    "0",            "connect points to form curves (0: none, 1: distance, 2: check, 3: period)");
    hestOptAdd(&hopt, "e",  "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "df", "div free",             airTypeBool,    0,  1,  &div_free,  "0",            "divergence free interpolation?");
    hestOptAdd(&hopt, "rx", "resolution factor",    airTypeInt,     0,  1,  &rx,        "10",           "output resolution upsample");
    hestOptAdd(&hopt, "ps", "point size",           airTypeFloat,   0,  1,  &pt_sz,     "2",            "point size for display");
    hestOptAdd(&hopt, "lw", "line width",           airTypeFloat,   0,  1,  &ln_w,      "1",            "line width for display");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Poincare map visualization",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef std::list<nvis::vec2> orbit_type;
struct colored_orbit : public std::pair<orbit_type, nvis::fvec3> {
    colored_orbit() {
        this->first.clear();
        this->second = nvis::fvec3();
    }
    colored_orbit(const colored_orbit& co) {
        this->first = orbit_type(co.first.begin(), co.first.end());
        this->second = co.second;
    }
};

typedef colored_orbit colored_orbit_type;


typedef std::list<std::pair<nvis::vec2, nvis::vec2> >   orbit_edges_type;
typedef std::pair<orbit_edges_type, nvis::fvec3>        colored_orbit_edges_type;

std::list<colored_orbit_type>               pplot;      // poincare plot
std::list<colored_orbit_edges_type>         pplot_links;
std::list<nvis::vec2>                       loose_ends;

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

#if 0
typedef boost::rational<int>    rational_type;

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
        c += field.distance(steps[i], steps[i+period]);
    }
    return c / (double)(steps.size() - period);
}

inline double __angle(const nvis::vec2& x, const nvis::vec2& y)
{
    double cosine = nvis::inner(x, y) / (nvis::norm(x) * nvis::norm(y));
    return acos(cosine) / M_PI * 180.;
}

struct __vertex {
    __vertex() : id(-1), next(-1), prev(-1) {}
    __vertex(int i) : id(i), next(-1), prev(-1) {}
    bool empty() const {
        return next < 0;
    }
    bool full() const {
        return (next >= 0 && prev >= 0);
    }
    void add(int i, const std::vector<nvis::vec2>& pts) {
        if (full()) {
            throw std::runtime_error("attempt to add edge to fully connected vertex");
        }
        if (!empty()) {
            prev = i;
            prev_dir = metric2d.displacement(pts[id], pts[i]);
        } else {
            next = i;
            next_dir = metric2d.displacement(pts[id], pts[i]);
        }
    }
    bool compatible(const nvis::vec2& dir) const {
        if (empty()) {
            return true;
        } else if (full()) {
            return false;
        } else {
            return (__angle(next_dir, dir) < max_angle);
        }
    }
    
    int id, next, prev;
    nvis::vec2 next_dir, prev_dir;
};

template<typename M>
void connect(const std::vector<nvis::vec2>& steps, const M& pmap)
{
    // PARAMS: max_angle, max_count;
    std::cerr << "There are " << steps.size() << " points in orbit\n";
    
    std::vector<__vertex> vertices(steps.size());
    for (int i=0 ; i<steps.size() ; ++i) {
        vertices[i].id = i;
    }
    
    const bool periodic[] = {true, false};
    std::ostringstream os;
    
    std::map<double, int> quality;
    for(unsigned int i = 1 ; i < steps.size() / 3 ; ++i) {
        double c = cost(steps, i, pmap.field());
        quality.insert(std::pair<double, int>(c, i));
        if(c == 0) {
            os << "WARNING: zero cost for " << steps.size() << " steps" << std::endl;
            std::cerr << os.str();
        }
    }
    if(quality.empty()) {
        return;
    }
    
    
    for (int i=0 ; i<steps.size() ; ++i) {
        if (vertices[i].full()) {
            std::cerr << "vertex #" << i << " is already fully connected\n";
            continue;
        }
        
        // consider a fixed number of samples from best period approximations
        int count = 0;
        std::map<double, int> loc_qual;
        
        os.clear();
        os.str("");
        os << "\nvertex #" << i << " at " << metric2d.modulo(steps[i]) << " is considering periods:\n";
        for (std::map<double, int>::const_iterator it=quality.begin() ;
                it!=quality.end() && count<max_count ; ++it, ++count) {
            int p = it->second;
            os << "p (" << it->first << ")\n";
            // apply period forward...
            if (i+p<steps.size()) {
                nvis::vec2 dir = metric2d.displacement(steps[i], steps[i+p]);
                double dist = nvis::norm(dir);
                os << "\tfwd: " << i+p << " at " << metric2d.modulo(steps[i+p]) << ", distance " << dist;
                if (vertices[i+p].compatible(dir) && vertices[i].compatible(-1.*dir)) {
                    loc_qual.insert(std::pair<double, int>(dist, i+p));
                    os << ", valid\n";
                } else {
                    os << ", invalid\n";
                }
            }
            // ...and backward
            if (i>p) {
                nvis::vec2 dir = metric2d.displacement(steps[i], steps[i-p]);
                double dist = nvis::norm(dir);
                os << "\tbwd: " << i-p << " at " << metric2d.modulo(steps[i-p]) << ", distance " << dist;
                if (vertices[i-p].compatible(dir) && vertices[i].compatible(-1*dir)) {
                    loc_qual.insert(std::pair<double, int>(dist, i-p));
                    os << ", valid\n";
                } else {
                    os << ", invalid\n";
                }
            }
        }
        
        std::cerr << os.str();
        
        if (loc_qual.empty()) {
            std::cerr << "empty list of potential neighbors!";
            continue;
        }
        
        vertices[i].add(loc_qual.begin()->second, steps);
        vertices[loc_qual.begin()->second].add(i, steps);
        
        std::cerr << "best match found is #" << loc_qual.begin()->second << " at "
                  << metric2d.modulo(steps[loc_qual.begin()->second]) << '\n';
                  
        // loose end check
        if (!vertices[i].full()) {
        
            std::cerr << "looking for another edge...\n";
            
            std::map<double, int>::const_iterator jt = loc_qual.begin();
            for (++jt ; jt!=loc_qual.end() ; ++jt) {
                // is this second edge compatible with the first one we selected?
                nvis::vec2 dir = metric2d.displacement(steps[i], steps[jt->second]);
                std::cerr << "considering #" << jt->second << " at distance " << nvis::norm(dir) << std::endl;
                if (vertices[i].compatible(-1*dir)) {
                    vertices[i].add(jt->second, steps);
                    vertices[jt->second].add(i, steps);
                    std::cerr << "valid\n";
                } else {
                    std::cerr << "invalid\n";
                    if (!vertices[i].compatible(-1.*dir)) {
                        std::cerr << "edge direction " << dir << " is incompatible with current edge " << vertices[i].next_dir << "\n";
                        std::cerr << "angle between the two is " << __angle(-1.*dir, vertices[i].next_dir) << "\n";
                    }
                }
            }
        }
    }
    
    pplot_links.push_back(colored_orbit_edges_type());
    colored_orbit_edges_type& coe = pplot_links.back();
    orbit_edges_type& orb = coe.first;
    for (int i=0 ; i<vertices.size() ; ++i) {
        nvis::vec2 xi = metric2d.modulo(steps[i]);
        if (!vertices[i].empty()) {
            int j = vertices[i].next;
            nvis::vec2 xj = xi + metric2d.displacement(steps[i], steps[j]);
            if (i<j) {
                std::vector<std::pair<nvis::vec2, nvis::vec2> > segments;
                metric2d.clip_segment(segments, xi, xj);
                for(int j = 0 ; j < segments.size() ; ++j) {
                    orb.push_back(segments[j]);
                }
            }
            if (vertices[i].full()) {
                j = vertices[i].prev;
                nvis::vec2 xj = xi + metric2d.displacement(steps[i], steps[j]);
                if (i<j) {
                    std::vector<std::pair<nvis::vec2, nvis::vec2> > segments;
                    metric2d.clip_segment(segments, xi, xj);
                    for(int j = 0 ; j < segments.size() ; ++j) {
                        orb.push_back(segments[j]);
                    }
                }
            }
        }
    }
    coe.second = hsv2rgb(drand48() * 360., 1, 0.5);
}
#endif

void connect(const std::vector<nvis::vec2>& steps)
{
    std::vector<std::vector<int> > curves;
    if (dolink == 1) {
        map_analysis::connect_by_distance(curves, steps, metric2d);
    } else if (dolink == 2) {
        map_analysis::connect_and_check(curves, steps, metric2d);
    } else if (dolink == 3) {
        map_analysis::connect_by_period(curves, steps, metric2d);
    } else {
        return;
    }
    
    pplot_links.push_back(colored_orbit_edges_type());
    colored_orbit_edges_type& coe = pplot_links.back();
    orbit_edges_type& orb = coe.first;
    
    for (int ci = 0 ; ci<curves.size() ; ++ci) {
        for (int i=0 ; i<curves[ci].size()-1 ; ++i) {
            nvis::vec2 x = metric2d.modulo(steps[curves[ci][i]]);
            nvis::vec2 y = metric2d.modulo(steps[curves[ci][i+1]]);
            std::vector<std::pair<nvis::vec2, nvis::vec2> > segments;
            metric2d.clip_segment(segments, x, y);
            for(int j = 0 ; j < segments.size() ; ++j) {
                orb.push_back(segments[j]);
            }
        }
    }
    coe.second = hsv2rgb(drand48() * 360., 1, 0.5);
}

void guiCallBack(int control);


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
            
            orbits[thread_id].push_back(colored_orbit_type());
            colored_orbit_type& orb = orbits[thread_id].back();
            for(unsigned int n = 0 ; n < steps.size() ; ++n) {
                orb.first.push_back(metric2d.modulo(steps[n]));
                data_point_type dp(metric2d.modulo(steps[n]), vertex_data(i, n));
                locator.insert(dp);
            }
            orb.second = seeds[i].second;
            
            if(dolink) {
                connect(steps);
            }
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
    nin = spurt::nrrd_utils::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> __array;
    spurt::nrrd_utils::to_vector(__array, nin);
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
    basic_field.verbose(false);
    
    double h = eps;
    
    trilinear_map_type basic_map(basic_field);
    basic_map.precision(eps);
    
    divfree_field<trilinear_field_type> divfree_field(domain, basic_field);
    divfree_field.verbose(false);
    
    divfree_map_type divfree_map(divfree_field);
    divfree_map.precision(eps);
    
    grid_type::bounds_type bbox(domain.bounds());
    _bounds = nvis::bounding_box<nvis::vec2> (nvis::vec2(bbox.min()[1], bbox.min()[0]),
              nvis::vec2(bbox.max()[1], bbox.max()[0]));
              
    metric2d.bounds() = _bounds;
    metric2d.periodic(0) = true;
    metric2d.periodic(1) = false;
    
    GLUT_helper::box = _bounds;
    std::cerr << "box set to " << GLUT_helper::box << std::endl;
    
    if(!div_free) {
        poincare_plot(pplot, basic_map, npts, niter);
    } else {
        poincare_plot(pplot, divfree_map, npts, niter);
    }
    
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
}

// --------------------------------------------------------------------------------

void idle(void)
{
    // switch context back to main window after GLUI activity
    glutSetWindow(main_window);
    glutPostRedisplay();
}

void draw(void)
{
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    
    if(dolink) {
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(ln_w);
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
    
    if (true) {
        // draw poincare plot - each orbit in a different color
        std::list<colored_orbit_type>::iterator it;
        unsigned int count = 0;
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_BLEND);
        glPointSize(pt_sz);
        glBegin(GL_POINTS);
        for(it = pplot.begin() ; it != pplot.end() ; ++it, ++count) {
            const nvis::fvec3& c = it->second;
            glColor3f(c[0], c[1], c[2]);
            for (orbit_type::const_iterator jt=it->first.begin () ; jt!=it->first.end() ; ++jt) {
                glVertex3f((*jt)[0], (*jt)[1], 0);
            }
        }
        glEnd();
    }
    
    // if (true) {
    //  glEnableClientState(GL_VERTEX_ARRAY);
    //
    //  glDisable(GL_POINT_SMOOTH);
    //  glDisable(GL_BLEND);
    //  glPointSize(4);
    //  glColor3f(1, 0, 0);
    //  glVertexPointer(2, GL_DOUBLE, 0, &loose_ends.front());
    //  glDrawArrays(GL_POINTS, 0, loose_ends.size());
    // }
}

void display(void)
{
    setup_display(draw);
    glutSwapBuffers();
}

void guiCallback(int)
{
    display();
}

void mouse(int button, int state, int x, int y)
{
    nvis::vec3 _wc = world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    std::ostringstream os;
    os << "Poincare plot - (" << wc[0] << ", " << GLUT_helper::box.max()[1] - wc[1] << ")";
    try {
        vertex_data pd = locator.find_close_point(wc);
        os << " - [" << pd.chain_id << ", " << pd.id << "]";
    } catch(...) {
    }
    glutSetWindowTitle(os.str().c_str());
    glut_helper_mouse(button, state, x, y);
}

void keyboard(unsigned char key, int x, int y)
{
    if(key == 'r') {
        resetCamera();
    } else if(key == 'x') {
        // output position coordinates
    }
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    initialize(argc, argv);
    
    width = 1200;
    height = 800;
    
    // initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(width, height);
    glutInitWindowPosition(20, 20);
    main_window = glutCreateWindow(argv[0]);
    
    // configure OpenGL for aesthetic results
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_LINE_SMOOTH, GL_NICEST);
    glHint(GL_POINT_SMOOTH, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    init();
    display();
    
    // set GLUT callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(glut_helper_reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(glut_helper_motion);
    glutKeyboardFunc(keyboard);
    
    // adjust camera to mesh
    resetCamera();
    resetCamera();
    
    // Enter GLUT event loop
    glutMainLoop();
    
    return 0;
}
