#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>

#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <math/rational.hpp>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

// NRRD interface
#include <image/nrrd_wrapper.hpp>
// data structure
#include <data/grid.hpp>
#include <data/raster_data.hpp>
// display
#include <graphics/colors.hpp>
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>
#include "logical2physical.hpp"
// poincare map API
#include "xmt_poincare_map.hpp"
#include "map_field_wrapper.hpp"
// topological analysis
#include "map_analysis.hpp"
#include "newton.hpp"
#include "fixpoints.hpp"
// math
#include <math/rational.hpp>

#if _OPENMP
#include <omp.h>
#endif


using namespace xavier;

// -------------------------
//
//      Data Structure
//
// -------------------------
typedef grid<double, 2>                                     plane_type;
typedef std::pair<nvis::vec2, int>                          map_data_type;
typedef raster_data<map_data_type, double, 2>               dataset_type;
typedef standard_map                                        map_type;

plane_type*                 _plane;
nvis::ivec2                 _plane_res;
nvis::bbox2                 _plane_bounds;
xavier::map_metric          _plane_metric;
nvis::vec2                  _plane_spacing;
dataset_type*               _dataset;


// -------------------------
//
//              UI
//
// -------------------------
char*    in, *out;
int     niter, maxp, width, height, show_vec, __res[2], jitt, nogfx;
double  eps, K;

struct iteration_data {
    iteration_data() : sf(0) {}
    
    iteration_data(const std::vector<nvis::vec2> _s) {
        std::copy(_s.begin(), _s.end(), std::back_inserter(steps));
        sf = safety_factor(steps, _plane_metric);
    }
    
    nvis::vec2 vector(int p) const {
        assert(p<niter);
        return _plane_metric.displacement(steps[0], steps[p]);
    }
    
    std::vector<nvis::vec2> steps;
    double sf;
};

std::vector<nvis::vec2>     vectors;

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
    hestOptAdd(&hopt, "k",      "K",                    airTypeDouble,  0,  1,  &K,         "0.5",          "Constant controlling chaos");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1,  1,  &out,       NULL,           "output file name containig map vector approximation (NRRD)");
    hestOptAdd(&hopt, "n",      "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "e",      "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "r",      "resolution",           airTypeInt,     2,  2,  &__res,     NULL,           "plane resolution");
    hestOptAdd(&hopt, "j",      "# jitter",             airTypeInt,     0,  1,  &jitt,      "3",            "number of jittered samples");
    hestOptAdd(&hopt, "vec",    "show vectors",         airTypeInt,     0,  0,  &show_vec,  NULL,           "show vertex vectors");
    hestOptAdd(&hopt, "nogfx",  "no display",           airTypeInt,     0,  0,  &nogfx,     NULL,           "turn off display");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Standard map to vector field conversion",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

// -------------------------
//
//          Display
//
// -------------------------
typedef std::vector<nvis::vec2>                 orbit_type;
typedef std::pair<orbit_type, nvis::fvec3>      color_orbit_type;

std::vector<color_orbit_type>           _orbits;
xavier::discrete_color_map<int>*        _cmap;
xavier::logical2physical*               _converter;
xavier::map_analysis_param              _params;
std::vector<nvis::ivec2>                _saddle_cells, _center_cells;
nvis::vec2                              _last;
std::vector<nvis::vec2>                 _problematic_seeds;
std::vector<color_orbit_type>           _rational_surfaces;

bool verbose;
inline double __safety_factor(const std::vector<nvis::vec2>& steps,
                              const map_metric& metric)
{

    if (steps.size()<2) {
        return 0;
    }
    
    std::vector< double > q(steps.size()-1);
    std::vector< double > d(steps.size()-1);
    const nvis::vec2& x0 = steps[0];
    for (unsigned int i = 1 ; i < steps.size() ; ++i) {
        double dtheta = (steps[i]-x0)[0]/(steps[i]-x0)[1];
        if (dtheta < 1.0e-9) {
            return 10000;
        }
        q[i-1] = (double)i / dtheta * metric.width();
        d[i-1] = metric.distance(x0, steps[i]);
        if (d[i-1] == 0) {
            return q[i-1];
        }
        if (d[i-1] < 1.0e-6) {
            return q[i-1];
            std::cerr << "d(" << x0 << ") = " << d[i-1] << std::endl;
            d[i-1] = 1.0e-6;
        }
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));
    
    double total_q = 0;
    double total_w = 0;
    for (int i = 1 ; i < steps.size() ; ++i) {
        double w = min / fabs(d[i-1]);
        total_w += w;
        total_q += w * q[i-1];
    }
    
    return total_q/total_w;
}


double sf(const nvis::vec2& x0, const map_type& pmap, int niter)
{
    map_type* amap = pmap.clone();
    std::vector<nvis::vec2> tmp;
    nvis::vec2 y(x0[0], x0[1]+1);
    try {
        amap->map(y, tmp, niter);
    } catch(...) {
        return 0;
    }
    if (tmp.size() < niter) {
        return 0;
    }
    
    std::vector<nvis::vec2> steps;
    steps.push_back(x0);
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(steps));
    
    if (verbose) {
        std::cerr << "orbit at " << x0 << " contains:\n";
        for (int i=0 ; i<steps.size() ; ++i) {
            std::cerr << i << ": " << steps[i] << std::endl;
        }
    }
    
    return __safety_factor(steps, _plane_metric);
}

void compute_map(const grid<double, 2>& mesh, const map_type& pmap,
                 const int niter)
{
    typedef grid<double, 2>                         grid_type;
    typedef std::pair<nvis::vec2, int>              data_type;
    typedef raster_data<data_type, double, 2>       dataset_type;
    typedef grid_type::bounds_type                  bounds_type;
    typedef grid_type::vec_type                     vec_type;
    typedef grid_type::ivec_type                    ivec_type;
    
    const vec_type& step        = mesh.spacing();
    const bounds_type& bounds   = mesh.bounds();
    const ivec_type& resolution = mesh.dimensions();
    int npoints = resolution[0]*resolution[1];
    
    srand48(time(0));
    
    double* __data = (double*)calloc(npoints, sizeof(double));
    
    std::cerr << "bounds = " << bounds << std::endl;
    std::cerr << "step = " << step << std::endl;
    std::cerr << "resolution = " << resolution << std::endl;
    std::cerr << "jittering = " << jitt << " x " << jitt << std::endl;
    
    npoints = resolution[0]*resolution[1];
    vectors.resize(npoints);
    
    nvis::timer _timer;
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
// #pragma omp parallel
    {
// #pragma omp for schedule(dynamic,1)
        for(int n = 0 ; n < npoints ; ++n) {
            int thread_id = 0;
// #if _OPENMP
            // thread_id = omp_get_thread_num();
// #endif
            verbose = false;
            
            int j = n / resolution[0];
            int i = n % resolution[0];
            nvis::ivec2 id(i,j);
            
            if ((j==0 || j==resolution[1]-1) && i==100) {
                verbose = true;
            }
            
            if (thread_id == 0) {
                std::ostringstream os;
                os << "\rcompleted " << n << " / " << npoints << " ("
                   << 100*n/npoints << "%)      \r" << std::flush;
                std::cout << os.str();
            }
            
            nvis::vec2 x0 = (*_plane)(id);
            map_type* amap = pmap.clone();
            __data[n] = sf(x0, *amap, niter);
        }
    }
    std::cerr << "computation of " << npoints << " orbits took " << _timer.elapsed() << " s." << std::endl;
    
    size_t __size[2] = {resolution[0], resolution[1]};
    double __spc[2] = {_plane_spacing[0], _plane_spacing[1]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, __data, nrrdTypeDouble, 2, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave(out, nout, NULL);
    std::cerr << "values exported\n";
    nrrdNuke(nout);
}

static void init()
{
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    _plane_bounds = nvis::bbox2(nvis::vec2(0,0), nvis::vec2(1,1));
    _plane_res = nvis::ivec2(__res[0], __res[1]);
    _plane_metric.bounds() = _plane_bounds;
    _plane_metric.periodic(0) = true;
    _plane_metric.periodic(1) = true;
    int npoints = _plane_res[0] * _plane_res[1];
    
    GLUT_helper::box = _plane_bounds;
    
    map_type pmap(K);
    _plane = new plane_type(_plane_res, _plane_bounds);
    _dataset = new dataset_type(*_plane, map_data_type(nvis::vec2(0,0), 0));
    _params.nb_iterations = niter;
    _params.max_period = maxp;
    _params.metric = _plane_metric;
    _plane_spacing = _plane->spacing();
}

// --------------------------------------------------------------------------------

int main_window;
nvis::vec2 wc;

void idle(void)
{
    // switch context back to main window after GLUI activity
    glutSetWindow(main_window);
    glutPostRedisplay();
}

void draw(void)
{
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    const plane_type& mesh = *_plane;
    
    // draw vector associated with each vertex
    if (show_vec) {
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(1.0);
        const double cos_alpha = cos(M_PI / 12.);
        const double sin_alpha = sin(M_PI / 12.);
        glBegin(GL_LINES);
        for (int n=0 ; n<vectors.size() ; ++n) {
            int i = n % _plane_res[0];
            int j = n / _plane_res[0];
            nvis::ivec2 id(i,j);
            nvis::vec2 p0 = (*_plane)(id);
            const nvis::vec2& v = vectors[n];
            if (nvis::norm(v) == 0) {
                continue;
            }
            nvis::vec2 p1 = p0 + v;
            nvis::fvec3 c(0,0,1);
            glColor3f(c[0], c[1], c[2]);
            glVertex2f(p0[0], p0[1]);
            glVertex2f(p1[0], p1[1]);
            nvis::vec2 e0 = 0.2 * (p0 - p1);
            nvis::vec2 e1(-e0[1], e0[0]);
            glVertex2f(p1[0], p1[1]);
            nvis::vec2 y = p1 + cos_alpha * e0 + sin_alpha * e1;
            glVertex2f(y[0], y[1]);
            glVertex2f(p1[0], p1[1]);
            y = p1 + cos_alpha * e0 - sin_alpha * e1;
            glVertex2f(y[0], y[1]);
        }
        glEnd();
    }
}

void display(void)
{
    GLUT_helper::setup_display(draw);
    glutSwapBuffers();
}

void guiCallback(int)
{
    display();
}

bool save_to_file;
void mouse(int button, int state, int x, int y)
{
    map_type pmap(K);
    
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    display();
    
    std::ostringstream os;
    os << "Cursor at " << wc;
    glutSetWindowTitle(os.str().c_str());
    GLUT_helper::glut_helper_mouse(button, state, x, y);
}

void keyboard(unsigned char key, int x, int y)
{
    if(key == 'r') {
        GLUT_helper::resetCamera();
    }
    
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    _last = nvis::vec2(-1,-1);
    
    initialize(argc, argv);
    
    width = 1200;
    height = 800;
    
    /*
        // initialize GLUT
        if (!nogfx) {
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
        }
    */
    
    init();
    map_type pmap(K);
    compute_map(*_plane, pmap, niter);
    
    if (false) {
        display();
        
        // set GLUT callback functions
        glutDisplayFunc(display);
        glutReshapeFunc(GLUT_helper::glut_helper_reshape);
        glutMouseFunc(mouse);
        glutMotionFunc(GLUT_helper::glut_helper_motion);
        glutKeyboardFunc(keyboard);
        
        // adjust camera to mesh
        GLUT_helper::resetCamera();
        
        // Enter GLUT event loop
        glutMainLoop();
    }
    
    return 0;
}
