#include <iostream>
#include <map>
#include <list>
#include <vector>

#include <math/fixed_vector.hpp>

#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <tokamak/poincare_map.hpp>

#include "maps_lib/adaptive_map_sampling.hpp"
#include "maps_lib/definitions.hpp"
#include "maps_lib/orbits.hpp"
#include "maps_lib/manifold.hpp"

#include <teem/hest_helper.hpp>

// GLUT business
#include "graphics/GLUT_helper.hpp"
#include "graphics/GUI/GLUI_Wrapper.h"

using namespace map_analysis;
using namespace GLUT_helper;

// interface to the computation
adaptive_map_sampling_params params;
adaptive_map_sampling_output output;

const nvis::fvec3 red(1, 0, 0);
const nvis::fvec3 green(0, 1, 0);
const nvis::fvec3 blue(0, 0, 1);
const nvis::fvec3 yellow(1, 1, 0);
const nvis::fvec3 cyan(0, 1, 1);
const nvis::fvec3 white(1, 1, 1);
const nvis::fvec3 black(0, 0, 0);
const nvis::fvec3 orange(1, 0.5, 0);
const nvis::fvec3 magenta(1, 0, 1);
// some funky color names after Apple's color editor
const nvis::fvec3 cayenne(0.5, 0, 0);
const nvis::fvec3 midnight(0, 0, 0.5);
const nvis::fvec3 aqua(0, 0.5, 1);
const nvis::fvec3 sea_green(0, 1, 0.5);
const nvis::fvec3 lime(0.5, 1, 0);
const nvis::fvec3 framboise(1, 0, 0.5);
const nvis::fvec3 carnation(1, 0.5, 1);

const nvis::fvec3 rainbow[] = {
    black, midnight, blue, aqua,
    cyan, sea_green, green, lime,
    yellow, orange, red, framboise,
    magenta, carnation, white
};

template<typename T>
struct adaptive_color_map {
    adaptive_color_map(const std::vector<T>& vals, const std::vector<nvis::fvec3>& scale)
        : t(scale.size()), colors(scale) {
        std::vector<T> tmp(vals);
        std::sort(tmp.begin(), tmp.end());
        unsigned int step = tmp.size() / (scale.size() - 1);
        for (int i = 0 ; i < t.size() ; ++i) {
            t[i] = tmp[i*step];
        }
        t.back() = tmp.back();
    }
    
    nvis::fvec3 operator()(const T& val) const {
        unsigned int bin = std::distance(t.begin(),
                                         std::lower_bound(t.begin(), t.end(), val));
        if (bin > t.size() - 2) {
            bin = t.size() - 2;
        }
        T u = (val - t[bin]) / (t[bin+1] - t[bin]);
        return (1. - u)*colors[bin] + u*colors[bin+1];
    }
    
    std::vector<T>          t;
    std::vector<nvis::fvec3> colors;
};

// add new graphical object types here
enum { ALL_MESH, CAND_MESH, SING_MESH, PROB_MESH, PROB_EDGES,
       HEDGEHOG, PPLOT, INDEX_VECTORS, FIXPOINTS, DEGENERATE_POINTS,
       SEPARATRICES, MANIFOLD
     };

struct gfx_ID {
    gfx_ID() : obj_type(-1), period(-1), list_id(0) {}
    gfx_ID(int type, unsigned int p)
        : obj_type(type), period(p), list_id(0) {}
    gfx_ID(const gfx_ID& gID)
        : obj_type(gID.obj_type), period(gID.period), list_id(gID.list_id) {}
        
    bool operator<(const gfx_ID& gID) const {
        if (obj_type < gID.obj_type) {
            return true;
        } else if (obj_type > gID.obj_type) {
            return false;
        }
        return (period < gID.period);
    }
    
    int obj_type;
    unsigned int period;
    GLuint list_id;
};

std::set<gfx_ID> all_primitives;

gfx_ID register_id(gfx_ID& gID)
{
    std::set<gfx_ID>::iterator it = all_primitives.find(gID);
    if (it == all_primitives.end()) {
        // this display list does not already exist - create it
        gID.list_id = glGenLists(1);
        all_primitives.insert(gID);
    } else {
        // this display list exists - overwrite it by giving its ID away
        gID.list_id = it->list_id;
    }
}

// computation parameters
int n_samples, n_iter, max_tri, max_per, min_per;
char* fout, *fin, *ts;
double eps, dq, min_area, max_area, max_ratio, max_err, dx, cl;
void init(int argc, char* argv[])
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
    hestOptAdd(&hopt, "f",      "file",                 airTypeString,  1, 1, &fin,         NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",                 airTypeString,  1, 1, &ts,          NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &fout,        NULL,       "output nrrd file name");
    hestOptAdd(&hopt, "eps",    "eps",                  airTypeDouble,  0, 1, &eps,         "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "n",      "# samples",            airTypeInt,     0, 1, &n_samples,   "5000",     "number of 1D samples");
    hestOptAdd(&hopt, "m",      "# rotations",          airTypeInt,     0, 1, &n_iter,      "50",       "number of rotations");
    hestOptAdd(&hopt, "dq",     "safety factor step",   airTypeDouble,  0, 1, &dq,          "0.05",     "maximum discrepancy between consecutive safety factors");
    hestOptAdd(&hopt, "mr",     "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,   "1.",       "maximum triangle aspect ratio");
    hestOptAdd(&hopt, "mt",     "max # triangles",      airTypeInt,     0, 1, &max_tri,     "150000",   "max number of triangles in adaptive sampling");
    hestOptAdd(&hopt, "ma",     "min triangle area",    airTypeDouble,  0, 1, &min_area,    "0.01",     "min triangle area in adaptive sampling");
    hestOptAdd(&hopt, "Ma",     "max triangle area",    airTypeDouble,  0, 1, &max_area,    "1.",       "max triangle area in adaptive sampling");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0, 1, &max_per,     "15",       "max considered period in fixed point search");
    hestOptAdd(&hopt, "pm",     "minx period",          airTypeInt,     0, 1, &min_per,     "1",        "min considered period in fixed point search");
    hestOptAdd(&hopt, "err",    "max tolerance",        airTypeDouble,  0, 1, &max_err,     "1.57",     "max tolerance of approximation quality criterion");
    hestOptAdd(&hopt, "dx",     "min dist",             airTypeDouble,  0, 1, &dx,          "1.0e-4",   "min allowable sampling distance in index computation");
    hestOptAdd(&hopt, "cd",     "close dist",           airTypeDouble,  0, 1, &cl,          "1.0",      "min allowable distance between chains of same type");
    
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Adaptively sample phase portrait along orbits to achieve accurate piecewise linear interpolation on a per-period basis",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T>
void draw_poincare_plot(const T& color_map)
{
    gfx_ID plot_id(PPLOT, 0);
    register_id(plot_id);
    glNewList(plot_id.list_id, GL_COMPILE);
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(1);
    glBegin(GL_POINTS);
    for (int i = 0 ; i < static_data::central_map_orbits.size() ; ++i) {
        const orbit& orb = static_data::central_map_orbits[i];
        nvis::fvec3 col = color_map(orb.period());
        glColor3f(col[0], col[1], col[2]);
        for (int j = 0 ; j < orb.size() ; ++j) {
            nvis::vec2 x = static_data::metric.modulo(orb[j]);
            glVertex2f(x[0], x[1]);
        }
    }
    glEnd();
    glEndList();
}

template<typename T, typename I>
void draw_mesh(const T& mesh, const nvis::fvec3& color, gfx_ID& id, const I& included)
{
    typedef typename T::triangle_type   triangle_type;
    
    // uniquify edges
    std::set<xavier::Edge> edges;
    for (int i = 0 ; i < mesh.get_nb_triangles() ; ++i) {
        if (!included(mesh, i)) {
            continue;
        }
        const triangle_type& tri = mesh.get_triangle_vertices(i);
        edges.insert(xavier::Edge(tri[0], tri[1]));
        edges.insert(xavier::Edge(tri[1], tri[2]));
        edges.insert(xavier::Edge(tri[2], tri[0]));
    }
    
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glBegin(GL_LINES);
    glColor3f(color[0], color[1], color[2]);
    for (std::set<xavier::Edge>::const_iterator it = edges.begin() ; it != edges.end() ; ++it) {
        const nvis::vec2& p0 = mesh.get_vertex(it->i0);
        const nvis::vec2& p1 = mesh.get_vertex(it->i1);
        glVertex3f(p0[0], p0[1], -1);
        glVertex3f(p1[0], p1[1], -1);
    }
    glEnd();
    glEndList();
}

template<typename T>
struct default_mesh_functor {
    bool operator()(const T&, unsigned int) const {
        return true;
    }
};

template<typename T>
struct selection_mesh_functor {
    selection_mesh_functor(const std::list<unsigned int>& included)
        : _included(included.begin(), included.end()) {}
        
    bool operator()(const T&, unsigned int tri) const {
        return _included.find(tri) != _included.end();
    }
    
    std::set<unsigned int> _included;
};

template<typename T>
void draw_mesh(const T& mesh, const nvis::fvec3& color, unsigned int period = 0)
{
    gfx_ID mesh_id(ALL_MESH, period);
    draw_mesh(mesh, color, mesh_id, default_mesh_functor<T>());
}

inline void draw_arrow(const nvis::vec2& p0, const nvis::vec2& p1)
{
    const double cos_alpha = cos(M_PI / 12.);
    const double sin_alpha = sin(M_PI / 12.);
    
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

template<typename T, typename F>
void draw_hedgehog(const T& mesh, const nvis::fvec3& color, const F& functor, double scaling)
{
    typedef typename T::data_type   data_type;
    
    size_t nb_pts = mesh.get_nb_vertices();
    
    gfx_ID id(HEDGEHOG, functor.period());
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
    for (int i = 0 ; i < nb_pts ; ++i) {
        const data_type& d = mesh.get_data(i);
        nvis::vec2 v = functor(d);
        if (nvis::norm(v) == 0) {
            continue;
        }
        v *= scaling;
        const nvis::vec2& p0 = mesh.get_vertex(i);
        nvis::vec2 p1 = p0 + v;
        draw_arrow(p0, p1);
    }
    glEnd();
    glEndList();
}

void draw_steps(const nvis::fvec3& color, unsigned int period, double scaling,
                const std::list<step_type>& steps)
{
    gfx_ID id(INDEX_VECTORS, period);
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
    for (list<step_type>::const_iterator it = steps.begin() ; it != steps.end() ; ++it) {
        draw_arrow(it->first, it->first + scaling*(it->second - it->first));
    }
    glEnd();
    glEndList();
}

template<typename T>
void draw_edges(const T& mesh, const nvis::fvec3& color, unsigned int period,
                const std::list<xavier::Edge>& edges)
{
    std::cerr << "drawing " << edges.size() << " " << period << "-edges with color"
              << color << "\n";
              
    gfx_ID id(PROB_EDGES, period);
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
    for (std::list<xavier::Edge>::const_iterator it = edges.begin() ; it != edges.end() ; ++it) {
        const typename T::point_type& p0 = mesh.get_vertex(it->i0);
        const typename T::point_type& p1 = mesh.get_vertex(it->i1);
        std::cerr << "drawing singular edge " << p0 << " - " << p1 << '\n';
        glVertex2f(p0[0], p0[1]);
        glVertex2f(p1[0], p1[1]);
    }
    glEnd();
    glEndList();
}

void draw_degenerate_points(const std::list<nvis::vec2>& points,
                            unsigned int period)
{
    gfx_ID id(DEGENERATE_POINTS, period);
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    
    std::cerr << "drawing " << points.size() << " degenerate points for period " << period << '\n';
    
    // draw points
    glColor3f(green[0], green[1], green[2]);
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(2);
    glBegin(GL_POINTS);
    for (std::list<nvis::vec2>::const_iterator it = points.begin();
            it != points.end(); ++it) {
        glVertex2f((*it)[0], (*it)[1]);
        std::cerr << *it << '\n';
    }
    glEnd();
    glEndList();
}

void draw_fixpoints(const std::list<std::list<xavier::fixpoint> >& chains, unsigned int size = 3)
{
    typedef std::list<xavier::fixpoint>     chain_type;
    typedef std::list<chain_type>           list_type;
    typedef list_type::const_iterator       it_type;
    
    gfx_ID id(FIXPOINTS, chains.front().front().K);
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(size);
    glBegin(GL_POINTS);
    glColor3f(1, 0, 0); // saddles in red
    for (it_type it = chains.begin() ; it != chains.end() ; ++it) {
        if (!it->front().saddle) {
            continue;
        }
        for (chain_type::const_iterator it2 = it->begin() ; it2 != it->end() ; ++it2) {
            glVertex3f(it2->pos[0], it2->pos[1], 0.001);    // makes fixpoint stick out
        }
    }
    glColor3f(0, 0, 1); // centers in blue
    for (it_type it = chains.begin() ; it != chains.end() ; ++it) {
        if (it->front().saddle) {
            continue;
        }
        for (chain_type::const_iterator it2 = it->begin() ; it2 != it->end() ; ++it2) {
            glVertex3f(it2->pos[0], it2->pos[1], 0.001);    // makes fixpoint stick out
        }
    }
    glEnd();
    
    // draw saddle's eigenvectors
    glBegin(GL_LINES);
    for (it_type it = chains.begin() ; it != chains.end() ; ++it) {
        if (!it->front().saddle) {
            continue;
        }
        for (chain_type::const_iterator it2 = it->begin() ; it2 != it->end() ; ++it2) {
            glColor3f(0, 0, 1); // stable eigenvectors in blue
            nvis::vec2 x = it2->pos + 0.5 * it2->evec[0];
            nvis::vec2 y = it2->pos - 0.5 * it2->evec[0];
            glVertex2f(x[0], x[1]);
            glVertex2f(y[0], y[1]);
            glColor3f(1, 0, 0); // unstable eigenvectors in red
            x = it2->pos + 0.5 * it2->evec[1];
            y = it2->pos - 0.5 * it2->evec[1];
            glVertex2f(x[0], x[1]);
            glVertex2f(y[0], y[1]);
        }
    }
    glEnd();
    glEndList();
}

template<typename Iterator>
void clip_curve(std::vector<std::vector<nvis::vec2> >& clipped,
                Iterator begin, Iterator end,
                const map_metric& metric)
{
    clipped.clear();
    clipped.push_back(std::vector<nvis::vec2>());
    nvis::vec2 prev = metric.modulo(*begin);
    clipped.back().push_back(prev);
    Iterator it = begin;
    ++it;
    for (; it != end ; ++it) {
        nvis::vec2 next = metric.modulo(*it);
        nvis::vec2 aux = prev + metric.displacement(prev, next);
        if (nvis::norm(aux - next) > 1.0e-6*metric.diameter()) {
            clipped.back().push_back(aux);
            clipped.push_back(std::vector<nvis::vec2>());
            clipped.back().push_back(next - metric.displacement(aux, next));
        }
        clipped.back().push_back(next);
        prev = next;
    }
}

void draw_separatrices(const std::list<std::pair<int, std::vector<separatrix_type> > >& sep,
                       unsigned int p)
{
    gfx_ID id(SEPARATRICES, p);
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glColor3f(0.8, 0.8, 0);
    for (std::list<separatrix_connection_type>::const_iterator it = sep.begin();
            it != sep.end() ; ++it) {
        for (int connection = 0 ; connection < it->second.size() ; ++connection) {
            std::vector<std::vector<nvis::vec2> > clipped;
            clip_curve(clipped, it->second[connection].begin(), it->second[connection].end(),
                       static_data::metric);
            for (int i = 0 ; i < clipped.size() ; ++i) {
                glBegin(GL_LINE_STRIP);
                for (int j = 0 ; j < clipped[i].size() ; ++j) {
                    glVertex2f(clipped[i][j][0], clipped[i][j][1]);
                }
                glEnd();
            }
        }
    }
    glEndList();
}

template<typename Map>
void compute_manifold(unsigned int period, double length, const Map& map)
{
    const map_metric& metric = map.metric();
    
    std::vector<fp_chain> all_p_chains;
    const std::list<std::list<xavier::fixpoint> >& chains = output.p_chains[period-1];
    
    typedef std::list<xavier::fixpoint> list_chain_type;
    for (std::list<list_chain_type>::const_iterator it = chains.begin() ;
            it != chains.end() ; ++it) {
        all_p_chains.push_back(fp_chain(it->begin(), it->end(), static_data::metric));
    }
    
    std::vector<std::vector<nvis::vec2> > manifolds, tmp;
    for (int i = 0 ; i < all_p_chains.size() ; ++i) {
        const fp_chain& ch = all_p_chains[i];
        if (!ch.saddle()) {
            continue;
        }
        
        manifold_progress pro1, pro2;
        manifold_type man1, man2;
        
        std::cerr << "drawing separatrix in forward direction\n";
        ManBVP(man1, pro1, map, ch[0], true, 0.001, length, 0.01*length, 0.3, 1.0e-5);
        std::cerr << man1.size() << " steps in output\n";
        
        std::cerr << "drawing separatrix in backward direction\n";
        ManBVP(man2, pro2, map, ch[0], false, 0.001, length, 0.01*length, 0.3, 1.0e-5);
        std::cerr << man2.size() << " steps in output\n";
        
        clip_curve(tmp, man1.begin(), man1.end(), metric);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(manifolds));
        clip_curve(tmp, man2.begin(), man2.end(), metric);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(manifolds));
    }
    
    gfx_ID id(MANIFOLD, period);
    register_id(id);
    glNewList(id.list_id, GL_COMPILE);
    glColor3f(0, 0.8, 0.8);
    for (int i = 0 ; i < manifolds.size() ; ++i) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0 ; j < manifolds[i].size() ; ++j) {
            glVertex2f(manifolds[i][j][0], manifolds[i][j][1]);
        }
        glEnd();
    }
    glEndList();
}

template<typename T>
struct p_map_functor {
    typedef typename T::data_type data_type;
    
    p_map_functor(unsigned int period, double dq)
        : _p(period), _dq(dq) {
        for (unsigned int d = 1 ; d <= _p ; ++d) {
            rational_type q(_p, d);
            if (q.numerator() != _p) {
                continue;
            } else {
                _qs.push_back(xavier::value(q));
            }
        }
    }
    
    unsigned int period() const {
        return _p;
    }
    
    nvis::vec2 operator()(const data_type& d) const {
        double q = d.period();
        std::vector<double> dist(_qs.size());
        for (int i = 0 ; i < _qs.size() ; ++i) {
            dist[i] = fabs(_qs[i] - q);
        }
        if (*std::min_element(dist.begin(), dist.end()) > _dq) {
            return nvis::vec2(0, 0);
        } else {
            return map_analysis::vector_value(d, _p, static_data::metric);
        }
    }
    
    std::vector<double> _qs;
    double              _dq;
    unsigned int        _p;
};

// GUI stuff
enum { CHANGE_HEDGEHOG_SCALING,
       CHANGE_MIN_PERIOD,
       CHANGE_MAX_PERIOD,
       CHANGE_FOCUS_PERIOD,
       CHANGE_SHOW_FOCUS_MESH,
       CHANGE_SHOW_CAND_MESH,
       CHANGE_SHOW_SING_MESH,
       CHANGE_SHOW_PROB_MESH,
       CHANGE_SHOW_PROB_EDGES,
       CHANGE_SHOW_DEGENERATE_POINTS,
       CHANGE_SHOW_HEDGEHOG,
       CHANGE_SHOW_SINGULARITIES,
       CHANGE_SHOW_SEPARATRICES,
       CHANGE_DRAW_MANIFOLD,
       CHANGE_MANIFOLD_LENGTH,
       CHANGE_SHOW_PPLOT,
       CHANGE_SHOW_INDEX_VECTORS,
       CHANGE_SHOW_FIXPOINTS,
       CHANGE_FIXPOINT_SIZE
     };

void guiCallBack(int control);

bool show_focus_mesh, show_cand_mesh, show_sing_mesh, show_prob_mesh, show_prob_edges;
bool show_hedgehog, show_plot, show_index_vectors, show_fixpoints, show_degenerate_points;
bool show_separatrices, draw_manifold;
float hedgehog_scale, manifold_length;
unsigned int min_period, max_period, focus_period, fp_size;
int main_window;
nvis::vec2 wc;

void idle()
{
    // switch context back to main window after GLUI activity
    glutSetWindow(main_window);
    std::ostringstream os;
    os << "Visual Analysis of Poincare Map - (" << wc[0] << ", " << wc[1] << ")";
    glutSetWindowTitle(os.str().c_str());
    glutPostRedisplay();
}

void guiCallback(int);

void init_glui()
{
    GLUI* glui = GLUI_Master.create_glui("", 0, GLUT_helper::width + 60, 50);
    
    show_hedgehog           = false;
    show_focus_mesh         = false;
    show_cand_mesh          = false;
    show_sing_mesh          = false;
    show_prob_mesh          = false;
    show_prob_edges         = false;
    show_fixpoints          = false;
    show_degenerate_points  = false;
    draw_manifold           = false;
    hedgehog_scale = 0.1;
    manifold_length = 10.;
    min_period = 1;
    max_period = max_per;
    focus_period = 1;
    fp_size = 3;
    show_plot = true;
    show_index_vectors = false;
    
    // poincare plot display control
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show poincare plot", &show_plot, CHANGE_SHOW_PPLOT, guiCallback);
    
    glui->add_separator();
    
    // focus period display control
    GLUI_Wrapper::addSpinner(glui, NULL, "Focus period", &focus_period, CHANGE_FOCUS_PERIOD, guiCallback);
    
    glui->add_separator();
    
    GLUI_Wrapper::clampLastAddedSpinner(1, max_per);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show focus mesh", &show_focus_mesh, CHANGE_SHOW_FOCUS_MESH, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show problem cells", &show_prob_mesh, CHANGE_SHOW_PROB_MESH, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show problem edges", &show_prob_edges, CHANGE_SHOW_PROB_EDGES, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show index vectors", &show_index_vectors, CHANGE_SHOW_INDEX_VECTORS, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show degenerate points", &show_degenerate_points, CHANGE_SHOW_DEGENERATE_POINTS, guiCallback);
    
    glui->add_separator();
    
    // period range display control
    GLUI_Wrapper::addSpinner(glui, NULL, "Min period", &min_period, CHANGE_MIN_PERIOD, guiCallback);
    GLUI_Wrapper::clampLastAddedSpinner(1, max_per);
    GLUI_Wrapper::addSpinner(glui, NULL, "Max period", &max_period, CHANGE_MAX_PERIOD, guiCallback);
    GLUI_Wrapper::clampLastAddedSpinner(1, max_per);
    
    glui->add_separator();
    
    // mesh + singularity display control
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show candidates", &show_cand_mesh, CHANGE_SHOW_CAND_MESH, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show singular", &show_sing_mesh, CHANGE_SHOW_SING_MESH, guiCallback);
    
    glui->add_separator();
    
    // hedgehog display control
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show hedgehog", &show_hedgehog, CHANGE_SHOW_HEDGEHOG, guiCallback);
    GLUI_Wrapper::addSpinner(glui, NULL, "Hedgehog scale", &hedgehog_scale, CHANGE_HEDGEHOG_SCALING, guiCallback);
    GLUI_Wrapper::clampLastAddedSpinner(0, 5);
    
    glui->add_separator();
    
    // topology display control
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show fixed points", &show_fixpoints, CHANGE_SHOW_FIXPOINTS, guiCallback);
    GLUI_Wrapper::addSpinner(glui, NULL, "Fixed point size", &fp_size, CHANGE_FIXPOINT_SIZE, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Show separatrices", &show_separatrices, CHANGE_SHOW_SEPARATRICES, guiCallback);
    GLUI_Wrapper::addCheckbox(glui, NULL, "Draw manifold", &draw_manifold, CHANGE_DRAW_MANIFOLD, guiCallback);
    GLUI_Wrapper::addSpinner(glui, NULL, "Manifold length", &manifold_length, CHANGE_MANIFOLD_LENGTH, guiCallback);
    
    glui->set_main_gfx_window(main_window);
    
    GLUI_Master.set_glutIdleFunc(idle);
}

// callback functions
void draw()
{
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    std::set<gfx_ID>::const_iterator it;
    for (it = all_primitives.begin(); it != all_primitives.end() ; ++it) {
        // period specific stuff
        if (it->period == focus_period) {
            bool display = true;
            switch (it->obj_type) {
                case ALL_MESH:
                    if (!show_focus_mesh) {
                        continue;
                    }
                    break;
                case INDEX_VECTORS:
                    if (!show_index_vectors) {
                        continue;
                    }
                    break;
                case PROB_MESH:
                    if (!show_prob_mesh) {
                        continue;
                    }
                    break;
                case PROB_EDGES:
                    if (!show_prob_edges) {
                        continue;
                    }
                    break;
                case DEGENERATE_POINTS:
                    if (!show_degenerate_points) {
                        continue;
                    }
                    break;
                default:
                    display = false;
                    break;
            }
            if (display) {
                glCallList(it->list_id);
            }
        }
        // range-based display
        if (!it->period || (it->period >= min_period && it->period <= max_period)) {
            // filter by type
            switch (it->obj_type) {
                case CAND_MESH:
                    if (!show_cand_mesh) {
                        continue;
                    }
                    break;
                case SING_MESH:
                    if (!show_sing_mesh) {
                        continue;
                    }
                    break;
                case HEDGEHOG:
                    if (!show_hedgehog) {
                        continue;
                    }
                    break;
                case PPLOT:
                    if (!show_plot) {
                        continue;
                    }
                    break;
                case FIXPOINTS:
                    if (!show_fixpoints) {
                        continue;
                    }
                    break;
                case SEPARATRICES:
                    if (!show_separatrices) {
                        continue;
                    }
                    break;
                case MANIFOLD:
                    if (!draw_manifold) {
                        continue;
                    }
                    break;
                default:
                    continue;
            }
            glCallList(it->list_id);
        }
    }
}

void display()
{
    setup_display(draw);
    glutSwapBuffers();
}

void mouse(int button, int state, int x, int y)
{
    nvis::vec3 tmp = world_coordinates(x, y);
    wc = nvis::vec2(tmp[0], tmp[1]);
    std::ostringstream os;
    os << "Visual Analysis of Poincare Map - (" << wc[0] << ", " << height - wc[1] << ")";
    glutSetWindowTitle(os.str().c_str());
    glut_helper_mouse(button, state, x, y);
}

void updateFocusPeriodPrimitives()
{
    gfx_ID gID;
    // update all focus period dependent display lists
    if (show_focus_mesh) {
        gID = gfx_ID(ALL_MESH, focus_period);
        if (focus_period <= output.p_meshes.size() &&
                all_primitives.find(gID) == all_primitives.end()) {
            // we have not drawn this one yet
            draw_mesh(output.p_meshes[focus_period-1], nvis::fvec3(0.1, 0.1, 0.1), focus_period);
            display();
        }
    }
    if (show_prob_mesh) {
        gID = gfx_ID(PROB_MESH, focus_period);
        if (all_primitives.find(gID) == all_primitives.end()) {
            // we have not drawn this one yet
            draw_mesh(output.p_meshes[focus_period-1], yellow, gID,
                      selection_mesh_functor<mesh_type>(output.p_prob_tris[focus_period-1]));
        }
    }
    if (show_prob_edges) {
        gID = gfx_ID(PROB_EDGES, focus_period);
        if (all_primitives.find(gID) == all_primitives.end()) {
            // we have not drawn this one yet
            draw_edges(output.p_meshes[focus_period-1], magenta, focus_period,
                       output.p_prob_edges[focus_period-1]);
        }
    }
    if (show_degenerate_points) {
        gID = gfx_ID(DEGENERATE_POINTS, focus_period);
        if (all_primitives.find(gID) == all_primitives.end()) {
            // we have not drawn this one yet
            draw_degenerate_points(output.p_degenerate_points[focus_period-1], focus_period);
        }
    }
    if (show_index_vectors) {
        if (focus_period) {
            draw_steps(yellow, focus_period, hedgehog_scale, output.p_index_vectors[focus_period-1]);
        }
    }
    
    display();
}

poincare_map* pmap; // map needed for the construction of manifolds

void guiCallback(int control)
{
    gfx_ID gID;
    switch (control) {
        case CHANGE_MIN_PERIOD:
            break;
        case CHANGE_MAX_PERIOD:
            break;
        case CHANGE_FOCUS_PERIOD:
            return updateFocusPeriodPrimitives();
        case CHANGE_SHOW_FOCUS_MESH:
            if (!show_focus_mesh) {
                break;
            }
            gID = gfx_ID(ALL_MESH, focus_period);
            if (focus_period <= output.p_meshes.size() &&
                    all_primitives.find(gID) == all_primitives.end()) {
                // we have not drawn this one yet
                draw_mesh(output.p_meshes[focus_period-1], nvis::fvec3(0.1, 0.1, 0.1), focus_period);
                display();
            }
            break;
        case CHANGE_SHOW_CAND_MESH:
            if (!show_cand_mesh) {
                break;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_meshes.size()) ;
                    ++i) {
                gID = gfx_ID(CAND_MESH, i);
                if (all_primitives.find(gID) == all_primitives.end()) {
                    // we have not drawn this one yet
                    draw_mesh(output.p_meshes[i-1], blue, gID,
                              selection_mesh_functor<mesh_type>(output.p_cand_tris[i-1]));
                }
            }
            break;
        case CHANGE_SHOW_SING_MESH:
            if (!show_sing_mesh) {
                break;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_meshes.size()) ;
                    ++i) {
                gID = gfx_ID(SING_MESH, i);
                if (all_primitives.find(gID) == all_primitives.end()) {
                    // we have not drawn this one yet
                    draw_mesh(output.p_meshes[i-1], red, gID,
                              selection_mesh_functor<mesh_type>(output.p_sing_tris[i-1]));
                }
            }
            break;
        case CHANGE_SHOW_PROB_MESH:
            if (!show_prob_mesh) {
                return;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_meshes.size()) ;
                    ++i) {
                gID = gfx_ID(PROB_MESH, i);
                if (all_primitives.find(gID) == all_primitives.end()) {
                    // we have not drawn this one yet
                    draw_mesh(output.p_meshes[i-1], yellow, gID,
                              selection_mesh_functor<mesh_type>(output.p_prob_tris[i-1]));
                }
            }
            break;
        case CHANGE_SHOW_PROB_EDGES:
            if (!show_prob_edges) {
                return;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_meshes.size()) ;
                    ++i) {
                gID = gfx_ID(PROB_EDGES, i);
                if (all_primitives.find(gID) == all_primitives.end()) {
                    // we have not drawn this one yet
                    draw_edges(output.p_meshes[i-1], magenta, i, output.p_prob_edges[i-1]);
                }
            }
            break;
        case CHANGE_SHOW_DEGENERATE_POINTS:
            if (!show_degenerate_points) {
                return;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_degenerate_points.size()) ;
                    ++i) {
                gID = gfx_ID(DEGENERATE_POINTS, i);
                if (all_primitives.find(gID) == all_primitives.end()) {
                    // we have not drawn this one yet
                    draw_degenerate_points(output.p_degenerate_points[i-1], i);
                }
            }
            break;
        case CHANGE_SHOW_HEDGEHOG:
            if (!show_hedgehog) {
                break;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_meshes.size()) ;
                    ++i) {
                // we are inconditionally redrawing
                p_map_functor<mesh_type> pmap(i, params.dq);
                // pick a color in cyclic fashion
                nvis::fvec3 c = rainbow[1+(i%14)];
                draw_hedgehog(output.p_meshes[i-1], c, pmap, hedgehog_scale);
            }
            break;
        case CHANGE_HEDGEHOG_SCALING:
            if (show_hedgehog) {
                for (int i = std::max(1, (int)min_period) ;
                        i <= std::min((size_t)max_period, output.p_meshes.size()) ;
                        ++i) {
                    // we are inconditionally redrawing
                    p_map_functor<mesh_type> pmap(i, params.dq);
                    // pick a color in cyclic fashion
                    nvis::fvec3 c = rainbow[1+(i%14)];
                    draw_hedgehog(output.p_meshes[i-1], c, pmap, hedgehog_scale);
                }
            }
            if (show_index_vectors) {
                if (focus_period) {
                    draw_steps(yellow, focus_period, hedgehog_scale, output.p_index_vectors[focus_period-1]);
                }
            }
            break;
        case CHANGE_SHOW_PPLOT:
            break;
        case CHANGE_SHOW_INDEX_VECTORS:
            if (!show_index_vectors) {
                return;
            }
            if (focus_period) {
                draw_steps(yellow, focus_period, hedgehog_scale, output.p_index_vectors[focus_period-1]);
            }
            break;
        case CHANGE_SHOW_FIXPOINTS:
        case CHANGE_FIXPOINT_SIZE:
            if (!show_fixpoints) {
                return;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_chains.size()) ;
                    ++i) {
                gID = gfx_ID(FIXPOINTS, i);
                draw_fixpoints(output.p_chains[i-1], fp_size);
            }
            break;
        case CHANGE_SHOW_SEPARATRICES:
            if (!show_separatrices) {
                return;
            }
            for (int i = std::max(1, (int)min_period) ;
                    i <= std::min((size_t)max_period, output.p_separatrices.size()) ;
                    ++i) {
                gID = gfx_ID(SEPARATRICES, i);
                if (all_primitives.find(gID) == all_primitives.end()) {
                    draw_separatrices(output.p_separatrices[i-1], i);
                }
            }
            break;
        case CHANGE_DRAW_MANIFOLD:
        case CHANGE_MANIFOLD_LENGTH: {
            if (!draw_manifold) {
                return;
            }
            typedef newton::cached_map<poincare_map>    cached_type;
            typedef newton::rhs_map<cached_type>        rhs_type;
            
            cached_type cmap(*pmap, focus_period);
            rhs_type    rhs(cmap);
            compute_manifold(focus_period, manifold_length, rhs);
            break;
        }
        default:
            break;
    }
    
    display();
}

// GLUT code
void keyboard(unsigned char key, int x, int y)
{
    if (key == 'r') {
        resetCamera();
    } else if (key == 'x') {
        // output position coordinates
    }
    glutPostRedisplay();
}

int main(int argc, char* argv[])
{
    typedef mesh_type::data_type    data_type;
    
    // read computation parameters from command line
    init(argc, argv);
    params.hdf_name     = fin;
    params.time_step    = ts;
    params.out_name     = fout;
    params.eps          = eps;
    params.dq           = dq;
    params.mr           = max_ratio;
    params.ma           = min_area;
    params.Ma           = max_area;
    params.err          = max_err;
    params.n            = n_samples;
    params.m            = n_iter;
    params.mt           = max_tri;
    params.mp           = max_per;
    params.pm           = min_per;
    params.dx           = dx;
    params.close_d      = cl;
    
    // launch computation
    map_analysis::adaptive_map_sampling(output, params);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(params.hdf_name),
                                          params.time_step);
    field->periodic_coordinates(false);
    pmap = new poincare_map(field);
    pmap->precision(1.0e-7);
    
    bool per[2] = {true, false};
    metric_type metric(field->bounds(), per);
    
    static_data::metric = metric;
    
    // save topological information to file
    std::fstream outf(fout, std::ios::out);
    outf << output.p_chains.size() << " periods\n";
    for (int i = 0 ; i < output.p_chains.size() ; ++i) {
        outf << i + 1 << '\n';
        outf << output.p_chains[i].size() << " chains or period " << i + 1 << "\n";
        for (std::list<std::list<xavier::fixpoint> >::iterator it = output.p_chains[i].begin() ;
                it != output.p_chains[i].end() ; ++it) {
            const std::list<xavier::fixpoint>& chain = *it;
            if (it->front().saddle) {
                for (std::list<xavier::fixpoint>::const_iterator it2 = chain.begin() ;
                        it2 != chain.end() ; ++it2) {
                    outf << "saddle" << " " << it2->pos << " " << it2->evec[0] << " "
                         << it2->evec[1] << '\n';
                }
            } else {
                for (std::list<xavier::fixpoint>::const_iterator it2 = chain.begin() ;
                        it2 != chain.end() ; ++it2) {
                    outf << "center" << " " << it2->pos << '\n';
                }
            }
        }
        outf << output.p_separatrices[i].size() << " groups of separatrices for period " << i + 1 << '\n';
        for (std::list<separatrix_connection_type>::const_iterator it = output.p_separatrices[i].begin() ;
                it != output.p_separatrices[i].end() ; ++it) {
            const std::vector<separatrix_type>& seps = it->second;
            outf << seps.size() << " separatrices in group\n";
            for (int j = 0 ; j < seps.size() ; ++j) {
                outf << seps[j].size() << " points in separatrix\n";
                for (int k = 0 ; k < seps[j].size() ; ++k) {
                    outf << seps[j][k][0] << " " << seps[j][k][1] << '\n';
                }
            }
        }
    }
    outf.close();
    
    box = output.base_mesh.get_bounds();
    nvis::vec2 total_diag = box.size();
    
    // loop over orbits to determine period range
    std::vector<double> periods;
    for (int i = 0 ; i < static_data::central_map_orbits.size() ; ++i) {
        const orbit& orb = static_data::central_map_orbits[i];
        periods.push_back(orb.period());
    }
    
    // create corresponding color map
    std::vector<nvis::fvec3> scale(rainbow, &rainbow[15]);
    adaptive_color_map<double> col_map(periods, scale);
    
    // ensure that the viewport does not distort the mesh
    width = 1200;
    screen_ratio = total_diag[0] / total_diag[1];
    height = (int)((float)width / screen_ratio);
    std::cerr << "height = " << height << '\n';
    
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
    
    // initialize GUI
    init_glui();
    
    // draw poincare plot
    draw_poincare_plot(col_map);
    
    // set GLUT callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(glut_helper_reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(glut_helper_motion);
    glutKeyboardFunc(keyboard);
    
    // adjust camera to mesh
    resetCamera();
    
    // Enter GLUT event loop
    glutMainLoop();
    
    return 0;
}
























