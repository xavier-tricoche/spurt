#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>
#include <functional>

#include <boost/format.hpp>
#include <boost/limits.hpp>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

// math
#include <math/rational.hpp>

// NRRD interface
#include <image/nrrd_wrapper.hpp>

// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
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
#include "invariant_manifold.hpp"

#if _OPENMP
#include <omp.h>
#endif


using namespace xavier;

const double invalid_double = std::numeric_limits<double>::max();

int nbthreads;

// -------------------------
//
//      Data Structure
//
// -------------------------
typedef grid<double, 2>                                     plane_type;
typedef raster_data<orbit_data, double, 2>                  dataset_type;
typedef standard_map                                        map_type;

map_metric  orbit_data::metric;
int         orbit_data::max_period;

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
char*    in, *phys, *topo, *sepfile;
int     niter, width, height, maxp, minsep, maxsep, show_orb, show_vec, show_edges, cell_analysis;
int     nogfx, fast_angle, __res[2];
double  _k;
float   pt_sz, ups;
bool    logical;

std::vector<std::vector<nvis::vec2> > xavier::broken_manifolds;

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
    hestOptAdd(&hopt, "k",      "K constant",           airTypeDouble,  1,  1,  &_k,        NULL,           "constant controlling chaos");
    hestOptAdd(&hopt, "topo",   "topology data",        airTypeString,  0,  1,  &topo,      "none",         "topology file name (TXT)");
    hestOptAdd(&hopt, "sep",    "separatrix data",      airTypeString,  0,  1,  &sepfile,   "none",         "separatrix file name (base name)");
    hestOptAdd(&hopt, "n",      "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0,  1,  &maxp,      "15",           "max considered period");
    hestOptAdd(&hopt, "msep",   "min separatrixperiod",         airTypeInt,     0,  1,  &minsep,        "1",            "min considered period for separatrix construction");
    hestOptAdd(&hopt, "Msep",   "max separatrixperiod",         airTypeInt,     0,  1,  &maxsep,        "15",           "max considered period for separatrix construction");
    hestOptAdd(&hopt, "r",      "resolution",           airTypeInt,     2,  2,  &__res,     NULL,           "plane resolution");
    hestOptAdd(&hopt, "ps",     "point size",           airTypeFloat,   0,  1,  &pt_sz,     "2",            "point size for display");
    hestOptAdd(&hopt, "cell",   "do cell analysis",     airTypeInt,     0,  0,  &cell_analysis, NULL,       "carry out cell analysis");
    hestOptAdd(&hopt, "orb",    "show orbits",          airTypeInt,     0,  0,  &show_orb,  NULL,           "show automatically generated orbits");
    hestOptAdd(&hopt, "fast",   "fast angle",           airTypeInt,     0,  0,  &fast_angle,NULL,           "skip most tricky cases in angle computation");
    hestOptAdd(&hopt, "nogfx",  "no graphics",          airTypeInt,     0,  0,  &nogfx,     NULL,           "no display");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Safety factor probe",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline bool equal(bool a, bool b)
{
    return ((a && b) || (!a && !b));
}

template<typename InputIterator>
InputIterator epsilon_find(const InputIterator first, const InputIterator last,
                           const nvis::vec2& x, const map_metric& metric, double eps)
{
    // return element corresponding to min distance under epsilon to
    // reference position
    std::map<double, InputIterator> norm_to_it;
    for (InputIterator i=first ; i!=last ; ++i) {
        const nvis::vec2& y = i->pos;
        norm_to_it[metric.distance(x, y)] = i;
    }
    if (norm_to_it.begin()->first < eps) {
        return norm_to_it.begin()->second;
    }
    
    return last;
}

// -------------------------
//
//          Analysis
//
// -------------------------
typedef rational_surface_found::edge_type               segment_type;
typedef rational_surface_found::edge_point              edge_point;
typedef boost::rational<int>                            rational_type;
typedef edge<nvis::ivec2, nvis::lexicographical_order>  edge_type;
typedef std::vector<nvis::vec2>                         orbit_type;
typedef std::pair<orbit_type, nvis::fvec3>              color_orbit_type;

template<typename T, typename Compare = std::less<T> >
struct tagged {
    typedef tagged<T, Compare>  self_type;
    
    tagged() : t(), p(0) {}
    tagged(T _t, int _p) : t(_t), p(_p) {}
    
    int& tag() {
        return p;
    }
    int tag() const {
        return p;
    }
    T& object() {
        return t;
    }
    const T& object() const {
        return t;
    }
    
    bool operator<(const self_type& hp) const {
        if (p < hp.p) {
            return true;
        } else if (p > hp.p) {
            return false;
        }
        Compare Lt;
        return Lt(t, hp.t);
    }
    
    T   t;
    int p;
};
typedef tagged<edge_type>                                   p_edge_type;
typedef tagged<nvis::vec2, nvis::lexicographical_order>     p_cell_type;

struct Lt_bbox {
    nvis::lexicographical_order Lt;
    
    bool operator()(const nvis::bbox2& b0, const nvis::bbox2& b1) const {
        if (Lt(b0.min(), b1.min())) {
            return true;
        } else if (Lt(b1.min(), b0.min())) {
            return false;
        }
        return Lt(b0.max(), b1.max());
    }
};

typedef tagged<nvis::bbox2, Lt_bbox>            p_box_type;

std::map<p_edge_type, double>                   _edge_angles;
std::set<p_edge_type>                           _unresolved_edges;
std::vector<tagged<p_cell_type> >               _cells_indices;
std::vector<tagged<p_box_type> >                _boxes_indices;
std::vector<p_edge_type>                        _failed_edges;
std::vector<std::vector<fixpoint> >             _chains;
std::vector<std::vector<separatrix> >           _separatrices;

// -------------------------
//
//          Display
//
// -------------------------
std::vector<color_orbit_type>                   _orbits;
xavier::discrete_color_map<int>*                _cmap;
xavier::map_analysis_param                      _params, _debug_params;
std::vector<nvis::ivec2>                        _saddle_cells, _center_cells;
nvis::vec2                                      _last;
std::vector<nvis::vec2>                         _problematic_seeds;
std::vector<color_orbit_type>                   _rational_surfaces;
std::vector<nvis::vec2>                         _saddles;
std::vector<nvis::vec2>                         _centers;

struct rational_segment {
    edge_point pt[2];
    rational_type sf;
};

std::vector<rational_segment>                   _rational_segments;
std::vector<fixpoint>                           _fixpoints;
std::vector<int>                                _tagged_cells;
double                                          _index_dx;
std::vector<rational_type>                      _valid_rationals;

enum {
    FIXED_POINT = 0,
    CONVERGED,
    DEGENERATE
};

static void init()
{
    nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    _plane_bounds = nvis::bbox2(nvis::vec2(0, 0),
                                nvis::vec2(1, 1));
                                
    int __nx = __res[0];
    int __ny = __res[1];
    _plane_res = nvis::ivec2(__nx, __ny);
    _plane_metric.bounds() = _plane_bounds;
    _plane_metric.periodic(0) = true;
    _plane_metric.periodic(1) = true;
    int npoints = _plane_res[0] * _plane_res[1];
    
    std::cerr << "plane resolution = " << _plane_res << std::endl;
    
    orbit_data::metric = _plane_metric;
    orbit_data::max_period = maxp;
    
    std::set<rational_type> tmp_r;
    for (int n=1 ; n<=maxp ; ++n) {
        for (int d=1 ; d<=n ; ++d) {
            tmp_r.insert(rational_type(n,d));
        }
    }
    std::copy(tmp_r.begin(), tmp_r.end(), std::back_inserter(_valid_rationals));
    
    GLUT_helper::box = _plane_bounds;
    
    _plane = new plane_type(_plane_res, _plane_bounds);
    _dataset = new dataset_type(*_plane, orbit_data());
    _params.nb_iterations = niter;
    _params.max_period = maxp;
    _params.metric = _plane_metric;
    _plane_spacing = _plane->spacing();
    
    _debug_params.verbose = true;
    _debug_params.metric = _plane_metric;
    
    std::cerr << "plane spacing = " << _plane_spacing << std::endl;
}

// per cell period range
void compute_cell_periods(std::vector<p_cell_type>& cells)
{
    nvis::ivec2 cellres(_plane_res[0]-1, _plane_res[1]-1);
    int ncells = cellres[0]*cellres[1];
    
    std::vector<p_cell_type> cached_cells[nbthreads];
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<ncells ; ++n) {
            int i = n % cellres[0];
            int j = n / cellres[0];
            nvis::ivec2 cell_id(i,j);
            
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            std::vector<int> periods;
            best_periods(periods, *_dataset, cell_id, _plane_metric, maxp, 5);
            for (int k=0 ; k<periods.size() ; ++k) {
                if (periods[k]<1 || periods[k]>maxp) {
                    std::cerr << "wrong period (" << periods[k] << ") in cell #" << cell_id << std::endl;
                }
                cached_cells[thread_id].push_back(p_cell_type(cell_id, periods[k]));
            }
            
            if (thread_id == 0) {
                std::cerr << "\rcomputed period of " << n << " / " << ncells
                          << " (" << 100*n/ncells << "%)          \r" << std::flush;
            }
        }
    }
    std::cerr << "\nperiod range computation took " << _timer.elapsed() << '\n';
    
    for (int i=0 ; i<nbthreads ; ++i) {
        std::copy(cached_cells[i].begin(), cached_cells[i].end(),
                  std::back_inserter(cells));
    }
}

void __edge_rotation(std::vector<p_edge_type>& failed,
                     const std::vector<p_edge_type>& edges, double lmin,
                     int nb_sub=4)
{

    typedef std::pair<p_edge_type, double> pair_type;
    std::vector<pair_type> cached_angles[nbthreads];
    std::vector<p_edge_type> cached_edges[nbthreads];
    
    standard_map pmap(_k);
    
    nvis::timer _timer;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<edges.size() ; ++n) {
            int thread_id = 0;
            _debug_params.verbose = false;
            std::ostringstream os;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            
            if (thread_id == 0) {
                std::cerr << "\rcomputed " << n << " edge angles / "
                          << edges.size() << " (" << 100*n/edges.size()
                          << "%)          \r" << std::flush;
            }
            
            const p_edge_type& e = edges[n];
            nvis::ivec2 i0 = e.object()[0];
            nvis::ivec2 i1 = e.object()[1];
            // if (nvis::all(i0 == nvis::ivec2(86, 69)) &&
            //  nvis::all(i1 == nvis::ivec2(86, 70))) {
            //  _debug_params.verbose = true;
            //  os << "processing edge " << i0 << " - " << i1 << std::endl;
            // }
            
            
            if ((*_dataset)(i0).sf == 0 || (*_dataset)(i1).sf == 0) {
                continue;
            }
            int p = e.tag();
            nvis::vec2 v0 = (*_dataset)(i0).vector(p);
            nvis::vec2 v1 = (*_dataset)(i1).vector(p);
            if (_debug_params.verbose) {
                os << "v0 = " << v0 << ", v1 = " << v1 << std::endl;
            }
            
            map_type* amap = pmap.clone();
            
            double theta = 0;
            if (lmin > 0) {
                nvis::vec2 x0 = (*_plane)(i0);
                nvis::vec2 x1 = (*_plane)(i1);
                double du = 1./(double)(nb_sub+1);
                nvis::vec2 x[nb_sub+2], v[nb_sub+2];
                double u = du;
                try {
                    for (int i=1; i<=nb_sub ; ++i, u+=du) {
                        x[i] = (1.-u)*x0 + u*x1;
                        v[i] = rhs(x[i], pmap, p, _plane_metric);
                    }
                } catch(...) {
                    continue;
                }
                x[0] = x0;
                x[nb_sub+1] = x1;
                v[0] = v0;
                v[nb_sub+1] = v1;
                try {
                    for (int i=0 ; i<=nb_sub ; ++i) {
                        theta += adaptive_rotation_angle(x[i], v[i], x[i+1], v[i+1], *amap, p, lmin, _debug_params);
                    }
                    if (_debug_params.verbose) {
                        os << "total dtheta for edge = " << theta << std::endl;
                    }
                    cached_angles[thread_id].push_back(pair_type(e, theta));
                } catch (index_step_size_underflow& err) {
                    if (_debug_params.verbose) {
                        os << "step size underflow caught" << std::endl;
                    }
                    cached_edges[thread_id].push_back(e);
                } catch(...) {
                }
            } else {
                theta = signed_angle(v0, v1);
                if (_debug_params.verbose) {
                    os << "direct angle = " << theta << std::endl;
                }
                if (fabs(theta) < MAX_ANGLE_VARIATION) {
                    cached_angles[thread_id].push_back(std::pair<p_edge_type, double>(edges[n], theta));
                    if (_debug_params.verbose) {
                        os << "valid angle. done" << std::endl;
                    }
                } else {
                    if (_debug_params.verbose) {
                        os << "angle is too large." << std::endl;
                    }
                    cached_edges[thread_id].push_back(edges[n]);
                }
            }
            if (_debug_params.verbose) {
                std::cerr << os.str();
            }
        }
    }
    std::cerr << "processing of " << edges.size() << " edges took "
              << _timer.elapsed() << " s.                          \n";
              
    failed.clear();
    for (int i=0 ; i<nbthreads ; ++i) {
        for (int j=0 ; j<cached_angles[i].size() ; ++j) {
            const pair_type& p = cached_angles[i][j];
            _edge_angles[p.first] = p.second;
        }
        std::copy(cached_edges[i].begin(), cached_edges[i].end(),
                  std::back_inserter(failed));
    }
    std::cerr << "We were unable to process " << failed.size() << " edges out of "
              << edges.size() << " input edges ("
              << 100*failed.size()/edges.size() << "%)\n";
}

// per edge vector rotation
void compute_edge_rotation(const std::vector<p_cell_type>& cells)
{
    // assume cells sorted by period
    
    std::cerr << "there are " << cells.size() << " cells in input\n";
    
    std::set<p_edge_type> unique_edges;
    std::cerr << "identifying edges... " << std::flush;
    edge_type e;
    for (int n=0 ; n<cells.size() ; ++n) {
        nvis::ivec2 id = cells[n].object();
        int p = cells[n].tag();
        e = edge_type(id, id + nvis::ivec2(1,0));
        unique_edges.insert(p_edge_type(e, p));
        e = edge_type(id + nvis::ivec2(1,0), id + nvis::ivec2(1,1));
        unique_edges.insert(p_edge_type(e, p));
        e = edge_type(id + nvis::ivec2(1,1), id + nvis::ivec2(0,1));
        unique_edges.insert(p_edge_type(e, p));
        e = edge_type(id + nvis::ivec2(0,1), id);
        unique_edges.insert(p_edge_type(e, p));
    }
    std::cerr << ": " << unique_edges.size() << " edges found\n";
    
    for (std::set<p_edge_type>::const_iterator it=unique_edges.begin() ;
            it!=unique_edges.end() ; ++it) {
        // initial value is invalid angle
        _edge_angles[*it] = invalid_double;
    }
    std::vector<p_edge_type> all_edges(unique_edges.begin(), unique_edges.end());
    std::cerr << "edge map created\n";
    
    std::vector<p_edge_type> difficult_edges, problematic_edges;
    // easy
    // __edge_rotation(difficult_edges, all_edges, 0);
    // difficult
    double lmin = 0.99*std::min(_plane_spacing[0], _plane_spacing[1])/64.;
    __edge_rotation(problematic_edges, all_edges, lmin);
    // problematic
    if (!fast_angle) {
        lmin = 0.99*std::min(_plane_spacing[0], _plane_spacing[1])/512;
        __edge_rotation(difficult_edges, problematic_edges, lmin);
    }
    _unresolved_edges.insert(difficult_edges.begin(), difficult_edges.end());
}

inline double edge_angle(const nvis::ivec2& cell_id, int edge_id, int p)
{

    const nvis::ivec2 _coords[] = {
        nvis::ivec2(0,0),
        nvis::ivec2(1,0),
        nvis::ivec2(1,1),
        nvis::ivec2(0,1),
        nvis::ivec2(0,0)
    };
    
    nvis::ivec2 i0 = cell_id + _coords[edge_id];
    nvis::ivec2 i1 = cell_id + _coords[edge_id+1];
    
    p_edge_type e(edge_type(i0, i1), p);
    if (_unresolved_edges.find(e) != _unresolved_edges.end()) {
        throw std::runtime_error("no valid value available for this edge");
    }
    std::map<p_edge_type, double>::iterator it = _edge_angles.find(e);
    if (it == _edge_angles.end()) {
        throw std::runtime_error("no value available for this edge");
    } else {
        return (edge_id > 1 ? -1. : 1.)*it->second;
    }
}

int __box_index_1D(const std::vector<nvis::ivec2>& cell_ids, int p, bool verbose=false)
{
    double theta = 0;
    // first
    theta += edge_angle(cell_ids.front(), 0, p);
    theta += edge_angle(cell_ids.front(), 2, p);
    theta += edge_angle(cell_ids.front(), 3, p);
    // last
    theta += edge_angle(cell_ids.back(), 0, p);
    theta += edge_angle(cell_ids.back(), 1, p);
    theta += edge_angle(cell_ids.back(), 2, p);
    // intermediate
    for (int i=1 ; i<cell_ids.size()-1 ; ++i) {
        theta += edge_angle(cell_ids[i], 0, p);
        theta += edge_angle(cell_ids[i], 2, p);
    }
    theta /= 2.*M_PI;
    return lround(theta);
}

int __box_index_2D(const std::vector<std::vector<nvis::ivec2> >& cell_ids, int p, bool verbose=false)
{
    double theta = 0;
    // left  and right borders
    for (int row=0 ; row<cell_ids.size() ; ++row) {
        const nvis::ivec2& left = cell_ids[row].front();
        theta += edge_angle(left, 3, p);
        const nvis::ivec2& right = cell_ids[row].back();
        theta += edge_angle(right, 1, p);
    }
    // top and bottom borders
    for (int col=0 ; col<cell_ids[0].size() ; ++col) {
        const nvis::ivec2& bottom = cell_ids[0][col];
        theta += edge_angle(bottom, 0, p);
        const nvis::ivec2& top = cell_ids.back()[col];
        theta += edge_angle(top, 2, p);
    }
    theta /= 2.*M_PI;
    return lround(theta);
}

void __multi_cell_index(const nvis::ivec2& cell_id, int p, bool verbose=false)
{
    int missed=0;
    for (int i=0 ; i<4 ; ++i) {
        try {
            double theta = edge_angle(cell_id, i, p);
        } catch(...) {
            missed += (1 << i);
        }
    }
    
    try {
        int idx;
        nvis::bbox2 bounds;
        switch (missed) {
            case 1: {
                if (cell_id[1] == 0) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id - nvis::ivec2(0, 1));
                box[1].push_back(cell_id);
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(0,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 2: {
                if (cell_id[0] == _plane_res[0]-1) {
                    return;
                }
                std::vector<nvis::ivec2> box(2);
                box[0] = cell_id;
                box[1] = cell_id + nvis::ivec2(1,0);
                bounds.min() = (*_plane)(cell_id);
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,1));
                idx = __box_index_1D(box, p, verbose);
                break;
            }
            case 3: {
                if (cell_id[0] == _plane_res[0]-1 || cell_id[1] == 0) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id - nvis::ivec2(0,1));
                box[0].push_back(cell_id + nvis::ivec2(-1,1));
                box[1].push_back(cell_id);
                box[1].push_back(cell_id + nvis::ivec2(1,0));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(0,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,1));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 4: {
                if (cell_id[1] == _plane_res[1]-1) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id);
                box[1].push_back(cell_id + nvis::ivec2(0,1));
                bounds.min() = (*_plane)(cell_id);
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 5: {
                if (cell_id[1] == _plane_res[1]-1 || cell_id[1] == 0) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(3);
                box[0].push_back(cell_id - nvis::ivec2(0,1));
                box[1].push_back(cell_id);
                box[2].push_back(cell_id + nvis::ivec2(0,1));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(0,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 6: {
                if (cell_id[0] == _plane_res[0]-1 || cell_id[1] == _plane_res[1]-1) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id);
                box[0].push_back(cell_id + nvis::ivec2(1,0));
                box[1].push_back(cell_id + nvis::ivec2(0,1));
                box[1].push_back(cell_id + nvis::ivec2(1,1));
                bounds.min() = (*_plane)(cell_id);
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 7: {
                if (cell_id[0] == _plane_res[0]-1 || cell_id[1] == _plane_res[1]-1 ||
                        cell_id[1] == 0) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(3);
                box[0].push_back(cell_id + nvis::ivec2(-1,0));
                box[0].push_back(cell_id + nvis::ivec2(-1,1));
                box[1].push_back(cell_id);
                box[1].push_back(cell_id + nvis::ivec2(1,0));
                box[2].push_back(cell_id + nvis::ivec2(0,1));
                box[2].push_back(cell_id + nvis::ivec2(1,1));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(0,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 8: {
                if (cell_id[0] == 0) {
                    return;
                }
                std::vector<nvis::ivec2> box(2);
                box[0] = cell_id - nvis::ivec2(1,0);
                box[1] = cell_id;
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,0));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
                idx = __box_index_1D(box, p, verbose);
                break;
            }
            case 9: {
                if (cell_id[0] == 0 || cell_id[1] == 0) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id - nvis::ivec2(1,1));
                box[0].push_back(cell_id - nvis::ivec2(0,1));
                box[1].push_back(cell_id - nvis::ivec2(1,0));
                box[1].push_back(cell_id);
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 10: {
                if (cell_id[0] == 0 || cell_id[0] == _plane_res[0]-1) {
                    return;
                }
                std::vector<nvis::ivec2> box(3);
                box[0] = cell_id - nvis::ivec2(1,0);
                box[1] = cell_id;
                box[2] = cell_id + nvis::ivec2(1,0);
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,0));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,1));
                idx = __box_index_1D(box, p, verbose);
                break;
            }
            case 11: {
                if (cell_id[0] == 0 || cell_id[0] == _plane_res[0]-1 ||
                        cell_id[1] == 0) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id - nvis::ivec2(1,1));
                box[0].push_back(cell_id - nvis::ivec2(0,1));
                box[0].push_back(cell_id + nvis::ivec2(1,-1));
                box[1].push_back(cell_id - nvis::ivec2(1,0));
                box[1].push_back(cell_id);
                box[1].push_back(cell_id + nvis::ivec2(1,0));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 12: {
                if (cell_id[0] == 0 || cell_id[1] == _plane_res[1]-1) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id - nvis::ivec2(1,0));
                box[0].push_back(cell_id);
                box[1].push_back(cell_id + nvis::ivec2(-1,1));
                box[1].push_back(cell_id + nvis::ivec2(0,1));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,0));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 13: {
                if (cell_id[0] == 0 || cell_id[1] == 0 ||
                        cell_id[1] == _plane_res[1]-1) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(3);
                box[0].push_back(cell_id + nvis::ivec2(-1,-1));
                box[0].push_back(cell_id + nvis::ivec2(0,-1));
                box[1].push_back(cell_id + nvis::ivec2(-1,0));
                box[1].push_back(cell_id);
                box[2].push_back(cell_id + nvis::ivec2(-1,1));
                box[2].push_back(cell_id + nvis::ivec2(0,1));
                bounds.min() = (*_plane)(cell_id + nvis::ivec2(-1,-1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 14: {
                if (cell_id[0] == 0 || cell_id[0] == _plane_res[0]-1 ||
                        cell_id[1] == _plane_res[1]-1) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(2);
                box[0].push_back(cell_id - nvis::ivec2(1,0));
                box[0].push_back(cell_id);
                box[0].push_back(cell_id + nvis::ivec2(1,0));
                box[1].push_back(cell_id + nvis::ivec2(-1,1));
                box[1].push_back(cell_id + nvis::ivec2(0,1));
                box[1].push_back(cell_id + nvis::ivec2(1,1));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,0));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            case 15: {
                if (cell_id[0] == 0 || cell_id[0] == _plane_res[0]-1 ||
                        cell_id[1] == 0 || cell_id[1] == _plane_res[1]-1) {
                    return;
                }
                std::vector<std::vector<nvis::ivec2> > box(3);
                box[0].push_back(cell_id + nvis::ivec2(-1,-1));
                box[0].push_back(cell_id + nvis::ivec2(0,-1));
                box[0].push_back(cell_id + nvis::ivec2(1,-1));
                box[1].push_back(cell_id + nvis::ivec2(-1,0));
                box[1].push_back(cell_id);
                box[1].push_back(cell_id + nvis::ivec2(1,0));
                box[2].push_back(cell_id + nvis::ivec2(-1,1));
                box[2].push_back(cell_id + nvis::ivec2(0,1));
                box[2].push_back(cell_id + nvis::ivec2(1,1));
                bounds.min() = (*_plane)(cell_id - nvis::ivec2(1,1));
                bounds.max() = (*_plane)(cell_id + nvis::ivec2(2,2));
                idx = __box_index_2D(box, p, verbose);
                break;
            }
            default: {
                std::ostringstream os;
                os << "invalid code (" << missed << ")" << std::endl;
                std::cerr << os.str();
                return;
            }
        }
        if (idx != 0) {
            std::ostringstream os;
            os << "found index " << idx << " in " << p << "-multicell" << std::endl;
            std::cerr << os.str();
            _boxes_indices.push_back(tagged<p_box_type>(p_box_type(bounds, p), idx));
        }
        // else {
        //  std::ostringstream os;
        //  os << "index zero found in multicell" << std::endl;
        //  std::cerr << os.str();
        // }
    } catch(...) {
        // std::ostringstream os;
        // os << "index computation in multicell threw an exception" << std::endl;
        // std::cerr << os.str();
    }
}

int __index(const nvis::ivec2& cell_id, int p, bool verbose=false)
{
    nvis::ivec2 pt[5];
    pt[0] = cell_id;
    pt[1] = cell_id + nvis::ivec2(1,0);
    pt[2] = cell_id + nvis::ivec2(1,1);
    pt[3] = cell_id + nvis::ivec2(0,1);
    pt[4] = cell_id;
    
    double theta = 0;
    bool valid = true;
    std::ostringstream os;
    for (int i=0 ; i<4 ; ++i) {
        edge_type e(pt[i], pt[i+1]);
        p_edge_type et(e, p);
        if (_unresolved_edges.find(et) != _unresolved_edges.end()) {
            __multi_cell_index(cell_id, p, verbose);
            return 0;
        }
        std::map<p_edge_type, double>::iterator it = _edge_angles.find(et);
        if (it == _edge_angles.end()) {
            // should not happen: all edges are initialized with invalid value
            if (verbose) {
                os << "edge " << pt[i] << "-" << pt[i+1]
                   << " was not found in database\n";
            }
            valid = false;
            break;
        }
        double dtheta = it->second;
        if (dtheta == invalid_double) {
            if (verbose) {
                os << "edge " << pt[i] << "-" << pt[i+1]
                   << " has an invalid angle\n";
            }
            valid = false;
            break;
        }
        if (verbose) {
            os << "edge " << pt[i] << "-" << pt[i+1]  << " has angle " << dtheta << '\n';
        }
        if (i>=2) {
            if (verbose) {
                os << "edge " << pt[i] << "-" << pt[i+1]
                   << " was computed the other way around\n";
            }
            dtheta *= -1;
        }
        theta += dtheta;
    }
    if (verbose) {
        os << "total angle is " << theta << std::endl;
        std::cerr << os.str();
    }
    
    if (valid) {
        theta /= 2.*M_PI;
        return lround(theta);
    } else {
        return 0;
    }
}

void compute_cell_index(const std::vector<p_cell_type>& cells)
{
    std::vector<tagged<p_cell_type> > cached_indices[nbthreads];
    
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<cells.size() ; ++n) {
            int thread_id = 0;
            bool verbose = false;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            nvis::ivec2 cell_id = cells[n].object();
            int p = cells[n].tag();
            int index = __index(cell_id, p);
            if (nvis::all(cell_id == nvis::ivec2(49, 68))) {
                std::ostringstream os;
                os << std::endl << std::endl << "index in cell " << cell_id << " is " << index << std::endl << std::endl;
                std::cerr << os.str();
            }
            if (index != 0) {
                cached_indices[thread_id].push_back(tagged<p_cell_type>(cells[n], index));
            }
            
            if (thread_id == 0) {
                std::cerr << "\rcomputed the index of " << n << " cells / "
                          << cells.size() << " (" << 100*n/cells.size() << "%)         \r"
                          << std::flush;
            }
        }
    }
    std::cerr << "\ncomputation of cell indices took " << _timer.elapsed() << " s.                          \n";
    
    int counter = 0;
    for (int i=0 ; i<nbthreads ; ++i) {
        std::copy(cached_indices[i].begin(), cached_indices[i].end(),
                  std::back_inserter(_cells_indices));
    }
    counter = _cells_indices.size();
    std::cerr << "There are " << counter << " cells containing fixed points "
              << "(" << 100*counter/cells.size() << "%)\n";
    std::cerr << "these cells are: " << std::endl;
    for (int i=0 ; i<_cells_indices.size() ; ++i) {
        std::cerr << "id = " << _cells_indices[i].object().object() << ", period = " << _cells_indices[i].object().tag()
                  << ", index = " << _cells_indices[i].tag() << std::endl;
    }
    std::vector<int> nbsaddles(maxp), nbcenters(maxp);
    for (int i=0 ; i<_cells_indices.size() ; ++i) {
        int idx = _cells_indices[i].tag();
        int p = _cells_indices[i].object().tag();
        if (idx < 0) {
            ++nbsaddles[p-1];
            nvis::ivec2 id = _cells_indices[i].object().object();
            nvis::vec2 x = (*_plane)(id) + 0.5*_plane_spacing;
            _saddles.push_back(x);
        } else if (idx > 0) {
            ++nbcenters[p-1];
            nvis::ivec2 id = _cells_indices[i].object().object();
            nvis::vec2 x = (*_plane)(id) + 0.5*_plane_spacing;
            _centers.push_back(x);
        }
    }
    std::cerr << "statistics by period:\n";
    for (int i=0 ; i<maxp ; ++i) {
        std::cerr << "period " << i+1 << ": " << nbsaddles[i] << " saddles, "
                  << nbcenters[i] << " centers\n";
    }
}

double mindistance(const nvis::vec2& x, const std::vector<fixpoint>& chain)
{
    std::vector<double> dist(chain.size());
    for (int i=0 ; i<chain.size() ; ++i) {
        dist[i]= _plane_metric.distance(x, chain[i].pos);
    }
    return *std::min_element(dist.begin(), dist.end());
}

double chain_distance(const std::vector<fixpoint>& chain0,
                      const std::vector<fixpoint>& chain1)
{
    if (chain0.size() != chain1.size() )  {
        return invalid_double;
    }
    if (chain0.front().saddle && !chain1.front().saddle ||
            !chain0.front().saddle && chain1.front().saddle) {
        return invalid_double;
    }
    return mindistance(chain0.front().pos, chain1);
}

void link_chains(const std::vector<fixpoint>& fps)
{
    if (!fps.size()) {
        return;
    }
    int period = fps.front().K;
    bool saddle = fps.front().saddle;
    std::vector<bool> removed(fps.size(), false);
    
    standard_map pmap(_k);
    
    // sort fixed points by norm of the associated p map
    std::map<double, int> norm_to_id;
    for (int i=0 ; i<fps.size() ; ++i) {
        try {
            nvis::vec2 y = pmap.map(fps[i].pos, period);
            norm_to_id[_plane_metric.distance(fps[i].pos, y)] = i;
        } catch(...) { // discard faulty positions
            removed[i] = true;
        }
    }
    
    // std::cerr << "filtering " << (saddle ? "saddles" : "centers") << " of period " << period << std::endl;
    // std::cerr << "there are " << norm_to_id.size() << " valid fixed points in input\n";
    
    for (std::map<double, int>::const_iterator it=norm_to_id.begin() ; it!=norm_to_id.end() ; ++it) {
        int n = it->second;
        if (removed[n]) {
            continue;
        }
        removed[n] = true;
        _chains.push_back(std::vector<fixpoint>());
        std::vector<fixpoint>& chain = _chains.back();
        chain.push_back(fps[n]);
        bool failed = false;
        nvis::vec2 x = fps[n].pos;
        // std::cerr << "starting at " << x << std::endl;
        for (int j=1 ; j<period ; ++j) {
            try {
                x = pmap.map(x, 1);
                // std::cerr << "next point is at " << x << std::endl;
                std::vector<fixpoint>::const_iterator iter =
                    epsilon_find(fps.begin(), fps.end(), x, _plane_metric,
                                 0.25*nvis::norm(_plane_spacing));
                if (iter == fps.end()) {
                    // std::cerr << "no fixed point found there\n";
                    failed = true;
                    break;
                } else {
                    int k = std::distance(fps.begin(), iter);
                    // std::cerr << "found fixed point has index " << k << std::endl;
                    if (removed[k]) {
                        // std::cerr << "this index was already processed. failed\n";
                        failed = true;
                        break;
                    } else if (!equal(fps[k].saddle, saddle)) {
                        // std::cerr << "inconsistent types. failed\n";
                        failed = true;
                        break;
                    } else {
                        // std::cerr << "ok\n";
                        chain.push_back(fps[k]);
                        removed[k] = true;
                    }
                }
            } catch(...) {
                failed = true;
                break;
            }
        }
        if (failed) {
            _chains.pop_back();
        }
    }
}

void __fixed_points(std::vector<fixpoint>& fixpoints,
                    std::map<p_cell_type, nvis::vec2>& cell_to_seed,
                    std::map<p_cell_type, nvis::vec2>& completed_cells,
                    bool prefound = false)
{

    standard_map pmap(_k);
    
    std::map<p_cell_type, nvis::vec2> detected_cells;
    
    std::vector<p_cell_type> cells;
    typedef std::map<p_cell_type, nvis::vec2>::iterator iterator;
    for (iterator it=cell_to_seed.begin() ; it!=cell_to_seed.end() ; ++it) {
        cells.push_back(it->first);
    }
    
    std::cerr << "processing " << cells.size() << " cells in input\n";
    
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<cells.size() ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            
            const p_cell_type& pcell = cells[n];
            nvis::ivec2 cell_id = pcell.object();
            int period = pcell.tag();
            
            std::cerr << "looking for fixed point in cell " << cell_id << std::endl;
            // std::cerr << "index make-up for this cell is: " << std::endl;
            // int idx = __index(cell_id, period, true);
            
            nvis::vec2 x = cell_to_seed[pcell];
            std::vector<nvis::vec2> chain;
            map_type* amap = pmap.clone();
            
            fixpoint fp;
            nvis::bbox2 bounds;
            bounds.min() = (*_plane)(cell_id);
            bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
            bool found = meta_newton_stdmap<map_type>(*amap, _plane_metric, bounds, x, 5, period,
                         fp, chain, 1.0e-5, true, 50, prefound);
                         
            if (found) {
#pragma openmp atomic
                completed_cells[pcell] = fp.pos;
#pragma openmp atomic
                fixpoints.push_back(fp);
                for (int i=1 ; i<chain.size() ; ++i) {
                    p_cell_type id(_plane->local_coordinates(chain[i]).first, period);
                    iterator it = cell_to_seed.find(id);
                    if (it!=cell_to_seed.end())
#pragma openmp atomic
                        it->second = chain[i]; // this cell was already identified
                        
#pragma openmp atomic
                    detected_cells[id] = chain[i]; // this cell is new
                }
                std::ostringstream os;
                os << fp << " found" << std::endl;
                std::cerr << os.str();
            }
            
            if (thread_id == 0) {
#pragma openmp atomic
                int nb_found = fixpoints.size();
                std::cerr << "processed " << n << " candidate cells from " << cells.size()
                          << " (" << 100*n/cells.size() << "%), found " << nb_found << " fixed points\n";
            }
        }
    }
    std::cerr << "this round of fixpoint extraction took " << _timer.elapsed() << '\n';
    
    cell_to_seed.clear();
    int nb_found_already = 0;
    for (iterator it=detected_cells.begin() ; it != detected_cells.end() ; ++it) {
        iterator jt = completed_cells.find(it->first);
        if (jt != completed_cells.end()) {
            ++nb_found_already;
        } else {
            cell_to_seed[it->first] = it->second;
        }
    }
    
    std::cerr << completed_cells.size() << " fixed points were found\n";
    std::cerr << cell_to_seed.size() << " additional fixed points were discovered\n";
}

void __fixed_points_from_boxes(std::vector<fixpoint>& fixpoints,
                               std::map<p_cell_type, nvis::vec2>& cell_to_seed,
                               std::map<p_cell_type, nvis::vec2>& completed_cells)
{

    standard_map pmap(_k);
    
    std::map<p_cell_type, nvis::vec2> detected_cells;
    
    std::cerr << "processing " << _boxes_indices.size() << " boxes in input\n";
    
    nvis::timer _timer;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<_boxes_indices.size() ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            
            const p_box_type b = _boxes_indices[n].object();
            int period = b.tag();
            const nvis::bbox2& bounds = b.object();
            
            nvis::vec2 x = bounds.center();
            std::vector<nvis::vec2> chain;
            map_type* amap = pmap.clone();
            
            fixpoint fp;
            bool found = meta_newton(*amap, _plane_metric, bounds, x, 6, period,
                                     fp, chain, 1.0e-5, 0, false, 50, false);
                                     
            if (found) {
                nvis::ivec2 cell_id = _plane->local_coordinates(chain[0]).first;
                p_cell_type pcell(cell_id, period);
                
#pragma openmp atomic
                completed_cells[pcell] = fp.pos;
#pragma openmp atomic
                fixpoints.push_back(fp);
                for (int i=1 ; i<chain.size() ; ++i) {
                    p_cell_type id(_plane->local_coordinates(chain[i]).first, period);
                    std::map<p_cell_type, nvis::vec2>::iterator it = cell_to_seed.find(id);
                    if (it!=cell_to_seed.end())
#pragma openmp atomic
                        it->second = chain[i]; // this cell was already identified
                        
#pragma openmp atomic
                    detected_cells[id] = chain[i]; // this cell is new
                }
                std::ostringstream os;
                os << fp << " found" << std::endl;
                std::cerr << os.str();
            }
            
//             if (thread_id == 0) {
// #pragma openmp atomic
//                 int nb_found = fixpoints.size();
//                 std::cerr << "processed " << n << " candidate cells from " << _boxes_indices.size()
//                           << " (" << 100*n/_boxes_indices.size() << "%), found " << nb_found << " fixed points\n";
//             }
        }
    }
    std::cerr << "this round of fixpoint extraction took " << _timer.elapsed() << '\n';
    
    int nb_found_already = 0;
    for (std::map<p_cell_type, nvis::vec2>::iterator it=detected_cells.begin() ; it != detected_cells.end() ; ++it) {
        std::map<p_cell_type, nvis::vec2>::iterator jt = completed_cells.find(it->first);
        if (jt != completed_cells.end()) {
            ++nb_found_already;
        } else {
            cell_to_seed[it->first] = it->second;
        }
    }
    
    std::cerr << completed_cells.size() << " fixed points were found\n";
    std::cerr << cell_to_seed.size() << " additional fixed points were discovered\n";
}

void compute_fixed_points()
{
    std::sort(_cells_indices.begin(), _cells_indices.end());
    
    standard_map pmap(_k);
    
    std::map<p_cell_type, nvis::vec2>   cell_to_seed;
    for (int i=0 ; i<_cells_indices.size() ; ++i) {
        const nvis::ivec2& cell_id = _cells_indices[i].object().object();
        cell_to_seed[_cells_indices[i].object()] = (*_plane)(cell_id) + 0.5*_plane_spacing;
    }
    
    std::map<p_cell_type, nvis::vec2>                   completed_cells;
    std::vector<fixpoint>&                              fixpoints = _fixpoints;
    
    for (int i=0 ; cell_to_seed.size() ; ++i) {
        std::cerr << "\n\n\n\nSTARTING ROUND " << i << " of fixed point extraction\n";
        __fixed_points(fixpoints, cell_to_seed, completed_cells, (i>0));
        if (!i) {
            __fixed_points_from_boxes(fixpoints, cell_to_seed, completed_cells);
        }
        std::cerr << "COMPLETED ROUND " << i << " of fixed point extraction\n\n\n\n";
    }
    std::cerr << fixpoints.size() << " total fixpoints found\n";
    
    if (std::strcmp(topo, "none")) {
        std::fstream _file(topo, std::ios::out);
        _file << fixpoints.size() << '\n';
        for (int i=0 ; i<fixpoints.size() ; ++i) {
            const fixpoint& fp = fixpoints[i];
            _file << (fp.saddle ? "1" : "0") << " " << fp.K << " "
                  << fp.pos[0] << " " << fp.pos[1];
            if (fp.saddle) {
                _file << " " << fp.evec[0][0] << " " << fp.evec[0][1]
                      << " " << fp.evec[1][0] << " " << fp.evec[1][1];
            }
            _file << '\n';
        }
        _file.close();
    }
}

void filter_chains()
{
    // group fixpoints by period and by type
    std::vector<fixpoint>   saddles[maxp+1];
    std::vector<fixpoint>   centers[maxp+1];
    for (int i=0 ; i<_fixpoints.size() ; ++i) {
        const fixpoint& fp = _fixpoints[i];
        if (fp.saddle) {
            saddles[fp.K].push_back(fp);
        } else {
            centers[fp.K].push_back(fp);
        }
    }
    std::cerr << "Summary:\n";
    for (int i=1 ; i<=maxp ; ++i) {
        std::cerr << "period " << i << ": " << saddles[i].size() << " saddles and "
                  << centers[i].size() << " centers\n";
    }
    
    for (int i=1 ; i<=maxp ; ++i) {
        link_chains(saddles[i]);
        link_chains(centers[i]);
    }
    // discard redundant chains
    std::vector<std::vector<fixpoint> > unique_chains;
    
    std::set<int> removed;
    for (int i=0 ; i<_chains.size() ; ++i) {
        if (removed.find(i) != removed.end()) {
            continue;
        }
        const std::vector<fixpoint>& chain = _chains[i];
        int p = chain.size();
        for (int j=i+1 ; j<_chains.size() ; ++j) {
            if (removed.find(j) != removed.end()) {
                continue;
            }
            if (_chains[j].size() != p) {
                continue;
            }
            double d = chain_distance(chain, _chains[j]);
            if (d < nvis::norm(_plane_spacing)) {
                removed.insert(j);
            }
        }
        unique_chains.push_back(chain);
    }
    std::cerr << "of the " << _chains.size() << " chains that were identified "
              << unique_chains.size() << " chains were unique\n";
    _chains.swap(unique_chains);
    
    std::vector<fixpoint> fps;
    for (int i=0 ; i<_chains.size() ; ++i) {
        for (int j=0 ; j<_chains[i].size() ; ++j) {
            fps.push_back(_chains[i][j]);
        }
    }
    _fixpoints.swap(fps);
}

void compute_manifolds()
{

    std::vector<fp_chain> all_p_chains[maxp+1];
    for (int i=0 ; i<_chains.size() ; ++i) {
        int p = _chains[i].size();
        if (!_chains[i].front().saddle) {
            continue;
        }
        all_p_chains[p].push_back(fp_chain(_chains[i], _plane_metric));
    }
    
    standard_map pmap(_k);
    
    for (int p=minsep ; p<=maxsep ; ++p) {
        std::cerr << "processing period " << p << std::endl;
        
        std::set<int> invalid;
// #pragma openmp parallel
        {
// #pragma omp for schedule(dynamic,1)
            for (int i=0 ; i<all_p_chains[p].size() ; ++i) {
                std::cerr << "processing saddle chain #" << i << std::endl;
                
                // if (invalid.find(i) != invalid.end()) continue;
                std::vector<separatrix> seps;
                
                map_type* amap = pmap.clone();
                
                std::cerr << "calling compute separatrix\n";
                bool ok = compute_separatrix(all_p_chains[p], seps,
                                             i, *amap, _plane_metric, 0.1*nvis::norm(_plane_spacing),
                                             _params);
                                             
                std::cerr << "separatrices computation completed\n";
                
                if (ok && seps.size()) {
// #pragma openmp atomic
                    {
                        invalid.insert(i);
                        invalid.insert(seps[0].end.first);
                    }
// #pragma openmp atomic
                    _separatrices.push_back(seps);
                } else {
                    std::cerr << "separatrix computation was invalid\n";
                }
            }
        }
    }
    
    // // downsample manifolds
    // std::vector<nvis::vec2> down;
    // for (int i=0 ; i<_separatrices.size() ; ++i) {
    //  for (int j=0 ; j<_separatrices[i].size() ; ++j) {
    //      if (_separatrices[i].size() < 500) continue;
    //      down.clear();
    //      const std::vector<nvis::vec2>& steps = _separatrices[i][j].manifold;
    //      double length = _separatrices[i][j].length;
    //      double dl = std::min(0.5, length/100.);
    //      down.push_back(_plane_metric.modulo(steps[0]));
    //      for (int k=1; k<steps.size()-1 ; ++k) {
    //          const nvis::vec2& x = steps[k];
    //          if (_plane_metric.distance(down.back(), x) > dl) {
    //              down.push_back(_plane_metric.modulo(x));
    //          }
    //          if (!_plane_bounds.inside(down.back())) {
    //              std::cerr << "WARNING: last inserted position: " << down.back()
    //                  << " is invalid\n";
    //          }
    //      }
    //      down.push_back(_plane_metric.modulo(steps.back()));
    //      _separatrices[i][j].manifold.swap(down);
    //  }
    // }
}

// void save_chains(const std::string& name) {
//     std::fstream file(name.c_str(), std::ios::out);
//
//     file << _chains.size() << '\n';
//     for (int i=0 ; i<_chains.size() ; ++i) {
//         file << _chains[i].size();
//         if (_chains[i].front().saddle)
//             file << " 1" << '\n';
//         else file << " 0" << '\n';
//         for (int j=0 ; j<_chains[i].size() ; ++j) {
//             nvis::vec2 x = _chains[i][j].pos;
//             file << x[0] << " " << x[1];
//             if (_chains[i].front().saddle) {
//                 nvis::vec2 e0 = _chains[i][j].evec[0];
//                 file << " " << e0[0] << " " << e0[1];
//                 nvis::vec2 e1 = _chains[i][j].evec[1];
//                 file << " " << e1[0] << " " << e1[1];
//             }
//             file << '\n';
//         }
//     }
//
//     file.close();
// }

bool load_fixpoints(const std::string& name)
{
    std::fstream file(name.c_str(), std::ios::in);
    if (!file) {
        return false;
    }
    
    _chains.clear();
    int N;
    file >> N;
    
    std::cerr << "there are " << N << " fixed points\n";
    
    _fixpoints.clear();
    
    for (int i=0 ; i<N ; ++i) {
        fixpoint fp;
        file >> fp.saddle >> fp.K >> fp.pos[0] >> fp.pos[1];
        if (fp.saddle) {
            file >> fp.evec[0][0] >> fp.evec[0][1] >> fp.evec[1][0] >> fp.evec[1][1];
        }
        std::cerr << "loaded fixpoint at " << fp.pos << std::endl;
        _fixpoints.push_back(fp);
    }
    file.close();
    
    return true;
}

bool load_chains(const std::string& name)
{
    std::fstream file(name.c_str(), std::ios::in);
    if (!file) {
        return false;
    }
    
    _chains.clear();
    int N;
    file >> N;
    
    std::cerr << "there are " << N << " chains\n";
    
    for (int i=0 ; i<N ; ++i) {
        _chains.push_back(std::vector<fixpoint>());
        int p;
        file >> p;
        std::cerr << "chain #" << i << " is of period " << p << std::endl;
        std::vector<fixpoint>& chain = _chains.back();
        chain.resize(p);
        int issaddle;
        file >> issaddle;
        if (issaddle) {
            std::cerr << "chain #" << i << " is a saddle\n";
        } else {
            std::cerr << "chain #" << i << " is a center\n";
        }
        for (int j=0 ; j<p ; ++j) {
            fixpoint& fp = chain[j];
            fp.K = p;
            fp.saddle = issaddle;
            file >> fp.pos[0] >> fp.pos[1];
            std::cerr << "read " << fp.pos << std::endl;
            if (issaddle) {
                file >> fp.evec[0][0] >> fp.evec[0][1] >> fp.evec[1][0] >> fp.evec[1][1];
                std::cerr << "read eigenvectors " << fp.evec[0] << " and "
                          << fp.evec[1] << std::endl;
            }
        }
    }
    
    file.close();
    return true;
}

bool save_chains(const std::string& name)
{
    std::ostringstream os;
    os << name << "-chains.txt";
    
    std::fstream file(os.str().c_str(), std::ios::out);
    if (!file) {
        return false;
    }
    
    file << _chains.size() << std::endl;
    
    std::cerr << "there are " << _chains.size() << " chains\n";
    
    for (int i=0 ; i<_chains.size() ; ++i) {
        bool saddle = _chains[i][0].saddle;
        
        file << "c " << _chains[i].size() << " " << (saddle ? "1" : "0") << std::endl;
        for (int n=0 ; n<_chains[i].size() ; ++n) {
            const fixpoint& fp = _chains[i][n];
            file << "p " << fp.pos[0] << " " << fp.pos[1];
            if (saddle)
                file << " " << fp.evec[0][0] << " "
                     << fp.evec[0][1] << " "
                     << fp.evec[1][0] << " "
                     << fp.evec[1][1];
            file << std::endl;
        }
    }
    
    file.close();
    return true;
}

bool save_separatrices(const std::string& name)
{

    std::ostringstream os;
    os << name << "-separatrices.txt";
    
    std::fstream file(os.str().c_str(), std::ios::out);
    if (!file) {
        return false;
    }
    
    file << _separatrices.size() << std::endl;
    
    std::cerr << "there are " << _separatrices.size() << " chains of separatrices\n";
    for (int i=0 ; i<_separatrices.size() ; ++i) {
        const std::vector<separatrix>& seps = _separatrices[i];
        file << "s " << seps.size() << std::endl;
        for (int j=0 ; j<seps.size() ; ++j) {
            const separatrix& s = seps[j];
            file << "m " << s.start.first << " " << s.start.second << " " << s.end.first << " " << s.end.second << " " << s.manifold.size() << " " << s.length << std::endl;
            for (int k=0 ; k<s.manifold.size() ; ++k) {
                file << "p " << s.manifold[k][0] << " " << s.manifold[k][1] << std::endl;
            }
        }
    }
    
    file.close();
}

void check_cell(const nvis::ivec2& cell_id)
{
    _params.vectors.clear();
    _params.record = true;
    _params.verbose = true;
    _params.lmin = 0.99 * std::min(_plane_spacing[0], _plane_spacing[1])/16.;
    
    std::vector<int> periods;
    period_range(periods, *_dataset, cell_id, _valid_rationals);
    std::cerr << "periods to be considered in this cell are:\n";
    for (int i=0 ; i<periods.size() ; ++i) {
        std::cerr << periods[i] << ", ";
    }
    std::cerr << std::endl;
    
    map_type pmap(_k);
    std::vector<std::pair<int, int> > indices;
    try {
        process_cell<map_type>(indices, *_dataset, *_plane, cell_id, pmap, _index_dx, _params);
    } catch(rational_surface_found& e) {}
    
    std::vector<nvis::vec2> chain;
    map_type* amap = pmap.clone();
    
    fixpoint fp;
    nvis::bbox2 bounds;
    bounds.min() = (*_plane)(cell_id);
    bounds.max() = (*_plane)(cell_id + nvis::ivec2(1,1));
    nvis::vec2 x = 0.5*(bounds.min() + bounds.max());
    bool found = meta_newton_stdmap<map_type>(*amap, _plane_metric, bounds, x, 5, periods[0],
                 fp, chain, 1.0e-5, true, 50, false);
                 
    _params.record = false;
    _params.verbose = false;
    _params.lmin = 0;
}

void check_orbit(const nvis::vec2& x)
{
    standard_map pmap(_k);
    
    _orbits.push_back(std::pair<orbit_type, nvis::fvec3>());
    iterate<map_type>(_last, _orbits.back().first, pmap, niter);
    int period = best_period(_orbits.back().first, maxp, _plane_metric);
    
    _orbits.back().second = (*_cmap)(period);
    std::cerr << "best period at " << x << " is " << period << '\n';
    
    std::vector<nvis::vec2> steps;
    std::cerr << "steps:\n";
    for (int i=0 ; i<_orbits.back().first.size() ; ++i) {
        steps.push_back(_orbits.back().first[i]);
        const nvis::vec2& _x = steps.back();
        std::cerr << "\n" << _x[0] << ", " << _x[1];
    }
    std::cerr << '\n';
    
    std::cerr << "steps when lifted once:\n";
    nvis::vec2 y = x + nvis::vec2(0,1);
    std::vector<nvis::vec2> _steps;
    for (int i=0 ; i<niter+1 ; ++i) {
        std::cerr << y[0] << ", " << y[1] << '\n';
        _steps.push_back(y);
        y = pmap.map(y);
    }
    std::cerr << "corresponding safety factor: " << safety_factor(_steps, _plane_metric);
    
    std::cerr << "steps when lifted twice:\n";
    y = x + nvis::vec2(0,2);
    _steps.clear();
    for (int i=0 ; i<niter+1 ; ++i) {
        std::cerr << y[0] << ", " << y[1] << '\n';
        _steps.push_back(y);
        y = pmap.map(y);
    }
    std::cerr << "corresponding safety factor: " << safety_factor(_steps, _plane_metric);
    
    std::cerr << "steps when lifted minus once:\n";
    y = x + nvis::vec2(0,-1);
    _steps.clear();
    for (int i=0 ; i<niter+1 ; ++i) {
        std::cerr << y[0] << ", " << y[1] << '\n';
        _steps.push_back(y);
        y = pmap.map(y);
    }
    std::cerr << "corresponding safety factor: " << safety_factor(_steps, _plane_metric);
    
    double sf = safety_factor(steps, _plane_metric);
    std::cerr << "safety factor at " << x << " is " << sf << '\n';
    
    double sf2 = safety_factor_bounded(steps, _plane_metric);
    std::cerr << "bounded safety factor at " << x << " is " << sf2 << '\n';
    
    double sf3 = safety_factor_stdmap(steps, _plane_metric);
    std::cerr << "STDMAP safety factor at " << x << " is " << sf3 << '\n';
    
    boost::rational<int> q = rational_approx(sf, maxp);
    std::cerr << "best rational approximation within prescribed bounds is " << q << std::endl;
    
    std::cerr << "approximation error = " << fabs(sf-value<int, double>(q)) << std::endl;
    
    try {
        nvis::vec2 y = pmap.map(x, period);
        _params.vectors.push_back(std::pair<nvis::vec2, map_analysis_param::tagged_vector_type>(x, map_analysis_param::tagged_vector_type(_plane_metric.displacement(x, y), period)));
    } catch(...) {}
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
    
    // draw mesh
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(0.2, 0.2, 0.2);
    glLineWidth(1.0);
    for (int i=0 ; i<_params.edges.size()/2 ; ++i) {
        glBegin(GL_LINES);
        nvis::vec2 x = _params.edges[2*i  ];
        nvis::vec2 y = _params.edges[2*i+1];
        glVertex2f(x[0], x[1]);
        glVertex2f(y[0], y[1]);
        glEnd();
    }
    
    // draw cells containing a saddle
    glLineWidth(2.0);
    glColor3f(1, 0, 0);
    for (int i=0 ; i<_saddle_cells.size() ; ++i) {
        const nvis::ivec2& id = _saddle_cells[i];
        glBegin(GL_LINE_STRIP);
        nvis::vec2 x = mesh(id);
        glVertex2f(x[0], x[1]);
        x = mesh(id + nvis::ivec2(1,0));
        glVertex2f(x[0], x[1]);
        x = mesh(id + nvis::ivec2(1,1));
        glVertex2f(x[0], x[1]);
        x = mesh(id + nvis::ivec2(0,1));
        glVertex2f(x[0], x[1]);
        x = mesh(id);
        glVertex2f(x[0], x[1]);
        glEnd();
    }
    
    // draw cells containing a center
    glColor3f(0, 0, 1);
    for (int i=0 ; i<_center_cells.size() ; ++i) {
        const nvis::ivec2& id = _center_cells[i];
        glBegin(GL_LINE_STRIP);
        nvis::vec2 x = mesh(id);
        glVertex2f(x[0], x[1]);
        x = mesh(id + nvis::ivec2(1,0));
        glVertex2f(x[0], x[1]);
        x = mesh(id + nvis::ivec2(1,1));
        glVertex2f(x[0], x[1]);
        x = mesh(id + nvis::ivec2(0,1));
        glVertex2f(x[0], x[1]);
        x = mesh(id);
        glVertex2f(x[0], x[1]);
        glEnd();
    }
    
    // draw orbits computed during automatic analysis
    if (show_orb) {
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_BLEND);
        glPointSize(pt_sz);
        glBegin(GL_POINTS);
        for (int i=0 ; i<_params.orbits.size() ; ++i) {
            const nvis::fvec3& c = (*_cmap)(_params.orbits[i].second);
            glColor3f(c[0], c[1], c[2]);
            for (int j=0 ; j<_params.orbits[i].first.size() ; ++j) {
                nvis::vec2 x = _plane_metric.modulo(_params.orbits[i].first[j]);
                glVertex2f(x[0], x[1]);
            }
        }
        glEnd();
    }
    
    // draw interactively computed orbits
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(pt_sz);
    glBegin(GL_POINTS);
    for (int i=0 ; i<_orbits.size() ; ++i) {
        const nvis::fvec3& c = _orbits[i].second;
        glColor3f(c[0], c[1], c[2]);
        for (int j=0 ; j<_orbits[i].first.size() ; ++j) {
            nvis::vec2 x = _plane_metric.modulo(_orbits[i].first[j]);
            glVertex2f(x[0], x[1]);
        }
    }
    glEnd();
    
    if (false) {
        // draw fixed points
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_BLEND);
        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(1, 0, 0);
        for (int i=0 ; i<_saddles.size() ; ++i) {
            nvis::vec2 x = _saddles[i];
            glVertex2f(x[0], x[1]);
        }
        glColor3f(0, 0, 1);
        for (int i=0 ; i<_centers.size() ; ++i) {
            nvis::vec2 x = _centers[i];
            glVertex2f(x[0], x[1]);
        }
        glEnd();
    }
    
    // draw vector associated with each vertex
    if (show_vec) {
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(1.0);
        const double cos_alpha = cos(M_PI / 12.);
        const double sin_alpha = sin(M_PI / 12.);
        glBegin(GL_LINES);
        for (int n=0 ; n<_dataset->size() ; ++n) {
            int i = n % _plane_res[0];
            int j = n / _plane_res[0];
            nvis::ivec2 id(i,j);
            const orbit_data& d = (*_dataset)(id);
            if (d.sf == 0) {
                continue;
            }
            nvis::vec2 p0 = (*_plane)(id);
            nvis::vec2 v = d.vector(d.p);
            if (nvis::norm(v) == 0) {
                continue;
            }
            int p = d.p;
            nvis::vec2 p1 = p0 + v;
            nvis::fvec3 c = (*_cmap)(p);
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
    
    if (show_edges) {
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glColor3f(1, 0, 0);
        glLineWidth(2.0);
        glBegin(GL_LINES);
        for (int i=0 ; i<_failed_edges.size() ; ++i) {
            nvis::ivec2 i0 = _failed_edges[i].object()[0];
            nvis::ivec2 i1 = _failed_edges[i].object()[1];
            nvis::vec2 x = (*_plane)(i0);
            nvis::vec2 y = (*_plane)(i1);
            glVertex2f(x[0], x[1]);
            glVertex2f(y[0], y[1]);
        }
        glEnd();
    }
    
    // draw vectors measured during poincare index computation
    if (true) {
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(1.0);
        const double cos_alpha = cos(M_PI / 12.);
        const double sin_alpha = sin(M_PI / 12.);
        glBegin(GL_LINES);
        for (int n=0 ; n<_params.vectors.size() ; ++n) {
            nvis::vec2 p0 = _params.vectors[n].first;
            const std::pair<nvis::vec2, int>& tv = _params.vectors[n].second;
            const nvis::vec2& v = tv.first;
            if (nvis::norm(v) == 0) {
                continue;
            }
            int p = tv.second;
            nvis::vec2 p1 = p0 + v;
            nvis::fvec3 c = (*_cmap)(p);
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
    
    // draw fixpoints
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(5);
    glBegin(GL_POINTS);
    for (int i=0 ; i<_fixpoints.size() ; ++i) {
        fixpoint& fp = _fixpoints[i];
        if (fp.saddle) {
            glColor3f(1, 0, 0);
        } else {
            glColor3f(0, 0, 1);
        }
        fp.pos = _plane_metric.modulo(fp.pos);
        glVertex3f(fp.pos[0], fp.pos[1], 0.5);
    }
    // }
    glEnd();
    
    // // draw chains
    //     glEnable(GL_BLEND);
    //     glEnable(GL_LINE_SMOOTH);
    //     glLineWidth(1.0);
    //     glColor3f(1,1,0);
    //     glBegin(GL_LINES);
    //     for (int i=0 ; i<_chains.size() ; ++i) {
    //         for (int j=0 ; j < _chains[i].size()-1 ; ++j) {
    //      const nvis::vec2& x0 = _chains[i][j].pos;
    //      const nvis::vec2& x1 = _chains[i][j+1].pos;
    //      glVertex2f(x0[0], x0[1]);
    //      glVertex2f(x1[0], x1[1]);
    //         }
    //     }
    // glEnd();
    
    // draw eigenvectors
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);
    glColor3f(1,1,1);
    glBegin(GL_LINES);
    // for (int i=0 ; i<_chains.size() ; ++i) {
    //      for (int j=0 ; j<_chains[i].size() ; ++j) {
    //          const fixpoint& fp = _chains[i][j];
    for (int i=0 ; i<_fixpoints.size() ; ++i) {
        fixpoint& fp = _fixpoints[i];
        if (fp.saddle) {
            nvis::vec2 x = fp.pos + 0.0025*fp.evec[0];
            glVertex2f(x[0], x[1]);
            x -= 0.005*fp.evec[0];
            glVertex2f(x[0], x[1]);
            x = fp.pos + 0.0025*fp.evec[1];
            glVertex2f(x[0], x[1]);
            x -= 0.005*fp.evec[1];
            glVertex2f(x[0], x[1]);
        }
        // }
    }
    glEnd();
    
    // draw separatrices
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINES);
    for (int i=0 ; i<_separatrices.size() ; ++i) {
        for (int j=0 ; j<_separatrices[i].size() ; ++j) {
            const separatrix& sep = _separatrices[i][j];
            const std::vector<nvis::vec2>& pts = sep.manifold;
            for (int n=0 ; n<pts.size()-1 ; ++n) {
                glVertex2f(pts[n][0], pts[n][1]);
                nvis::vec2 x = pts[n] + _plane_metric.displacement(pts[n], pts[n+1]);
                glVertex2f(x[0], x[1]);
            }
        }
    }
    glEnd();
    
    
    // draw locations where period identification failed
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(2);
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_POINTS);
    for (int i=0 ; i<_problematic_seeds.size() ; ++i) {
        const nvis::vec2& x = _problematic_seeds[i];
        glVertex3f(x[0], x[1], 0.05);
    }
    glEnd();
    
    // draw rational surfaces
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    for (int i=0 ; i<_rational_segments.size() ; ++i) {
        int period = _rational_segments[i].sf.numerator();
        nvis::fvec3 c = (*_cmap)(period);
        glColor3f(c[0], c[1], c[2]);
        nvis::vec2 x = _plane_metric.modulo(_rational_segments[i].pt[0].second);
        nvis::vec2 y = x + _plane_metric.displacement(x, _rational_segments[i].pt[1].second);
        glVertex2f(x[0], x[1]);
        glVertex2f(y[0], y[1]);
    }
    glEnd();
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
    standard_map pmap(_k);
    
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    
    if (button == GLUT_MIDDLE_BUTTON) {
    
        nvis::ivec2 cell_id = _plane->local_coordinates(wc).first;
        
        std::cerr << "computing index in cell " << cell_id << std::endl;
        check_cell(cell_id);
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
    } else if (key == 's') {
        // save
        save_to_file = true;
    } else if (key == 'i') {
        // iterate
        standard_map pmap(_k);
        
        _orbits.push_back(std::pair<orbit_type, nvis::fvec3>());
        iterate<map_type>(_last, _orbits.back().first, pmap, 200);
        _orbits.back().second = (*_cmap)(best_period(_orbits.back().first, maxp, _plane_metric));
        display();
    } else if (key == 'x') {
        // erase
        _orbits.clear();
        _params.vectors.clear();
        display();
    } else if (key == 'q') {
        check_orbit(wc);
    } else if (key = 'j') {
        standard_map pmap(_k);
        orbit_type tmp;
        iterate<map_type>(_last, tmp, pmap, 200);
        int period = best_period(tmp, maxp, _plane_metric);
        std::cerr << "best period at " << _last << " is " << period << std::endl;
        std::pair<nvis::vec2, nvis::mat2> res = pmap.map_and_jacobian(_last, period);
        nvis::mat2 r = pmap.fjacobian(_last[0], _last[1]);
        std::cerr << "1-jacobian at " << _last << " is " << r << std::endl;
        nvis::mat2& J = res.second;
        fixpoint fp;
        fp.pos = x;
        linear_analysis(J, period, _last, fp);
        std::cerr << "linear analysis at " << _last << ":" << std::endl
                  << "J = " << J << std::endl
                  << (fp.saddle ? "saddle" : "center") << std::endl;
    }
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    _last = nvis::vec2(-1,-1);
    
    initialize(argc, argv);
    save_to_file = false;
    
    width = 1200;
    height = 1200;
    
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
    
    std::vector<nvis::fvec3> colors;
    xavier::spiral_scale(colors, maxp, 1);
    std::vector<int> reference_values(maxp);
    for (int i=1 ; i<=maxp ; ++i) {
        reference_values[i-1] = i;
    }
    _cmap = new xavier::discrete_color_map<int>(reference_values, colors);
    
    init();
    
    if (!nogfx) {
        // draw mesh
        double x = _plane_bounds.min()[0];
        for (int i=0 ; i<_plane_res[0] ; ++i, x+=_plane_spacing[0]) {
            _params.edges.push_back(nvis::vec2(x, _plane_bounds.min()[1]));
            _params.edges.push_back(nvis::vec2(x, _plane_bounds.max()[1]));
        }
        double y = _plane_bounds.min()[1];
        for (int j=0 ; j<_plane_res[1] ; ++j, y+=_plane_spacing[1]) {
            _params.edges.push_back(nvis::vec2(_plane_bounds.min()[0], y));
            _params.edges.push_back(nvis::vec2(_plane_bounds.max()[0], y));
        }
    }
    
    int npoints = _plane_res[0] * _plane_res[1];
    
    standard_map pmap(_k);
    
    bool done = false;
    if (strcmp(topo, "none")) {
        done = load_fixpoints(topo);
        if (!done) {
            std::cerr << "topological information will be saved to " << topo << std::endl;
        }
    }
    
    nvis::timer _timer;
    if (!done || show_vec) {
        sample_raster<map_type>(*_dataset, *_plane, pmap, _params, false);
        double dt = _timer.elapsed();
        std::cerr << "map sampling at grid vertices took " << dt << " seconds "
                  << "(" << _plane_res[0]*_plane_res[1]/dt << "Hz)\n";
                  
        /*
        // check resulting values
        double *qs = (double*)calloc(_plane_res[0] * _plane_res[1], sizeof(double));
        for (int n=0 ; n<npoints ; ++n) {
            int i = n % _plane_res[0];
            int j = n / _plane_res[0];
            nvis::ivec2 id(i,j);
            if ((*_dataset)(id).sf == 0) {
                _problematic_seeds.push_back((*_plane)(id));
            }
            qs[n] = (*_dataset)(id).sf;
        }
        std::vector<size_t> s(2);
        s[0] = _plane_res[0];
        s[1] = _plane_res[1];
        std::vector<double> sp(2);
        sp[0] = _plane_spacing[0];
        sp[1] = _plane_spacing[1];
        xavier::writeNrrd(qs, "qs.nrrd", nrrdTypeDouble, s, sp);
        */
    }
    
    if (!done && cell_analysis) {
        std::vector<p_cell_type> cells;
        compute_cell_periods(cells);
        compute_edge_rotation(cells);
        compute_cell_index(cells);
        compute_fixed_points();
        
        std::cerr << "total processing time was " << _timer.elapsed() << std::endl;
    }
    filter_chains();
    
    compute_manifolds();
    if (strcmp(sepfile, "none")) {
        save_chains(sepfile);
        std::cerr << "fixpoint chains saved" << std::endl;
        save_separatrices(sepfile);
        std::cerr << "separatrices chains saved" << std::endl;
    }
    
    if (!nogfx) {
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
