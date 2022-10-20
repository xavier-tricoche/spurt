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
#include <teem/hest_helper.hpp>
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


using namespace spurt;

// -------------------------
//
//      Data Structure
//
// -------------------------
typedef grid<double, 3>                                     volume_type;
typedef grid<double, 2>                                     plane_type;
typedef raster_data<nvis::vec3, double, 3>                  field_type;
typedef std::pair<nvis::vec2, int>                          map_data_type;
typedef raster_data<map_data_type, double, 2>               dataset_type;
typedef xmt_poincare_map<spurt::map::wrapper<field_type> > map_type;

field_type*                 _field;
std::vector<nvis::vec3>     _vectors;
volume_type*                _volume;
nvis::ivec3                 _volume_res;
nvis::bbox3                 _volume_bounds;
nvis::vec3                  _volume_spacing;
plane_type*                 _plane;
nvis::ivec2                 _plane_res;
nvis::bbox2                 _plane_bounds;
spurt::map_metric          _plane_metric;
nvis::vec2                  _plane_spacing;
dataset_type*               _dataset;


// -------------------------
//
//              UI
//
// -------------------------
char*    in, *phys, *outmap;
int     niter, width, height, maxp, show_orb, show_vec, __res[2];
double  eps;
float   pt_sz;
bool    logical;

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
std::vector<iteration_data>     _iteration_data;

void check_data_structure(const field_type& field)
{
    nvis::ivec3 res = field.get_grid().dimensions();
    nvis::bbox3 bounds = field.get_grid().bounds();
    nvis::vec3 step = field.get_grid().spacing();
    
    nvis::ivec3 did[] = {
        nvis::ivec3(0,0,0), nvis::ivec3(1,0,0),
        nvis::ivec3(1,1,0), nvis::ivec3(0,1,0),
        nvis::ivec3(0,0,1), nvis::ivec3(1,0,1),
        nvis::ivec3(1,1,1), nvis::ivec3(0,1,1),
    };
    
    double errf = 0;
    double errJ = 0;
    for (int n=0 ; n<10 ; ++n) {
    
        // pick a cell at random
        int i = floor(drand48()*res[0]);
        int j = floor(drand48()*res[1]);
        int k = floor(drand48()*res[2]);
        nvis::ivec3 cell_id(i,j,k);
        
        // get values associated with vertices
        nvis::vec3 val[8];
        for (int k=0 ; k<8 ; ++k) {
            val[k] = field(cell_id + did[k]);
        }
        
        nvis::vec3 x0 = field.get_grid()(cell_id);
        
        // probe field at a number of random locations in that cell
        for (int k=0 ; k<1000 ; ++k) {
            double u = drand48();
            double v = drand48();
            double w = drand48();
            
            nvis::vec3 x = x0 + nvis::vec3(u,v,w)*step;
            nvis::vec3 f =
                (1.-u)*(1.-v)*(1.-w)*val[0] +
                u*(1.-v)*(1.-w)*val[1] +
                u*v*(1.-w)*val[2] +
                (1.-u)*v*(1.-w)*val[3] +
                (1.-u)*(1.-v)*w*val[4] +
                u*(1.-v)*w*val[5] +
                u*v*w*val[6] +
                (1.-u)*v*w*val[7];
                
            nvis::vec3 ddx = (
                                 -(1.-v)*(1.-w)*val[0] +
                                 (1.-v)*(1.-w)*val[1] +
                                 v*(1.-w)*val[2] -
                                 v*(1.-w)*val[3] -
                                 (1.-v)*w*val[4] +
                                 (1.-v)*w*val[5] +
                                 v*w*val[6] -
                                 v*w*val[7] ) / step[0];
                                 
            nvis::vec3 ddy = (
                                 -(1.-u)*(1.-w)*val[0] -
                                 u*(1.-w)*val[1] +
                                 u*(1.-w)*val[2] +
                                 (1.-u)*(1.-w)*val[3] -
                                 (1.-u)*w*val[4] -
                                 u*w*val[5] +
                                 u*w*val[6] +
                                 (1.-u)*w*val[7] ) / step[1];
                                 
            nvis::vec3 ddz = (
                                 -(1.-u)*(1.-v)*val[0] -
                                 u*(1.-v)*val[1] -
                                 u*v*val[2] -
                                 (1.-u)*v*val[3] +
                                 (1.-u)*(1.-v)*val[4] +
                                 u*(1.-v)*val[5] +
                                 u*v*val[6] +
                                 (1.-u)*v*val[7] ) / step[2];
                                 
            nvis::mat3 J;
            for (int r=0 ; r<3 ; ++r) {
                J(r,0) = ddx[r];
                J(r,1) = ddy[r];
                J(r,2) = ddz[r];
            }
            
            nvis::vec3 ff = field.interpolate(x);
            nvis::mat3 JJ = nvis::transpose(nvis::mat3(field.derivative(x)));
            
            errf += nvis::norm(ff-f);
            errJ += nvis::norm(JJ-J);
        }
    }
    
    std::cerr << "cumulative interpolation error = " << errf << '\n'
              << "cumulative derivation error = " << errJ << '\n';
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
    hestOptAdd(&hopt, "i",      "input",                airTypeString,  1,  1,  &in,        NULL,           "input file name (NRRD)");
    hestOptAdd(&hopt, "map",    "map file",             airTypeString,  1,  1,  &outmap,    NULL,           "input/output file name containig map approximation (NRRD)");
    hestOptAdd(&hopt, "n",      "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0,  1,  &maxp,      "50",           "max considered period");
    hestOptAdd(&hopt, "e",      "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "r",      "resolution",           airTypeInt,     2,  2,  &__res,     NULL,           "plane resolution");
    hestOptAdd(&hopt, "ps",     "point size",           airTypeFloat,   0,  1,  &pt_sz,     "2",            "point size for display");
    hestOptAdd(&hopt, "phys",   "phys space",           airTypeString,  0,  1,  &phys,      "none",         "HDF5 reference file for mapping");
    hestOptAdd(&hopt, "orb",    "show orbits",          airTypeInt,     0,  0,  &show_orb,  NULL,           "show automatically generated orbits");
    hestOptAdd(&hopt, "vec",    "show vectors",         airTypeInt,     0,  0,  &show_vec,  NULL,           "show vertex vectors");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "approximated map analysis",
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
spurt::discrete_color_map<int>*        _cmap;
spurt::logical2physical*               _converter;
spurt::map_analysis_param              _params;
std::vector<nvis::ivec2>                _saddle_cells, _center_cells;
nvis::vec2                              _last;
std::vector<nvis::vec2>                 _problematic_seeds;
std::vector<color_orbit_type>           _rational_surfaces;

// -------------------------
//
//          Analysis
//
// -------------------------
typedef rational_surface_found::edge_type   edge_type;
typedef rational_surface_found::edge_point  edge_point;

struct rational_segment {
    edge_point pt[2];
    boost::rational<int> sf;
};
std::vector<rational_segment>           _rational_segments;
std::vector<fixpoint>                   _fixpoints;
std::vector<int>                        _tagged_cells;
double                                  _index_dx;

enum {
    FIXED_POINT = 0,
    CONVERGED,
    DEGENERATE
};

// void construct_rational_curves() {
//  std::multimap<boost::rational<int>, rational_segment> per_period;
//  typedef std::pair<boost::rational<int>, rational_segment> pair_type;
//
//  for (int i=0 ; i<_rational_segments.size() ; ++i) {
//      per_period.insert(pair_type(_rational_segments[i].sf, _rational_segments[i]));
//  }
//
//
// }

void compute_map(std::vector<orbit_data>& dataset,
                 const grid<double, 2>& mesh, const map_type& pmap,
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
    
    std::fstream input(outmap, std::ios::in);
    if (input) {
        input.close();
        Nrrd* nin = nrrdNew();
        nin = spurt::readNrrd(outmap);
        
        if (nin->dim != 3 || nin->axis[1].size != resolution[0] ||
                nin->axis[2].size != resolution[1] ||
                nin->axis[0].size != 2*(niter+1)+1) {
            std::cerr << "invalid file for given parameters\n";
            exit(0);
        }
        
        const double* data = (const double*)nin->data;
        dataset.resize(npoints);
        for (int j=0 ; j<resolution[1] ; ++j) {
            for (int i=0 ; i<resolution[0] ; ++i) {
                int n = i + j*resolution[0];
                dataset[n].steps.resize(niter+1);
                for (int k=0 ; k<niter+1 ; ++k) {
                    dataset[n].steps[k][0] = data[n*(2*niter + 3)+2*k  ];
                    dataset[n].steps[k][1] = data[n*(2*niter + 3)+2*k+1];
                    dataset[n].steps[k][2] = data[n*(2*niter + 3)+2*k+2];
                }
                dataset[n].set();
            }
        }
        
        nrrdNuke(nin);
        return;
    }
    
    input.close();
    
    int orbit_length = 2*(niter+1);
    int row_length = orbit_length + 1;
    double* __data = (double*)calloc(npoints*row_length, sizeof(double));
    
    std::cerr << "bounds = " << bounds << std::endl;
    std::cerr << "step = " << step << std::endl;
    std::cerr << "resolution = " << resolution << std::endl;
    
    npoints = resolution[0]*resolution[1];
    
    nvis::timer _timer;
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<std::pair<int, orbit_type> > cached_orbits[nbthreads];
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n = 0 ; n < npoints ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            int j = n / resolution[0];
            int i = n % resolution[0];
            nvis::ivec2 id(i,j);
            
            if (thread_id == 0) {
                std::ostringstream os;
                os << "\rcompleted " << n << " / " << npoints << " ("
                   << 100*n/npoints << "%)      \r" << std::flush;
                std::cout << os.str();
            }
            
            nvis::vec2 x0 = bounds.min() + step * nvis::vec2(i,j);
            
            map_type* amap = pmap.clone();
            std::vector<nvis::vec2> tmp;
            try {
                amap->map(x0, tmp, niter);
            } catch(...) {
                continue;
            }
            if (tmp.size() < niter) {
                continue;
            } else {
                cached_orbits[thread_id].push_back(std::pair<int, orbit_type>());
                cached_orbits[thread_id].back().first = n;
                orbit_type& steps = cached_orbits[thread_id].back().second;
                steps.push_back(x0);
                std::copy(tmp.begin(), tmp.end(), std::back_inserter(steps));
            }
        }
    }
    std::cout << '\n';
    std::cerr << "computation of " << npoints << " orbits took " << _timer.elapsed() << std::endl;
    
    dataset.resize(npoints);
    // gather
    for (int i=0 ; i<nbthreads ; ++i) {
        for (int j=0 ; j<cached_orbits[i].size() ; ++j) {
            int n = cached_orbits[i][j].first;
            dataset[n].steps.swap(cached_orbits[i][j].second);
            dataset[n].sf = safety_factor(dataset[n].steps, _plane_metric);
            int offset = n*row_length;
            for (int k=0 ; k<niter+1 ; ++k) {
                const nvis::vec2& y = dataset[n].steps[k];
                __data[offset+2*k  ] = y[0];
                __data[offset+2*k+1] = y[1];
            }
            __data[offset+row_length-1] = dataset[n].sf;
        }
    }
    size_t __size[3] = {row_length, resolution[0], resolution[1]};
    double __spc[3] = {airNaN(), _plane_spacing[0], _plane_spacing[1]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, __data, nrrdTypeDouble, 3, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave(outmap, nout, NULL);
    std::cerr << "map exported for reuse\n";
    nrrdNuke(nout);
}

static void init()
{
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    Nrrd* nin = nrrdNew();
    nin = spurt::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> _array;
    spurt::to_vector(_array, nin);
    _volume_res = nvis::ivec3(nin->axis[1].size,
                              nin->axis[2].size,
                              nin->axis[3].size);
    _volume_spacing = nvis::vec3(nin->axis[1].spacing,
                                 nin->axis[2].spacing,
                                 nin->axis[3].spacing);
    _volume = new volume_type(_volume_res, _volume_spacing,
                              nvis::fixed_vector<bool, 3>(false, true, true));
    _vectors.resize(_array.size() / 3);
    for(int i = 0 ; i < _array.size() / 3 ; ++i) {
        _vectors[i][0] = _array[3*i  ];
        _vectors[i][1] = _array[3*i+1];
        _vectors[i][2] = _array[3*i+2];
    }
    _field = new field_type(*_volume, _vectors);
    _field->verbose(false);
    
    // check_data_structure(*_field);
    
    double h = eps;
    
    _volume_bounds = _volume->bounds();
    _plane_bounds = nvis::bbox2(nvis::vec2(_volume_bounds.min()[1], _volume_bounds.min()[0]),
                                nvis::vec2(_volume_bounds.max()[1], _volume_bounds.max()[0]));
    _plane_res = nvis::ivec2(__res[0], __res[1]);
    _plane_metric.bounds() = _plane_bounds;
    _plane_metric.periodic(0) = true;
    _plane_metric.periodic(1) = false;
    int npoints = _plane_res[0] * _plane_res[1];
    
    if (logical) {
        GLUT_helper::box = _plane_bounds;
    } else {
        GLUT_helper::box = _converter->bounds();
    }
    std::cerr << "box set to " << GLUT_helper::box << std::endl;
    
    map_type pmap(*_field);
    pmap.precision(eps);
    _plane = new plane_type(_plane_res, _plane_bounds);
    _dataset = new dataset_type(*_plane, map_data_type(nvis::vec2(0,0), 0));
    _params.nb_iterations = niter;
    _params.max_period = maxp;
    _params.metric = _plane_metric;
    _plane_spacing = _plane->spacing();
    
    nvis::timer _timer;
    compute_map(_iteration_data, *_plane, pmap, niter);
    double dt = _timer.elapsed();
    std::cerr << "map sampling at grid vertices took " << dt << " seconds "
              << "(" << _plane_res[0]*_plane_res[1]/dt << "Hz)\n";
              
    // check resulting values
    for (int n=0 ; n<npoints ; ++n) {
        if (_iteration_data[n].sf == 0) {
            _problematic_seeds.push_back((*_plane)(id));
        }
    }
    
    nvis::ivec2 cellres(_plane_res[0]-1, _plane_res[1]-1);
    int ncells = cellres[0]*cellres[1];
    
    // identify cells that require further processing
    std::vector<std::pair<nvis::ivec2, int> > cached_sing_cells[nbthreads];
    std::vector<rational_segment> cached_rational_segments[nbthreads];
    
    _timer.restart();
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n=0 ; n<ncells ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            int j = n/cellres[0];
            int i = n%cellres[0];
            nvis::ivec2 cell_id(i,j);
            
            // _params.verbose = true;
            _params.record = false;
            _params.verbose = false;
            
            map_type* amap = pmap.clone();
            
            std::vector<std::pair<int, int> > indices;
            
            try {
                process_cell<map_type>(indices, *_dataset, *_plane, cell_id, *amap, _index_dx, _params);
            } catch(rational_surface_found& e) {
                cached_rational_segments[thread_id].push_back(rational_segment());
                rational_segment& rs = cached_rational_segments[thread_id].back();
                rs.pt[0] = e.xings[0];
                rs.pt[1] = e.xings[1];
                rs.sf = e.sf;
            }
            
            for (int i=0 ; i<indices.size() ; ++i) {
                const std::pair<int, int>& idx = indices[i];
                cached_sing_cells[thread_id].push_back(std::pair<nvis::ivec2, int>(cell_id, idx.first * idx.second));
            }
            
            if (thread_id == 0) {
                std::ostringstream os;
                os << "processed " << n << " out of " << ncells << " cells "
                   << "(" << 100*n/ncells << "%)         \r";
                std::cerr << os.str();
            }
        }
    }
    std::cerr << '\n';
    dt = _timer.elapsed();
    std::cerr << "cell index computation " << dt << " seconds "
              << "(" << cellres[0]*cellres[1]/dt << "Hz)\n";
              
    std::vector<std::pair<nvis::ivec2, int> > singular_cells;
    for (int i=0 ; i<nbthreads ; ++i) {
        std::copy(cached_sing_cells[i].begin(), cached_sing_cells[i].end(),
                  std::back_inserter(singular_cells));
        std::copy(cached_rational_segments[i].begin(), cached_rational_segments[i].end(),
                  std::back_inserter(_rational_segments));
    }
    int nsingular = singular_cells.size();
    std::cerr << "there are " << nsingular << " singular cells out of "
              << ncells << " (" << 100*nsingular/ncells << "%)\n";
    if (!nsingular) {
        return;
    }
    
    // construct_rational_curves();
    
    std::random_shuffle(singular_cells.begin(), singular_cells.end());
    
    // now process them in parallel
    std::vector<int> touched(ncells, false);
    
    double total_t = 0;
    int sing_step, nbrounds;
    // catch special cases
    if (nsingular < 100) {
        nbrounds = 1;
        sing_step = nsingular;
    } else {
        sing_step = nsingular / 20;
        nbrounds = ceil(nsingular/sing_step);
    }
    for (int r=0 ; r<nbrounds ; ++r) {
        int nstart = r*sing_step;
        int nend = std::min((r+1)*sing_step, nsingular);
        
        std::vector<fixpoint> cached_fp[nbthreads];
        std::vector<int> cached_touched[nbthreads];
        
        _timer.restart();
        
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for(int n=nstart ; n<nend ; ++n) {
                int thread_id = 0;
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                
                nvis::ivec2 cell_id = singular_cells[n].first;
                int period = abs(singular_cells[n].second);
                int index = singular_cells[n].second / period;
                int i = cell_id[0];
                int j = cell_id[1];
                int cell_n = i + cellres[0]*j;
                
                if (touched[cell_n]) {
                    continue;
                }
                cached_touched[thread_id].push_back(cell_n);
                
                nvis::vec2 x = _plane_bounds.min() + nvis::vec2((float)i + 0.5, (float)j+0.5)*_plane_spacing;
                std::vector<fixpoint> chain;
                map_type* amap = pmap.clone();
                
                bool found = linear_chain_analysis(*amap, _plane_metric, x, period,
                                                   chain, 2*nvis::norm(_plane_spacing), 1.0e-3);
                                                   
                std::ostringstream os;
                if (found) {
                    os << "found " << period << "-chain of type "
                       << (chain[0].saddle ? "saddle" : "center")
                       <<  std::endl;
                    for (int k=0 ; k<chain.size() ; ++k) {
                        nvis::ivec2 id = _plane->local_coordinates(chain[k].pos).first;
                        int m = id[0] + id[1]*(_plane_res[0]-1);
                        cached_touched[thread_id].push_back(m);
                    }
                    std::copy(chain.begin(), chain.end(), std::back_inserter(cached_fp[thread_id]));
                }
                std::cerr << os.str();
            }
        }
        for (int k=0 ; k<nbthreads ; ++k) {
            std::copy(cached_fp[k].begin(), cached_fp[k].end(), std::back_inserter(_fixpoints));
            for (int i=0 ; i<cached_touched[k].size() ; ++i) {
                touched[cached_touched[k][i]] = true;
            }
        }
        
        dt = _timer.elapsed();
        std::cerr << "processing " << nend - nstart << " cells took " << dt << " seconds"
                  << " (" << (nend-nstart)/dt << " Hz)\n";
        total_t += dt;
    }
    std::cerr << "total compute time for cell processing was " << total_t
              << " (" << ncells/total_t << " Hz)\n";
              
    // filter out duplicates
    std::vector<fixpoint> tmp;
    std::set<nvis::ivec2, nvis::lexicographical_order> unique;
    for (int i=0 ; i<_fixpoints.size() ; ++i) {
        nvis::ivec2 id = _plane->local_coordinates(_fixpoints[i].pos).first;
        if (unique.find(id) == unique.end()) {
            unique.insert(id);
            tmp.push_back(_fixpoints[i]);
        }
    }
    std::cerr << "out of " << _fixpoints.size() << ", " << tmp.size() << " were unique\n";
    std::swap(_fixpoints, tmp);
    
    std::cerr << _fixpoints.size() << " fixed points detected\n";
}

void check_cell(const nvis::ivec2& cell_id)
{
    _params.vectors.clear();
    _params.record = true;
    _params.lmin = 0.01;
    
    map_type pmap(*_field);
    pmap.precision(eps);
    try {
        std::vector<std::pair<int, int> > indices;
        process_cell<map_type>(indices, *_dataset, *_plane, cell_id, pmap, _index_dx, _params);
    } catch(rational_surface_found& e) {}
    _params.record = false;
    _params.lmin = 0;
}

void check_orbit(const nvis::vec2& x)
{
    map_type pmap(*_field);
    pmap.precision(eps);
    
    _orbits.push_back(std::pair<orbit_type, nvis::fvec3>());
    iterate<map_type>(_last, _orbits.back().first, pmap, niter);
    int period = best_period(_orbits.back().first, maxp, _plane_metric);
    _orbits.back().second = (*_cmap)(period);
    std::cerr << "best period at " << x << " is " << period << '\n';
    
    std::vector<nvis::vec2> steps;
    for (int i=0 ; i<_orbits.back().first.size() ; ++i) {
        steps.push_back(_orbits.back().first[i]);
    }
    
    double sf = safety_factor(steps, _plane_metric);
    std::cerr << "safety factor at " << x << " is " << sf << '\n';
    
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
        if (!logical) {
            x = (*_converter)(x);
            y = (*_converter)(y);
        }
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
                if (!logical) {
                    x = (*_converter)(x);
                }
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
            if (!logical) {
                x = (*_converter)(x);
            }
            glVertex2f(x[0], x[1]);
        }
    }
    glEnd();
    
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
            nvis::vec2 p0 = (*_plane)(id);
            const map_data_type& tv = (*_dataset)(id);
            const nvis::vec2& v = tv.first;
            if (nvis::norm(v) == 0) {
                continue;
            }
            int p = tv.second;
            nvis::vec2 p1 = p0 + v;
            nvis::fvec3 c = (*_cmap)(p);
            glColor3f(c[0], c[1], c[2]);
            if (!logical) {
                p0 = (*_converter)(p0);
                p1 = (*_converter)(p1);
            }
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
            const map_data_type& tv = _params.vectors[n].second;
            const nvis::vec2& v = tv.first;
            if (nvis::norm(v) == 0) {
                continue;
            }
            int p = tv.second;
            nvis::vec2 p1 = p0 + v;
            nvis::fvec3 c = (*_cmap)(p);
            glColor3f(c[0], c[1], c[2]);
            if (!logical) {
                p0 = (*_converter)(p0);
                p1 = (*_converter)(p1);
            }
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
        const fixpoint& fp = _fixpoints[i];
        if (fp.saddle) {
            glColor3f(1, 0, 0);
        } else {
            glColor3f(0, 0, 1);
        }
        glVertex3f(fp.pos[0], fp.pos[1], 0.5);
    }
    glEnd();
    // draw separatrices
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);
    glColor3f(1,1,1);
    glBegin(GL_LINES);
    for (int i=0 ; i<_fixpoints.size() ; ++i) {
        const fixpoint& fp = _fixpoints[i];
        if (fp.saddle) {
            nvis::vec2 x = fp.pos + 0.5*fp.evec[0];
            glVertex2f(x[0], x[1]);
            x -= fp.evec[0];
            glVertex2f(x[0], x[1]);
            x = fp.pos + 0.5*fp.evec[1];
            glVertex2f(x[0], x[1]);
            x -= fp.evec[1];
            glVertex2f(x[0], x[1]);
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
        if (!logical) {
            x = (*_converter)(x);
            y = (*_converter)(y);
        }
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
    map_type pmap(*_field);
    pmap.precision(eps);
    
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    
    if (button == GLUT_MIDDLE_BUTTON) {
        int i = floor(wc[0]);
        int j = floor(wc[1]);
        nvis::ivec2 cell_id(i,j);
        
        std::cerr << "computing index in cell " << cell_id << std::endl;
        
        // _params.verbose = true;
        _params.record = true;
        _params.lmin = _index_dx;
        _params.verbose = true;
        
        try {
            std::vector<std::pair<int, int> > indices;
            process_cell<map_type>(indices, *_dataset, *_plane, cell_id, pmap, _index_dx, _params);
            for (int i=0 ; i<indices.size() ; ++i) {
                const std::pair<int, int>& idx = indices[i];
                if (idx.first > 0) {
                    _center_cells.push_back(cell_id);
                    std::cerr << "cell contains a center\n";
                } else if (idx.first < 0) {
                    _saddle_cells.push_back(cell_id);
                    std::cerr << "cell contains a saddle\n";
                }
            }
            if (indices.empty()) {
                std::cerr << "cell has zero index\n";
            }
            
            _params.lmin = 0;
            
            for (int i=0 ; i<indices.size() ; ++i) {
                const std::pair<int, int>& idx = indices[i];
                if (idx.first != 0) {
                    nvis::vec2 x = _plane_bounds.min() + nvis::vec2((float)i + 0.5, (float)j+0.5)*_plane_spacing;
                    std::vector<fixpoint> chain;
                    bool found = linear_chain_analysis(pmap, _plane_metric, x, idx.second,
                                                       chain, 2*nvis::norm(_plane_spacing), 1.0e-3);
                    if (found) {
                        std::cerr << idx.second << "-chain of " << (chain[0].saddle ? "saddles" : "centers")
                                  << " found.\n";
                        std::copy(chain.begin(), chain.end(), std::back_inserter(_fixpoints));
                    } else {
                        std::cerr << "Newton search failed\n";
                    }
                }
            }
        } catch(...) {}
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
        map_type pmap(*_field);
        pmap.precision(eps);
        
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
    }
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    _last = nvis::vec2(-1,-1);
    
    initialize(argc, argv);
    save_to_file = false;
    
    if(strcmp(phys, "none")) {
        logical = false;
        _converter = new spurt::logical2physical(phys);
    } else {
        logical = true;
    }
    
    width = 1200;
    height = 800;
    
    // initialize GLUT
    /*    glutInit(&argc, argv);
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
    */
    std::vector<nvis::fvec3> colors;
    spurt::spiral_scale(colors, maxp, 1);
    std::vector<int> reference_values(maxp);
    for (int i=1 ; i<=maxp ; ++i) {
        reference_values[i] = i;
    }
    _cmap = new spurt::discrete_color_map<int>(reference_values, colors);
    init();
    /*    display();
    
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
    */
    return 0;
}
