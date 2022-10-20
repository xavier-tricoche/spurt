#include <iostream>
#include <iomanip>
#include <list>
#include <map>
#include <vector>
#include <complex>
#include <sstream>
#include <math.h>

#include <boost/format.hpp>
#include <boost/limits.hpp>

// teem
#include <teem/hest_helper.hpp>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>

// spurt
#include <math/math.hpp>
#include "definitions.hpp"
#include "xmt_poincare_map.hpp"
#include <data/grid.hpp>
#include <data/raster_data.hpp>
#include <maps/period.hpp>
#include "map_field_wrapper.hpp"
#include <image/nrrd_wrapper.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;
using namespace map_analysis;

double eps;
int resx, resy, maxp, maxit;
char* outs, *file, *ts;
double minx, maxx, miny, maxy;
double _minx, _maxx, _miny, _maxy;
double hx, hy;
bool square;

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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &file,    NULL,       "input NRRD file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output name");
    hestOptAdd(&hopt, "rx",     "resolution x",     airTypeInt,     0, 1, &resx,    "1024",     "sampling resolution in X");
    hestOptAdd(&hopt, "ry",     "resolution y",     airTypeInt,     0, 1, &resy,    "1024",     "sampling resolution in Y");
    hestOptAdd(&hopt, "maxi",   "max iterations",   airTypeInt,     0, 1, &maxit,   "10",       "max number of map iterations");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  0, 1, &eps,     "1.0e-6",   "integration accuracy");
    hestOptAdd(&hopt, "minx",   "min x coord",      airTypeDouble,  0, 1, &minx,    "-10000",   "min x in bounding box");
    hestOptAdd(&hopt, "maxx",   "max x coord",      airTypeDouble,  0, 1, &maxx,    "10000",    "max x in bounding box");
    hestOptAdd(&hopt, "miny",   "min y coord",      airTypeDouble,  0, 1, &miny,    "-10000",   "min y in bounding box");
    hestOptAdd(&hopt, "maxy",   "max y coord",      airTypeDouble,  0, 1, &maxy,    "10000",    "max y in bounding box");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Compute FTLE value of discrete map after a given number of iterations. Intermediate steps are not saved to disk.",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

namespace {
double __mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}
}

typedef grid<double, 3>                         grid_type;
typedef raster_data<nvis::vec3, double, 3>      field_type;


typedef xmt_poincare_map<spurt::map::wrapper<field_type> > map_type;

inline int pos(const nvis::vec2& x, unsigned int res, const spurt::map_metric& metric)
{
    using namespace static_data;
    int i = floor(res * (x[0] - metric.bounds().min()[0]) / metric.width());
    int j = floor(res * (x[1] - metric.bounds().min()[1]) / metric.height());
    return i + j*res;
}

inline double lmax(int n, const std::vector< nvis::vec2 >& pos, const spurt::map_metric& metric)
{
    static double hx = (maxx - minx) / (double)resx;
    static double hy = (maxy - miny) / (double)resy;
    
    unsigned int i = n % resx;
    unsigned int j = n / resx;
    if (i == 0 || i == resx - 1 || j == 0 || j == resy - 1) {
        return -1.0;
    }
    
    // look for valid neighboring values to compute derivatives
    nvis::vec2 Jx = 1. / (2.*hx) * metric.displacement(pos[n-1], pos[n+1]);
    nvis::vec2 Jy = 1. / (2.*hy) * metric.displacement(pos[n-resx], pos[n+resx]);
    double a = nvis::inner(Jx, Jx);
    double b = nvis::inner(Jx, Jy);
    double c = nvis::inner(Jy, Jy);
    
    double lmaj = 0.5 * (a * a + 2 * b * b + c * c + (a + c) * sqrt((a - c) * (a - c) + 4 * b * b));
    
    return lmaj;
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    std::cout << "parameters: max period = " << maxp
              << ", max iterations = " << maxit << ", resolution = " << resx << " x " << resy
              << ", eps = " << eps << '\n';
              
    Nrrd* nin = nrrdNew();
    nin = spurt::readNrrd(file);
    
    // verify data type
    if (nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> __array;
    spurt::to_vector(__array, nin);
    nvis::ivec3 dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    std::vector<nvis::vec3> vectors(__array.size() / 3);
    for (int i = 0 ; i < __array.size() / 3 ; ++i) {
        vectors[i][0] = __array[3*i  ];
        vectors[i][1] = __array[3*i+1];
        vectors[i][2] = __array[3*i+2];
    }
    field_type __field(domain, vectors);
    __field.verbose(false);
    
    double h = eps;
    
    map_type __map(__field);
    __map.precision(eps);
    
    grid_type::bounds_type bbox(domain.bounds());
    nvis::bbox2 __bounds = nvis::bounding_box<nvis::vec2> (nvis::vec2(bbox.min()[1], bbox.min()[0]),
                           nvis::vec2(bbox.max()[1], bbox.max()[0]));
                           
    bool per[2] = {true, false};
    spurt::map_metric metric(__bounds, per);
    _minx = metric.bounds().min()[0];
    _miny = metric.bounds().min()[1];
    _maxx = metric.bounds().max()[0];
    _maxy = metric.bounds().max()[1];
    
    minx = std::max(_minx, minx);
    maxx = std::min(_maxx, maxx);
    miny = std::max(_miny, miny);
    maxy = std::min(_maxy, maxy);
    
    assert(minx < maxx && miny < maxy);
    
    hx = (maxx - minx) / (double)resx;
    hy = (maxy - miny) / (double)resy;
    
    std::vector< nvis::vec2 > pos_fwd(resx*resy);
    std::vector< nvis::vec2 > pos_bwd(resx*resy);
    
    unsigned int blocked = 0;
    float last_pct = 0.;
    
    std::cout << "minx = " << minx << ", miny = " << miny << ", maxx = " << maxx << ", maxy = " << maxy
              << ", hx = " << hx << ", hy = " << hy << std::endl;
              
    nvis::timer timer;
    
    std::cout << "Computing flow map..." << std::endl;
    
    // do that in parallel
    #pragma omp parallel
    {
        int thread_id = 0;
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < resx*resy ; ++n) {
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            if (thread_id == 0) {
                float pct = 100 * n / (resx * resy);
                std::cout << (int)pct << "% forward integration completed      \r";
            }
            
            const map_type* clone = __map.clone();
            
            unsigned int i = n % resx;
            unsigned int j = n / resx;
            nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);
            
            try {
                pos_fwd[n] = metric.modulo(clone->map(x, maxit));
            } catch (...) {
            }
        }
    }
    std::cerr << '\n';
    
    #pragma omp parallel
    {
        int thread_id = 0;
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < resx*resy ; ++n) {
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            if (thread_id == 0) {
                float pct = 100 * n / (resx * resy);
                std::cout << (int)pct << "% backward integration completed      \r";
            }
            
            const map_type* clone = __map.clone();
            
            unsigned int i = n % resx;
            unsigned int j = n / resx;
            nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);
            
            try {
                pos_bwd[n] = metric.modulo(clone->map(x, -maxit));
                // std::cerr << "distance to start = " << metric.distance(x, pos_bwd[n])
                // << '\n';
            } catch (...) {
                // std::cerr << "exception caught in backward integration\n";
            }
        }
    }
    std::cerr << '\n';
    
    
    float* _ftle = (float*)calloc(3 * resx * resy, sizeof(float));
    float* _lmax = (float*)calloc(3 * resx * resy, sizeof(float));
    
    #pragma omp parallel
    {
        int thread_id = 0;
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < resx*resy ; ++n) {
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            double lmax_f = 0., lmax_b = 0.;
            double ftle_f = 0., ftle_b = 0.;
            lmax_f = log(std::max(1., lmax(n, pos_fwd, metric)));
            ftle_f = lmax_f;
            lmax_b = log(std::max(1., lmax(n, pos_bwd, metric)));
            ftle_b = lmax_b;
            
            if (!std::isinf(ftle_f) && !std::isnan(ftle_f)) {
                _ftle[3*n] = ftle_f;
            }
            if (!std::isinf(ftle_b) && !std::isnan(ftle_b)) {
                // std::cerr << ftle_b << " in backward\n";
                _ftle[3*n+2] = ftle_b;
            }
        }
    }
    
    std::cerr << timer.elapsed() << " seconds needed for the computation\n";
    
    std::ostringstream os;
    os << outs << "-ftle-p=" << maxit << ".nrrd";
    Nrrd* nout = nrrdNew();
    if (nrrdWrap_va(nout, _ftle, nrrdTypeFloat, 3, 3, resx, resy)) {
        std::cout << "ERROR while wrapping data: " << biffGetDone(NRRD)
                  << std::endl;
        if (_ftle) {
            delete[] _ftle;
        } else {
            nrrdNuke(nout);
        }
        exit(-1);
    }
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoKind, nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell);
    nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, airNaN(), hx, hy);
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                  << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    std::cout << "exported " << os.str() << std::endl;
    
    return 0;
}




















