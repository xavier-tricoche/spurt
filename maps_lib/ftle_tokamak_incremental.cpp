#include <iostream>
#include <iomanip>

#include <vector>
#include <complex>
#include <sstream>
#include <math.h>

#include <boost/format.hpp>
#include <boost/limits.hpp>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>

// xavier
#include <math/math.hpp>

// christoph
#include <tokamak/poincare_map.hpp>
#include <tokamak/tokamak_nimrod.hpp>
#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <tokamak/tokamak_nrrd.hpp>

#include "definitions.hpp"

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;
using namespace map_analysis;

double eps;
int resx, resy, maxp, maxit, di;
char* outs, *file, *ts;
double minx, maxx, miny, maxy;
double _minx, _maxx, _miny, _maxy;
double hx, hy;
int square;

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
    hestOptAdd(&hopt, "f",      "file",             airTypeString,  1, 1, &file,    NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",             airTypeString,  1, 1, &ts,      NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output name");
    hestOptAdd(&hopt, "rx",     "resolution x",     airTypeInt,     0, 1, &resx,    "1024",     "sampling resolution in X");
    hestOptAdd(&hopt, "ry",     "resolution y",     airTypeInt,     0, 1, &resy,    "1024",     "sampling resolution in Y");
    hestOptAdd(&hopt, "maxi",   "max iterations",   airTypeInt,     0, 1, &maxit,   "10",       "max number of map iterations");
    hestOptAdd(&hopt, "di",     "delta iterations", airTypeInt,     0, 1, &di,      "0",        "number of map iterations between exports");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  0, 1, &eps,     "1.0e-6",   "integration accuracy");
    hestOptAdd(&hopt, "minx",   "min x coord",      airTypeDouble,  0, 1, &minx,    "-10000",   "min x in bounding box");
    hestOptAdd(&hopt, "maxx",   "max x coord",      airTypeDouble,  0, 1, &maxx,    "10000",    "max x in bounding box");
    hestOptAdd(&hopt, "miny",   "min y coord",      airTypeDouble,  0, 1, &miny,    "-10000",   "min y in bounding box");
    hestOptAdd(&hopt, "maxy",   "max y coord",      airTypeDouble,  0, 1, &maxy,    "10000",    "max y in bounding box");
    hestOptAdd(&hopt, "type",   "data type",        airTypeInt,     0, 1, &square,  "1",        "dataset type: 0: physical tokamak, 1: parametric tokamak, 2: raster");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE value of discrete map after a given number of iterations. Intermediate steps are not saved to disk.",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline int pos(const nvis::vec2& x, unsigned int res)
{
    using namespace static_data;
    int i = floor(res * (x[0] - metric.bounds().min()[0]) / metric.width());
    int j = floor(res * (x[1] - metric.bounds().min()[1]) / metric.height());
    return i + j*res;
}

inline double lmax(int n, const std::vector< nvis::vec2 >& pos, const xavier::map_metric& metric)
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
              
    tokamak_field* field;
    
    switch (square) {
        case 0:
            field = new tokamak_nimrod(std::string(file), std::string(ts));
            break;
        case 1:
            field = new tokamak_nimrod_parametric(std::string(file), std::string(ts));
            break;
        case 2:
            field = new tokamak_nrrd(std::string(file));
            break;
        default:
            std::cerr << "unknown data type\n";
            exit(-1);
    }
    
    poincare_map _map(field);
    _map.precision(eps);
    
    bool per[2] = {true, false};
    xavier::map_metric metric(field->bounds(), per);
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
    for (int n = 0 ; n < resx*resy ; ++n) {
        unsigned int i = n % resx;
        unsigned int j = n / resx;
        nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);
        pos_fwd[n] = pos_bwd[n] = x;
    }
    
    unsigned int blocked = 0;
    float last_pct = 0.;
    
    std::cout << "minx = " << minx << ", miny = " << miny << ", maxx = " << maxx << ", maxy = " << maxy
              << ", hx = " << hx << ", hy = " << hy << std::endl;
              
    nvis::timer timer;
    
    std::cout << "Computing flow map..." << std::endl;
    
    if (!di) {
        di = maxit;
    }
    
    for (int iter = di ; iter <= maxit ; iter += di) {
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
                
                const poincare_map* clone = _map.clone();
                
                try {
                    pos_fwd[n] = metric.modulo(clone->map(pos_fwd[n], di));
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
                
                const poincare_map* clone = _map.clone();
                
                try {
                    pos_bwd[n] = metric.modulo(clone->map(pos_bwd[n], -di));
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
        os << outs << "-ftle-p=" << iter << "_from_" << maxit << "-t=" << ts << ".nrrd";
        Nrrd* nout = nrrdNew();
        size_t dims[] = {3, resx, resy};
        if (nrrdWrap_nva(nout, _ftle, nrrdTypeFloat, 3, dims)) {
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
    }
    
    return 0;
}























