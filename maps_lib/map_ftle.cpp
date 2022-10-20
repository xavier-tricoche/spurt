#include <iostream>
#include <iomanip>

#include <vector>
#include <complex>
#include <sstream>
#include <math.h>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>

// spurt
#include <math/math.hpp>

#include <iostream>
#include <list>

#include "definitions.hpp"
#include "xmt_poincare_map.hpp"
#include <data/grid.hpp>
#include <data/raster.hpp>

#include <util/wall_timer.hpp>

#include <image/nrrd_wrapper.hpp>


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;
using namespace map_analysis;
// using namespace div_cleaning;

double eps;
int resx, resy, maxp, maxit, it_step;
char* outs, *file, *ts;
double minx, maxx, miny, maxy;
double _minx, _maxx, _miny, _maxy;
double hx, hy;

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
    hestOptAdd(&hopt, "i",      "file",             airTypeString,  1, 1, &file,    NULL,       "input NRRD file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &outs,    NULL,       "output name");
    hestOptAdd(&hopt, "rx",     "resolution x",     airTypeInt,     0, 1, &resx,    "1024",     "sampling resolution in X");
    hestOptAdd(&hopt, "ry",     "resolution y",     airTypeInt,     0, 1, &resy,    "1024",     "sampling resolution in Y");
    hestOptAdd(&hopt, "maxi",   "max iterations",   airTypeInt,     0, 1, &maxit,   "10",       "max number of map iterations");
    hestOptAdd(&hopt, "s",      "iterations step",  airTypeInt,     0, 1, &it_step, "0",        "iteration step size");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  0, 1, &eps,     "1.0e-8",   "integration accuracy");
    hestOptAdd(&hopt, "minx",   "min x coord",      airTypeDouble,  0, 1, &minx,    "-10000",   "min x in bounding box");
    hestOptAdd(&hopt, "maxx",   "max x coord",      airTypeDouble,  0, 1, &maxx,    "10000",    "max x in bounding box");
    hestOptAdd(&hopt, "miny",   "min y coord",      airTypeDouble,  0, 1, &miny,    "-10000",   "min y in bounding box");
    hestOptAdd(&hopt, "maxy",   "max y coord",      airTypeDouble,  0, 1, &maxy,    "10000",    "max y in bounding box");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE value of discrete map after a given number of iterations. Intermediate steps are saved to disk if requested.",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline double lmax(int n, const std::vector< nvis::vec2 >& pos, const spurt::default_metric_type& metric)
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

    if (lmaj < 1 && !(n%1000)) {
        std::ostringstream os;
        os << "Jx = " << Jx << ", Jy = " << Jy << std::endl;
        os << "pos[" << n+1 << "] = " << pos[n+1] << ", pos[" << n-1 << "] = " << pos[n-1] << std::endl;
        std::cerr << os.str();
    }

    return lmaj;
}
using namespace spurt;
using namespace map_analysis;

namespace {
double _mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}
}

typedef image<nvis::vec3, 3, double, size_t>  field_type;
typedef spurt::grid::uniform_grid<double, 3>         grid_type;
typedef nvis::ivec3                           ivec_type;

struct wrapper {
    wrapper(const field_type& field) :
        _field(field), _periodic(false), _grid(field.grid()) {}

    wrapper(const wrapper& w) :
        _field(w._field), _periodic(w._periodic), _grid(_field.grid()) {}

    void periodic(bool p) {
        _periodic = p;
    }

    nvis::vec3 unproject(const nvis::vec2& v) const {
        return nvis::vec3(v[1], v[0], 0.0);
    }

    nvis::vec2 project(const nvis::vec3& v) const {
        return (_periodic ?
                nvis::vec2(_mod(v[1], _grid.resolution()[1] - 1),
                           _mod(v[0], _grid.resolution()[0] - 1)) :
                nvis::vec2(v[1], v[0]));
    }

    nvis::vec2 modulo(const nvis::vec2& x) const {
        return project(_grid.dmodulo(unproject(x)));
    }

    bool intersect(const nvis::vec3& p0, const nvis::vec3& p1) const {
        int z0 = (int)trunc(p0[2] / (_grid.resolution()[2] - 1));
        int z1 = (int)trunc(p1[2] / (_grid.resolution()[2] - 1));

        return z0 != z1;
    }

    nvis::vec3 operator()(const double&, const nvis::vec3& p, bool verbose = false) const {
        try {
            nvis::vec3 f = _field.value(p);
            return f;
        } catch (std::runtime_error& e) {
            throw;
        }
    }

    const grid_type& grid() const {
        return _grid;
    }

    const field_type&   _field;
    grid_type           _grid;
    bool                _periodic;
};

typedef xmt_poincare_map<wrapper>  map_type;

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
    nin = spurt::nrrd_utils::readNrrd(file);

    // verify data type
    if (nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }

    std::vector<double> _array;
    spurt::nrrd_utils::to_vector(_array, nin);
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    std::vector<nvis::vec3> vectors(_array.size() / 3);
    for (int i = 0 ; i < _array.size() / 3 ; ++i) {
        vectors[i][0] = _array[3*i  ];
        vectors[i][1] = _array[3*i+1];
        vectors[i][2] = _array[3*i+2];
    }
    field_type field(domain, vectors);

    map_type map(field);
    map.precision(eps);

    grid_type::bounds_type bbox(domain.bounds());
    nvis::bbox2 _bounds(nvis::vec2(bbox.min()[1], bbox.min()[0]),
                        nvis::vec2(bbox.max()[1], bbox.max()[0]));

    std::cerr << "bounding box = " << _bounds << '\n';

    bool per[2] = {true, false};
    spurt::default_metric_type metric2d(_bounds, per);

    _minx = _bounds.min()[0];
    _miny = _bounds.min()[1];
    _maxx = _bounds.max()[0];
    _maxy = _bounds.max()[1];

    minx = std::max(_minx, minx);
    maxx = std::min(_maxx, maxx);
    miny = std::max(_miny, miny);
    maxy = std::min(_maxy, maxy);

    assert(minx < maxx && miny < maxy);

    hx = (maxx - minx) / (double)resx;
    hy = (maxy - miny) / (double)resy;

    unsigned int blocked = 0;
    float last_pct = 0.;

    std::cout << "minx = " << minx << ", miny = " << miny << ", maxx = " << maxx << ", maxy = " << maxy
              << ", hx = " << hx << ", hy = " << hy << std::endl;
    std::cout << "it_step = " << it_step << std::endl;

    nvis::timer timer;

    std::cout << "initializing coordinates\n";
    std::vector<nvis::vec2> pos_fwd(resx*resy), pos_bwd(resx*resy);
    for (int n = 0 ; n < resx*resy ; ++n) {
        unsigned int i = n % resx;
        unsigned int j = n / resx;
        nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);
        pos_fwd[n] = pos_bwd[n] = x;
    }


    for (int iteration = it_step ; iteration <= maxit ; iteration += it_step) {
        std::cout << "Computing flow map..." << std::endl;
        std::cout << "\niteration " << iteration << " from " << maxit << '\n';

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

                const map_type* clone = map.clone();

                nvis::vec2 prev = pos_fwd[n];

                try {
                    pos_fwd[n] = map.field().modulo(clone->map(pos_fwd[n], it_step));
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

                const map_type* clone = map.clone();

                unsigned int i = n % resx;
                unsigned int j = n / resx;
                nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);

                try {
                    pos_bwd[n] = map.field().modulo(clone->map(pos_bwd[n], -it_step));
                } catch (...) {
                    // std::cerr << "exception caught in backward integration\n";
                }
            }
        }
        std::cerr << '\n';

        std::cout << "computing FTLE...\n";
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

                double lm = lmax(n, pos_fwd, metric2d);
                // if (lm < 1) {
                //  std::ostringstream os;
                //  os << "lmax at #" << n << " was negative (" << lm << ")" << std::endl;
                //  std::cerr << os.str();
                // }

                lmax_f = log(std::max(1., lm));
                ftle_f = lmax_f;

                // if (!(n%1000)) {
                //  std::ostringstream os;
                //  os << "pos_fwd[" << n << "] = " << pos_fwd[n] << std::endl;
                //  std::cerr << os.str();
                // }

                lm = lmax(n, pos_fwd, metric2d);
                // if (lm < 1) {
                //  std::ostringstream os;
                //  os << "lmin at #" << n << " was negative (" << lm << ")" << std::endl;
                //  std::cerr << os.str();
                // }
                lmax_b = log(std::max(1., lm));
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
        os << outs << "-ftle-p=" << iteration << "_from_" << maxit
           << "-res=" << resx << "x" << resy
           << "-eps=" << eps << ".nrrd";
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
    }

    return 0;
}
