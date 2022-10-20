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
#include <misc/progress.hpp>

#include <util/wall_timer.hpp>

#include <image/nrrd_wrapper.hpp>

#include <Eigen/Eigenvalues>

#ifndef NO_TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>
#endif

#ifndef NO_TBB
tbb::mutex output_mutex;
tbb::mutex progress_mutex;
tbb::mutex orbit_mutex;
tbb::atomic<size_t> progress_counter;
#else
size_t progress_counter;
#endif

void write_to_ostream(std::ostream& os, const std::string& str) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(output_mutex);
    #endif
        os << str << '\n';
    }
}

void update_progress(spurt::ProgressDisplay& progress) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(progress_mutex);
    #endif
        progress.update(progress_counter);
    }
}

using namespace nvis;
using namespace map_analysis;

double eps;
unsigned int maxp, maxit, it_step, init_step, nparsed_r, nparsed_x, nparsed_y, nparsed_z, nparsed_w;
char* outs, *file, *ts;
double minx, maxx, miny, maxy, minz, maxz, minw, maxw, k1, k2;
double _minx, _maxx, _miny, _maxy;
double hx, hy;
double xb[2], yb[2], zb[2], wb[2];
size_t nx, ny, nz, nw;

typedef spurt::symplectic4D::state_type state_type;
typedef spurt::symplectic4D::deriv_type deriv_type;
typedef spurt::symplectic4D::bounds_type bounds_type;

typedef Eigen::Matrix<size_t, 4, 1> int4;
typedef Eigen::Matrix<double, 4, 4> matrix_type;

typedef spurt::map_metric<4, state_type, bounds_type> metric4D_type;

int4 res;
state_type spacing;
size_t nsamples;

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
    hestOptAdd(&hopt, "k1",   "K1",                 airTypeDouble, 0, 1, &k1,       "0.5",              "K1 parameter");
    hestOptAdd(&hopt, "k2",   "K2",                 airTypeDouble, 0, 1, &k2,       "0.5",              "K2 parameter");
    hestOptAdd(&hopt, "o",    "output",             airTypeString, 1, 1, &outs,      NULL,              "output name");
    // hestOptAdd(&hopt, "r",    "resolution",         airTypeUInt,   0, 4, &res4,      "256 256 256 256", "sampling resolution", &nparsed_r);
    hestOptAdd(&hopt, "maxi", "max iterations",     airTypeUInt,   0, 1, &maxit,     "10",              "max number of map iterations");
    hestOptAdd(&hopt, "s",    "iterations step",    airTypeUInt,   0, 1, &it_step,   "0",               "iteration step size");
    hestOptAdd(&hopt, "nx",    "nb samples",    airTypeULongInt,   0, 1, &nx,   "64",               "number of samples along x axis");
    hestOptAdd(&hopt, "ny",    "nb samples",    airTypeULongInt,   0, 1, &ny,   "64",               "number of samples along y axis");
    hestOptAdd(&hopt, "nz",    "nb samples",    airTypeULongInt,   0, 1, &nz,   "64",               "number of samples along z axis");
    hestOptAdd(&hopt, "nw",    "nb samples",    airTypeULongInt,   0, 1, &nw,   "64",               "number of samples along w axis");
    hestOptAdd(&hopt, "i",    "initial iterations", airTypeUInt,   0, 1, &init_step, "0",               "preliminary number of iterations");
    hestOptAdd(&hopt, "xb",   "x bounds",           airTypeDouble, 0, 2, &xb,        "-0.5 0.5",        "bounds in x coordinate", &nparsed_x);
    hestOptAdd(&hopt, "yb",   "y bounds",           airTypeDouble, 0, 2, &yb,        "-0.5 0.5",        "bounds in y coordinate", &nparsed_y);
    hestOptAdd(&hopt, "zb",   "z bounds",           airTypeDouble, 0, 2, &zb,        "-0.5 0.5",        "bounds in z coordinate", &nparsed_z);
    hestOptAdd(&hopt, "wb",   "w bounds",           airTypeDouble, 0, 2, &wb,        "-0.5 0.5",        "bounds in w coordinate", &nparsed_w);

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE value of 4D standard map after a given number of iterations. Intermediate steps can be saved to disk.",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline double lmax(size_t n, const std::vector< state_type >& pos, const metric4D_type& metric)
{
    size_t tmp(n);
    int4 id;
    for (size_t d=0; d<4; ++d) {
        id[d] = tmp % res[d];
        if (id[d] == 0 || id[d] == res[d]-1) { return -1; }
        tmp /= res[d];
    }

    // look for valid neighboring values to compute derivatives
    size_t stride=1;
    matrix_type J = matrix_type::Zero();
    for (size_t d=0; d<4; ++d) {
        state_type der = 1./(2.*spacing[d]) * metric.displacement(pos[n-stride], pos[n+stride]);
        for (size_t c=0; c<4; ++c) J(d,c) = der[c];
    }

    matrix_type C = J.transpose() * J;
    Eigen::EigenSolver<matrix_type> solver(C);
    auto eigenvals = solver.eigenvalues();
    double lambdas[4] = { eigenvals[0].real(), eigenvals[1].real(), eigenvals[2].real(), eigenvals[3].real() };
    return *std::max_element(lambdas, lambdas+4);
}
using namespace spurt;
using namespace map_analysis;

namespace {
double __mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}
}

typedef spurt::symplectic4D    map_type;

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    std::cout << "parameters: max period = " << maxp
              << ", max iterations = " << maxit << ", resolution = " << '\n' << res << '\n'
              << ", eps = " << eps << '\n';

    map_type map(k1, k2, 1.0);

    bounds_type bounds = map.bounds();

    std::cerr << "bounding box =\n" << bounds.min()  << "->\n" << bounds.max() << '\n';

    std::array<bool, 4> per({true, true, true, true});
    metric4D_type metric4d(bounds, per);

    std::cout << "sampling resolution is:\n" << res << '\n';
    std::cout << "there were " << nparsed_r << " values parsed\n";

    res[0] = nx;
    res[1] = ny;
    res[2] = nz;
    res[3] = nw;

    spacing = (bounds.max() - bounds.min()).array() / (state_type(res[0], res[1], res[2], res[3]) - state_type(1, 1, 1, 1)).array();

    unsigned int blocked = 0;
    float last_pct = 0.;

    if (it_step == 0) it_step = maxit;

    std::cout << "bounds: \n" << bounds.min() << "\n->\n" << bounds.max() << '\n';
    std::cout << "sampling distance:\n" << spacing << '\n';
    std::cout << "it_step = " << it_step << std::endl;
    std::cout << "size of double = " << sizeof(double) << std::endl;

    nvis::timer timer;
    nvis::timer subtimer;

    nsamples = res[0]*res[1]*res[2]*res[3];

    std::cout << "initializing coordinates\n";
    std::vector<state_type> pos_fwd(nsamples), pos_bwd(nsamples);

    spurt::ProgressDisplay progress;

    progress.start(nsamples, "Initializing flow map");
    progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
                      [&](tbb::blocked_range<size_t> r)
    {
        for (size_t n=r.begin() ; n!=r.end() ; ++n) {
            size_t tmp = n;
            state_type id;
            for (int i=0; i<4; ++i) {
                id[i] = tmp % res[i];
                tmp /= res[i];
            }
            state_type x = bounds.min().array() + id.array() * spacing.array();
            pos_fwd[n] = pos_bwd[n] = x;
            ++progress_counter;
            update_progress(progress);
        }
    });
    progress.end();

    timer.restart();

    int to=(init_step>0 ? init_step : it_step);

    for (int at=0; to<=maxit; to+=it_step) {
        std::cout << "\niteration " << to << " from " << maxit << '\n';

        subtimer.restart();

        // do that in parallel
        progress.start(nsamples, "Integrating forward flow map");
        progress_counter = 0;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
                          [&](tbb::blocked_range<size_t> r)
        {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                try {
                    // std::ostringstream os;
                    // os << "pos_fwd[" << n << "] before" << pos_fwd[n] << "->";
                    pos_fwd[n] = metric4d.modulo(map.map(pos_fwd[n], to-at));
                    // os << "after: " << pos_fwd[n] << '\n';
                    // std::cout << os.str();
                } catch (...) {}
                ++progress_counter;
                update_progress(progress);
            }
        });
        progress.end();

        progress.start(nsamples, "Integrating backward flow map");
        progress_counter = 0;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples), [&](tbb::blocked_range<size_t> r)
        {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                try {
                    pos_bwd[n] = metric4d.modulo(map.map(pos_bwd[n], at-to));
                } catch (...) {}
                ++progress_counter;
                update_progress(progress);
            }
        });
        progress.end();

        std::cout << "computing FTLE...\n";
        float* _ftle = (float*)calloc(nsamples*2, sizeof(float));

        progress.start(nsamples, "Computing FTLE");
        progress_counter = 0;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
                          [&](tbb::blocked_range<size_t> r)
        {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                double lmax_f = 0., lmax_b = 0.;
                double ftle_f = 0., ftle_b = 0.;

                double lm = lmax(n, pos_fwd, metric4d);
                lmax_f = log(std::max(1., lm));
                ftle_f = lmax_f;

                lm = lmax(n, pos_bwd, metric4d);
                lmax_b = log(std::max(1., lm));
                ftle_b = lmax_b;

                if (!std::isinf(ftle_f) && !std::isnan(ftle_f)) {
                    _ftle[2*n] = ftle_f;
                }
                if (!std::isinf(ftle_b) && !std::isnan(ftle_b)) {
                    // std::cerr << ftle_b << " in backward\n";
                    _ftle[2*n+1] = ftle_b;
                }
                ++progress_counter;
                update_progress(progress);
            }
        });
        progress.end();

        std::cout << "Last iteration took " << subtimer.elapsed() << " s.\n";
        std::cout << "Total computation time so far is " << timer.elapsed()
                  << " s.\n";

        std::ostringstream os;
        int w=std::to_string(maxit).size();
        os << outs << "-ftle-p=" << std::setw(w) << std::setfill('0') << to << "_from_" << maxit
           << "-K1=" << k1  << "-K2=" << k2 << "-res=" << res[0] << "x" << res[1] << "x" << res[2] << "x" << res[3]
           << ".nrrd";
        Nrrd* nout = nrrdNew();
        size_t sizes[5];
        sizes[0] = 2;
        for (int i=0; i<4; ++i) sizes[i+1] = res[i];
        if (nrrdWrap_nva(nout, _ftle, nrrdTypeFloat, 5, sizes)) {
            std::cerr << "ERROR while wrapping data: " << biffGetDone(NRRD)
                      << std::endl;
            if (_ftle) {
                delete[] _ftle;
            } else {
                nrrdNuke(nout);
            }
            exit(-1);
        }
        int kinds[5] = {nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace, nrrdKindSpace, nrrdKindSpace};
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoKind, kinds);
        int centers[5] = {nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell, nrrdCenterCell, nrrdCenterCell};
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, centers);
        double spc[5] = {AIR_NAN, spacing[0], spacing[1], spacing[2], spacing[3]};
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        if (nrrdSave(os.str().c_str(), nout, NULL)) {
            std::cerr << "ERROR while exporting file: " << biffGetDone(NRRD)
                      << std::endl;
            exit(-1);
        }
        nrrdNuke(nout);
        std::cout << "exported " << os.str() << std::endl;
        at=to;
    }

    return 0;
}
