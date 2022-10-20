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
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/progress.hpp>
#include <misc/option_parse.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>

typedef Eigen::Matrix<double, 2, 2> matrix_t;
typedef Eigen::Matrix<double, 2, 1> vector_t;
typedef Eigen::Matrix<int, 2, 1> ivector_t;
typedef Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver_t;
typedef spurt::raster_grid<2> raster_t;
typedef spurt::map_metric<2> metric_t;

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;
using namespace map_analysis;
// using namespace div_cleaning;

std::array<size_t, 2> res({1024, 1024});
nvis::bbox2 bounds;
size_t maxit=100, it_step=1, init_step=0;
std::string out_basename;
vector_t minp, maxp, h;
double K;
bool add_ridge=false;
bool add_ftle=false;
bool verbose, full_domain=true;

tbb::mutex progress_mutex;
tbb::atomic<size_t> progress_counter;

typedef tbb::blocked_range<int> tbb_block;

void update_progress(spurt::ProgressDisplay& progress) {
    {
        tbb::mutex::scoped_lock lock(progress_mutex);
        progress.update(progress_counter);
    }
}

spurt::ProgressDisplay partial_progress(true), total_progress(false);

void initialize(int argc, const char* argv[]) {
    namespace xcl = spurt::command_line;

    xcl::option_traits
        required(true, false, "Required Options"),
    optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
    "Compute standard map's FTLE and FTLE's ridge strength");

    maxit = 100;
    init_step = 0;
    it_step = 1;
    std::array<double, 4> bounds4({0, 0, 1, 1});
    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("k", K, "K parameter", required);
        parser.add_value("output", out_basename, "Output basename", required);
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional);
        parser.add_value("maxi", maxit, maxit, "Total nb. of iterations", optional);
        parser.add_value("step", it_step, it_step, "Nb. of iterations between outputs", optional);
        parser.add_value("init", init_step, init_step, "Nb. of iterations before 1st output", optional);
        parser.add_tuple<4>("bounds", bounds4, bounds4, "Sampling bounds (within unit square)", optional);
        parser.add_value("addr", add_ridge, add_ridge, "Accumulate ridge strength results", optional);
        parser.add_value("addf", add_ftle, add_ftle, "Accumulate ftle results", optional);
        parser.add_flag("verbose", verbose, "Toggle verbose output", optional);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
            << e.what() << "\n"
                << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    bounds.min() = nvis::vec2(bounds4[0], bounds4[1]);
    bounds.max() = nvis::vec2(bounds4[2], bounds4[3]);
}

inline size_t ij2n(size_t i, size_t j) {
    return i + res[0]*j;
}

inline double lmax(size_t n, const std::vector< nvis::vec2 >& pos, const spurt::map_metric<2>& metric)
{
    size_t i = n % res[0];
    size_t j = n / res[1];
    nvis::vec2 Jx, Jy;
    if (full_domain) {
        size_t i_plus, j_plus, i_minus, j_minus;
        i_plus = (i+1)%res[0];
        i_minus = (i+res[0]-1)%res[0];
        j_plus = (j+1)%res[1];
        j_minus = (j+res[1]-1)%res[1];
        Jx = 1. / (2.*h[0]) * metric.displacement(pos[ij2n(i_minus,j)], pos[ij2n(i_plus,j)]);
        Jy = 1. / (2.*h[1]) * metric.displacement(pos[ij2n(i,j_minus)], pos[ij2n(i,j_plus)]);
    }
    else {
        if (i==0) Jx = 1./h[0] * metric.displacement(pos[n], pos[n+1]);
        else if (i==res[0]-1) Jx = 1./h[0] * metric.displacement(pos[n-1], pos[n]);
        else Jx = 1./(2*h[0]) * metric.displacement(pos[n-1], pos[n+1]);
        if (j==0) Jy = 1./h[1] * metric.displacement(pos[n], pos[n+res[0]]);
        else if (j==res[1]-1) Jy = 1./h[1] * metric.displacement(pos[n-res[0]], pos[n]);
        else Jy = 1./(2*h[1]) * metric.displacement(pos[n-res[0]], pos[n+res[0]]);
    }
    double a = nvis::inner(Jx, Jx);
    double b = nvis::inner(Jx, Jy);
    double c = nvis::inner(Jy, Jy);

    double lmaj = 0.5 * (a*a + 2*b*b + c*c + (a+c)*sqrt((a-c)*(a-c) + 4*b*b));

    return lmaj;
}

inline double hlmin(size_t n, const std::vector<double>& values)
{
    size_t i = n % res[0];
    size_t j = n / res[0];
    double hxx, hxy, hyy;
    size_t i_plus, j_plus, i_minus, j_minus;
    if (full_domain || (i>0 && i<res[0]-1 && j>0 && j<res[1]-1)) {
        i_plus = (i+1)%res[0];
        i_minus = (i+res[0]-1)%res[0];
        j_plus = (j+1)%res[1];
        j_minus = (j+res[1]-1)%res[1];
    }
    else return 0;

    hxx = 1./(h[0]*h[0]) * (values[ij2n(i_minus,j)] - 2*values[n] + values[ij2n(i_plus,j)]);
    hyy = 1./(h[1]*h[1]) * (values[ij2n(i,j_minus)] - 2*values[n] + values[ij2n(i,j_plus)]);
    hxy = 1./(4.*h[0]*h[1]) * (values[ij2n(i_minus,j_minus)] +
                               values[ij2n(i_plus,j_plus)] -
                               values[ij2n(i_plus,j_minus)] - values[ij2n(i_minus,j_plus)]);

#if 0
    matrix_t H;
    H(0,0) = hxx;
    H(1,1) = hyy;
    H(0,1) = H(1,0) = hxy;
    eigensolver_t solver;
    solver.compute(H);
    return solver.eigenvalues()[0]; // min eigenvalue of the Hessian
#else
    double tr = hxx + hyy;
    double det = hxx*hyy - hxy*hxy;
    return (tr - sqrt(tr*tr - 4*det))/2; // min eigenvalue of the Hessian
#endif
}

using namespace spurt;
using namespace map_analysis;

namespace {
double __mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}
}

typedef spurt::standard_map map_type;

int main(int argc, const char* argv[])
{
    initialize(argc, argv);
    std::cout << "parameters: max iterations = " << maxit
              << ", resolution = " << res[0] << " x " << res[1] << '\n';

    map_type map(K);

    nvis::bbox2 _bounds(nvis::vec2(0,0), nvis::vec2(1,1));
    bool per[2] = {true, true};
    metric_t metric2d(_bounds, per);
    bounds.min()[0] = std::max(static_cast<double>(0), bounds.min()[0]);
    bounds.min()[1] = std::max(static_cast<double>(0), bounds.min()[1]);
    bounds.max()[0] = std::min(static_cast<double>(1), bounds.max()[0]);
    bounds.max()[1] = std::min(static_cast<double>(1), bounds.max()[1]);
    assert(nvis::all(bounds.max()>bounds.min()));
    full_domain = (bounds.min()[0] == 0) && (bounds.min()[1] == 0) &&
        (bounds.size()[0] == 1) && (bounds.size()[1] == 1);
    h[0] = bounds.size()[0]/double(res[0]);
    h[1] = bounds.size()[1]/double(res[1]);

    size_t npoints = res[0]*res[1];

    if (it_step == 0) it_step = maxit;

    std::cout << "Number of samples = " << npoints << '\n';
    std::cout << "selected bounding box = "
        << bounds << '\n';
    std::cout << "h[0] = " << h[0] << ", h[1] = " << h[1] << std::endl;
    std::cout << "it_step = " << it_step << std::endl;

    std::cout << "Initializing coordinates\n";
    std::vector<nvis::vec2> pos_fwd(npoints), pos_bwd(npoints);
    for (size_t n = 0 ; n < npoints ; ++n) {
        size_t i = n % res[0];
        size_t j = n / res[0];
        nvis::vec2 x(bounds.min()[0] + h[0]*(double)i, bounds.min()[1] + h[1]*(double)j);
        pos_fwd[n] = pos_bwd[n] = x;
    }

    std::vector<double> ftle_fwd(npoints, 0);
    std::vector<double> ftle_bwd(npoints, 0);
    std::vector<double> ridge_strength_fwd(npoints, 0);
    std::vector<double> ridge_strength_bwd(npoints, 0);

    size_t to=(init_step>0 ? init_step : it_step);

    total_progress.start();

    for (size_t at=0; to<=maxit; to+=it_step) {
        std::cout << "\niteration " << to << " from " << maxit << '\n';
        partial_progress.start(npoints, "Forward stdmap");
        progress_counter = 0;

        // do that in parallel
        tbb::parallel_for(tbb_block(0,npoints), [&](tbb_block r) {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                ++progress_counter;
                update_progress(partial_progress);
                const map_type* clone = map.clone();
                nvis::vec2 prev = pos_fwd[n];
                try {
                    pos_fwd[n] = metric2d.modulo(clone->map(pos_fwd[n], to-at));
                } catch (...) {
                }

                delete clone;
            }
        });
        partial_progress.stop();
        std::cerr << '\n';
        std::cout << partial_progress << '\n';

        progress_counter = 0;

        partial_progress.start(npoints, "Backward stdmap");
        tbb::parallel_for(tbb_block(0,npoints), [&](tbb_block r) {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                ++progress_counter;
                update_progress(partial_progress);
                const map_type* clone = map.clone();

                size_t i = n % res[0];
                size_t j = n / res[0];
                nvis::vec2 x(bounds.min()[0] + h[0]*(double)i, bounds.min()[1] + h[1]*(double)j);

                try {
                    pos_bwd[n] = metric2d.modulo(clone->map(pos_bwd[n], at-to));
                } catch (...) {
                }

                delete clone;
            }
        });
        partial_progress.stop();
        std::cout << '\n';
        std::cout << partial_progress << '\n';

        if (!add_ftle) {
            std::fill(ftle_fwd.begin(), ftle_fwd.end(), 0);
            std::fill(ftle_bwd.begin(), ftle_bwd.end(), 0);
        }
        if (!add_ridge) {
            std::fill(ridge_strength_fwd.begin(), ridge_strength_fwd.end(), 0);
            std::fill(ridge_strength_bwd.begin(), ridge_strength_bwd.end(), 0);
        }
        progress_counter = 0;
        partial_progress.start(npoints, "FTLE computation");
        tbb::parallel_for(tbb_block(0,npoints), [&](tbb_block r) {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                ++progress_counter;
                update_progress(partial_progress);

                double lmax_f = 0., lmax_b = 0.;
                double ftle_f = 0., ftle_b = 0.;

                double lm = lmax(n, pos_fwd, metric2d);
                lmax_f = log(std::max(1., lm));
                ftle_f = lmax_f;

                lm = lmax(n, pos_bwd, metric2d);
                lmax_b = log(std::max(1., lm));
                ftle_b = lmax_b;
                if (!std::isinf(ftle_f) && !std::isnan(ftle_f)) {
                    ftle_fwd[n] += ftle_f;
                }
                if (!std::isinf(ftle_b) && !std::isnan(ftle_b)) {
                    // std::cerr << ftle_b << " in backward\n";
                    ftle_bwd[n] += ftle_b;
                }
            }
        });
        partial_progress.stop();
        std::cout << '\n';
        std::cout << partial_progress << '\n';

        progress_counter = 0;
        partial_progress.start(npoints, "Ridge strength computation");
        tbb::parallel_for(tbb_block(0,npoints), [&](tbb_block r) {
            for (size_t n=r.begin() ; n!=r.end() ; ++n) {
                ++progress_counter;
                update_progress(partial_progress);

                double lmf = std::min(hlmin(n, ftle_fwd), (double)0);
                double lmb = std::min(hlmin(n, ftle_bwd), (double)0);

                if (!std::isinf(lmf) && !std::isnan(lmf)) {
                    ridge_strength_fwd[n] += fabs(lmf);
                }
                if (!std::isinf(lmb) && !std::isnan(lmb)) {
                    ridge_strength_bwd[n] += fabs(lmb);
                }
            }
        });
        partial_progress.stop();
        std::cout << '\n';
        std::cout << partial_progress << '\n';

        float *ftle = (float*)malloc(2*npoints*sizeof(float));
        float *strength = (float*)malloc(2*npoints*sizeof(float));
        for (size_t i=0; i<npoints; ++i) {
            ftle[2*i] = ftle_fwd[i];
            ftle[2*i+1] = ftle_bwd[i];
            strength[2*i] = ridge_strength_fwd[i];
            strength[2*i+1] = ridge_strength_bwd[i];
        }

        std::ostringstream os;
        int w=std::to_string(maxit).size();
        os << out_basename << "-ftle-";
        if (add_ftle) os << "cumulative-";
        os << "p=" << std::setw(w) << std::setfill('0') << to << "_from_" << maxit
           << "-K=" << K  << "-res=" << res[0] << "x" << res[1]
           << ".nrrd";
        Nrrd* nout = nrrdNew();
        size_t sizes[3] = { 2, res[0], res[1] };
        int kinds[3] = {nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace};
        int centers[3] = {nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell};
        double spc[3] = {std::numeric_limits<double>::quiet_NaN(), h[0], h[1]};
        if (nrrdWrap_nva(nout, ftle, nrrdTypeFloat, 3, sizes)) {
            std::cerr << "ERROR while wrapping data: " << biffGetDone(NRRD)
                      << std::endl;
            if (ftle) {
                delete[] ftle;
            } else {
                nrrdNuke(nout);
            }
            exit(-1);
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoKind, kinds);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, centers);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        spurt::nrrd_utils::writeNrrd(nout, os.str(), true);
        nrrdNuke(nout);
        std::cout << "exported " << os.str() << std::endl;

        os.clear();
        os.str("");
        os << out_basename << "-ftle_ridge_strength-";
        if (add_ridge) os << "cumulative-";
        os << "p=" << std::setw(w) << std::setfill('0') << to << "_from_"
           << maxit << "-K=" << K  << "-res=" << res[0] << "x" << res[1]
           << ".nrrd";
        nout = nrrdNew();
        if (nrrdWrap_nva(nout, strength, nrrdTypeFloat, 3, sizes)) {
            std::cerr << "ERROR while wrapping data: " << biffGetDone(NRRD)
                      << std::endl;
            if (ftle) {
                delete[] ftle;
            } else {
                nrrdNuke(nout);
            }
            exit(-1);
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoKind, kinds);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, centers);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        spurt::nrrd_utils::writeNrrd(nout, os.str(), true);
        nrrdNuke(nout);
        std::cout << "exported " << os.str() << std::endl;
        at=to;

        std::cout << "Total computation time after " << at << " iterations:\n";
        std::cout << total_progress << '\n';
    }

    return 0;
}
