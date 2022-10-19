#include "divergence_cleaning.hpp"
#include <vector>
#include <teem/hest.h>
#include <iostream>
#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>
#include <assert.h>
#include <gmm/gmm.h>
#include <gmm/gmm_least_squares_cg.h>

#include <complex>
#include <data/grid.hpp>
#include <data/raster.hpp>

char*   in, *out;
double  eps;
int     maxiter, solverid;

using namespace div_cleaning;
typedef xavier::raster_grid<3, double> raster_type;
typedef xavier::grid<double, 3> grid_type;
typedef xavier::raster_data<double, 3, double, size_t> sdata_type;
typedef xavier::raster_data<nvis::vec3, 3, double, size_t> vdata_type;
typedef double scalar_type;

void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",  "input",    airTypeString,  1,  1,  &in,        NULL,       "input file name");
    hestOptAdd(&hopt, "o",  "output",   airTypeString,  1,  1,  &out,       NULL,       "output file name");
    hestOptAdd(&hopt, "e",  "eps",      airTypeDouble,  0,  1,  &eps,       "1.0e-16",  "targeted accuracy");
    hestOptAdd(&hopt, "m",  "max iter", airTypeInt,     0,  1,  &maxiter,   "50",       "max number of solver iterations");
    hestOptAdd(&hopt, "s",  "solver id", airTypeInt,    0,  1,  &solverid,  "0",        "type of solver (0: CG, 1: singular CG, 2: LS CG, 3: GMRES, 4: QMR)");
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Apply divergence cleaning to a 3D vector field defined over a uniform grid",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc,  const char* argv[])
{
    using namespace xavier;
    
    initialize(argc, argv);
    
    Nrrd* nin = nrrdNew();
    nin = xavier::readNrrd(in);
    
    // verify data type
    if (nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        return -1;
    }
    
    std::vector<scalar_type> __scalar_in;
    xavier::to_vector(__scalar_in, nin);
    std::vector<vec_type> __vec_in(__scalar_in.size() / 3);
    for (int i = 0 ; i < __vec_in.size() ; ++i) {
        __vec_in[i] = vec_type(__scalar_in[3*i], __scalar_in[3*i+1], __scalar_in[3*i+2]);
    }
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    vdata_type input(domain, __vec_in);
    
    // compute divergence for all interior vertices
    raster_type raster(dims, nvis::vec3(0), spc);
    sdata_type div(raster);
    std::fill(div.begin(), div.end(), 0);
    divergence(div, input);
    
    double maxdiv = *std::max_element(div.begin(), div.end());
    double mindiv = *std::min_element(div.begin(), div.end());
    
    double initmaxdiv = std::max(fabs(maxdiv), fabs(mindiv));
    
    double sum = 0;
    for (int i = 0 ; i < div.size() ; ++i) {
        sum += fabs(div[i]);
    }
    sum /= (double)div.size();
    
    double initmeandiv = sum;
    
    // Poisson system with homogeneous Neumann boundary conditions along x-axis
    sparse_matrix A(domain.size(), domain.size());
    dense_vector RHS(domain.size(), 0);
    const ivec_type dc[] = { ivec_type(1, 0, 0),
                             ivec_type(0, 1, 0),
                             ivec_type(0, 0, 1)
                           };
    double invhsq[] = {1. / (spc[0]*spc[0]), 1 / (spc[1]*spc[1]), 1 / (spc[2]*spc[2])};
    
    for (int k = 0 ; k < dims[0] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 0 ; i < dims[0] ; ++i) {
                ivec_type c(i, j, k);
                int row = domain.index(c);
                if (i == 0) {
                    int col = domain.index(c + dc[0]);
                    A(row, row) = 1 / spc[0];
                    A(row, col) = -1 / spc[0];
                    RHS[row] = 0;
                } else if (i == dims[0] - 1) {
                    int col = domain.index(c - dc[0]);
                    A(row, row) = 1 / spc[0];
                    A(row, col) = -1 / spc[0];
                    RHS[row] = 0;
                } else {
                    A(row, row) = 0;
                    for (int dim = 0 ; dim < 3 ; ++dim) {
                        int col0 = domain.index(c - dc[dim]);
                        int col1 = domain.index(c + dc[dim]);
                        A(row, col0) = invhsq[dim];
                        A(row, col1) = invhsq[dim];
                        A(row, row) -= 2 * invhsq[dim];
                    }
                    RHS[row] = div(c);
                }
            }
        }
    }
    std::cerr << "matrix set\n";
    if (domain.size() < 100) {
        std::cerr << "matrix has size " << A.nrows() << " x " << A.ncols() << '\n';
        for (int i = 0 ; i < A.nrows() ; ++i) {
            std::cerr << i << " (";
            for (int j = 0 ; j < A.ncols() ; ++j) {
                std::cerr << '\t' << A(i, j);
            }
            std::cerr << ")\n";
        }
    }
    
    gmm::iteration iter(eps, 1, maxiter);
    dense_vector X(domain.size(), 0);
    gmm::ildlt_precond<sparse_matrix> PC(A);
    std::cout << "preconditioner set" << std::endl;
    size_t restart = 50;
    double hmin = *std::min(spc.begin(), spc.end());
    
    switch (solverid) {
        case 0: // CG
            gmm::cg(A, X, RHS, PC, iter);
            break;
        case 1: // singular CG
            singular_cg(hmin*hmin, A, X, RHS, PC, iter);
            break;
        case 2: // LS CG
            gmm::least_squares_cg(A, X, RHS, iter);
            break;
        case 3: // GMRES
            gmm::gmres(A, X, RHS, PC, restart, iter);
            break;
        case 4: // QMR
            gmm::qmr(A, X, RHS, PC, iter);
            break;
    }
    
    sdata_type    phi(raster, X.begin(), X.end());
    vdata_type    grad(raster);
    std::fill(grad.begin(), grad.end(), nvis::vec3(0));
    gradient(grad, phi);
    
    std::cerr << "after divergence cleaning:\n";
    for (int i = 0 ; i < raster.size() ; ++i) {
        input[i] -= grad[i];
    }
    
    sdata_type error(raster);
    std::fill(error.begin(), error.end(), 0);
    divergence(error, input);
    double _maxdiv = std::max(fabs(*std::max_element(error.begin(), error.end())),
                              fabs(*std::min_element(error.begin(), error.end())));
    double _avgdiv = 0;
    for (int i = 0 ; i < raster.size() ; ++i) {
        _avgdiv += fabs(error[i]);
    }
    _avgdiv /= raster.size();
    sum /= (double)raster.size();
    
    // double init_max = std::max(fabs(mindiv), fabs(maxdiv));
    std::cerr << "INPUT: max divergence was " << initmaxdiv
              << ", mean divergence was " << initmeandiv << '\n';
    std::cerr << "OUTPUT: max divergence is " << _maxdiv
              << ", mean divergence = " << _avgdiv << '\n';
    std::cerr << "max divergence reduced by " << 100*(1-_maxdiv/initmaxdiv) << "%, "
              << " mean divergence reduced by " << 100*(1-_avgdiv/initmeandiv) << "%.\n";
              
    double* res = (double*)calloc(3 * dims[0] * dims[1] * dims[2], sizeof(double));
    // double *_grad = (double*)calloc(3 * dims[0] * dims[1] * dims[2], sizeof(double));
    // double *_phi = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    // double *_rhs = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    // double *_RHS = (double*)calloc((dims[0] - 2) * (dims[1] - 2) * (dims[2] - 2), sizeof(double));
    // double *_div = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    for (int i = 0 ; i < raster.size() ; ++i) {
        vec_type v = input[i];
        // vec_type g = grad[i];
        for (int d=0 ; d<3 ; ++d) {
            res[3*i+d] = v[d];
            // _grad[3*i+d] = g[d];
        }
    }
    // for (int i = 0 ; i < domain.size() ; ++i) {
    //  _phi[i] = phi.data()[i];
    //  _rhs[i] = rhs[i];
    //  _div[i] = div.data()[i];
    // }
    
    size_t __size[4] = {3, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]};
    double __spc[4] = {airNaN(), spc[0], spc[1], spc[2]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, res, nrrdTypeDouble, 4, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave(out, nout, NULL);
    //
    // nrrdNuke(nout);
    // nout = nrrdNew();
    // nrrdWrap_nva(nout, _grad, nrrdTypeDouble, 4, __size);
    // nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    // nrrdSave("gradphi.nrrd", nout, NULL);
    //
    // nrrdNuke(nout);
    // nout = nrrdNew();
    // nrrdWrap_nva(nout, _phi, nrrdTypeDouble, 3, &__size[1]);
    // nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    // nrrdSave("phi.nrrd", nout, NULL);
    //
    // nrrdNuke(nout);
    // nout = nrrdNew();
    // nrrdWrap_nva(nout, _rhs, nrrdTypeDouble, 3, &__size[1]);
    // nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    // nrrdSave("div.nrrd", nout, NULL);
    //
    // nrrdNuke(nout);
    // nout = nrrdNew();
    // nrrdWrap_nva(nout, _div, nrrdTypeDouble, 3, &__size[1]);
    // nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    // nrrdSave("finaldiv.nrrd", nout, NULL);
    //
    // for (int i = 0 ; i < domain.size() ; ++i) {
    //  _RHS[i] = RHS[i];
    // }
    // nrrdNuke(nout);
    // nout = nrrdNew();
    // __size[1] -= 2;
    // nrrdWrap_nva(nout, _RHS, nrrdTypeDouble, 3, &__size[1]);
    // nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    // nrrdSave("RHS.nrrd", nout, NULL);
    
    return 0;
}











