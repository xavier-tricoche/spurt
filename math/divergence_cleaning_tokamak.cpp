#include "divergence_cleaning.hpp"
#include <teem/hest.h>
#include <image/nrrd_wrapper.hpp>
#include "poisson_gmm.hpp"

char*    in, *out;
double  eps;
int     maxiter, solverid;

using namespace div_cleaning;

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
    hestOptAdd(&hopt, "e",  "eps",      airTypeDouble,  0,  1,  &eps,       "1.0e-16",  "integration length for flow map computation");
    hestOptAdd(&hopt, "m",  "max iter", airTypeInt,     0,  1,  &maxiter,   "50",       "max number of solver iterations");
    hestOptAdd(&hopt, "s",  "solver id", airTypeInt,    0,  1,  &solverid,  "0",        "type of solver (0: CG, 1: singular CG, 2: LS CG, 3: GMRES, 4: QMR)");
    
    std::string comment = "Apply divergence cleaning to a 3D vector field defined over a uniform grid";
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, comment.c_str(),
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct index_converter {
    index_converter(const ivec_type& size) : sg(size) {
        si = sg - ivec_type(2, 0, 0);
    }
    
    ivec_type glob2in(const ivec_type& g) const {
        return g - ivec_type(1, 0, 0);
    }
    
    ivec_type in2glob(const ivec_type& i) const {
        return i + ivec_type(1, 0, 0);
    }
    
    int glob2in(int g) const {
        return in(glob2in(glob(g)));
    }
    
    int in2glob(int i) const {
        return glob(in2glob(in(i)));
    }
    
    int glob(const ivec_type& id) const {
        return id[0] + sg[0]*(id[1] + sg[1]*id[2]);
    }
    
    ivec_type glob(int n) const {
        int x = n % sg[0];
        n /= sg[0];
        int y = n % sg[1];
        int z = n / sg[1];
        return ivec_type(x, y, z);
    }
    
    int in(const ivec_type& id) const {
        return id[0] + si[0]*(id[1] + si[1]*id[2]);
    }
    
    ivec_type in(int n) const {
        int x = n % si[0];
        n /= si[0];
        int y = n % si[1];
        int z = n / si[1];
        return ivec_type(x, y, z);
    }
    
    ivec_type sg, si;
};

inline bool valid(size_type n) { 
    return n!=static_cast<size_type>(n);
}

int main(int argc, const char* argv[])
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
    
    std::vector<double> flow;
    xavier::to_vector(flow, nin);
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    
    // compute divergence for all interior vertices
    std::vector<double> rhs(domain.size(), 0);
    divergence(rhs, flow, domain);
    
    double maxdiv = *std::max_element(rhs.begin(), rhs.end());
    double mindiv = *std::min_element(rhs.begin(), rhs.end());
    double sum = 0;
    for (int i = 0 ; i < rhs.size() ; ++i) {
        sum += fabs(rhs[i]);
    }
    sum /= (double)rhs.size();
    std::cerr << "INPUT: min divergence is " << mindiv << ", max divergence is " << maxdiv
              << ", mean divergence is " << sum << '\n';
              
    std::vector<double> phi(domain.size(), 0);
    
    index_converter converter(domain.dimensions());
    size_t size = converter.si[0] * converter.si[1] * converter.si[2];
    double __maxdiv = 0;
    double __avgdiv = 0;
    for (int i = 0 ; i < size ; ++i) {
        double v = fabs(rhs[converter.in2glob(i)]);
        __avgdiv += v;
        __maxdiv = std::max(v, maxdiv);
    }
    __avgdiv /= size;
    std::cerr << "INPUT: interior max divergence = " << __maxdiv << ", mean divergence = " << __avgdiv << '\n';
    
    nvis::bounding_box<nvis::ivec3> valid_coords(nvis::ivec3(0, 0, 0), nvis::ivec3(dims[0] - 3, dims[1] - 1, dims[2] - 1));
    
    // set up matrix
    sparse_matrix A(size, size);
    dense_vector RHS(size, 0);
    double invhsq[] = {1. / (spc[0]*spc[0]), 1 / (spc[1]*spc[1]), 1 / (spc[2]*spc[2])};
    for (int row = 0 ; row < size ; ++row) {
        nvis::ivec3 ic = converter.in(row);
        nvis::ivec3 gc = converter.in2glob(ic);
        A(row, row) = 0;
        for (int dim = 0 ; dim < 3 ; ++dim) {
            ivec_type lowc(ic), hic(ic);
            lowc[dim] -= 1;
            hic[dim] += 1;
            if (dim > 0) {
                if (!valid(lowc[dim])) {
                    lowc[dim] = dims[dim] - 1;
                }
                if (hic[dim] == dims[dim]) {
                    hic[dim] = 0;
                }
                A(row, converter.in(lowc)) = invhsq[dim];
                A(row, converter.in(hic)) = invhsq[dim];
            } else {
                if (valid(lowc[dim]) && lowc[dim] < dims[dim] - 2) {
                    A(row, converter.in(lowc)) = invhsq[dim];
                }
                if (valid(hic[dim]) && hic[dim] < dims[dim] - 2) {
                    A(row, converter.in(hic)) = invhsq[dim];
                }
            }
            A(row, row) -= 2 * invhsq[dim];
        }
        RHS[row] = rhs[converter.glob(gc)];
    }
    std::cerr << "matrix set\n";
    if (size < 100) {
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
    dense_vector X(size, 0);
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
    for (int i = 0 ; i < size ; ++i) {
        phi[converter.in2glob(i)] = X[i];
    }
    
    
    std::cerr << "after divergence cleaning:\n";
    std::vector<double> grad(3*domain.size(), 0);
    gradient(grad, phi, domain);
    for (int i = 0 ; i < flow.size() ; ++i) {
        flow[i] -= grad[i];
    }
    std::vector<double> div(domain.size(), 0);
    divergence(div, flow, domain);
    double _maxdiv = 0;
    double _avgdiv = 0;
    for (int i = 0 ; i < size ; ++i) {
        double v = fabs(div[converter.in2glob(i)]);
        _avgdiv += v;
        _maxdiv = std::max(v, maxdiv);
    }
    _avgdiv /= size;
    
    maxdiv = *std::max_element(div.begin(), div.end());
    mindiv = *std::min_element(div.begin(), div.end());
    sum = 0;
    for (int i = 0 ; i < div.size() ; ++i) {
        sum += fabs(div[i]);
    }
    sum /= (double)rhs.size();
    std::cerr << "OUTPUT: min divergence is " << mindiv << ", max divergence is " << maxdiv
              << ", mean divergence = " << sum << '\n';
    std::cerr << "OUTPUT: max interior divergence is " << _maxdiv
              << ", mean divergence = " << _avgdiv << '\n';
              
    double* res = (double*)calloc(3 * dims[0] * dims[1] * dims[2], sizeof(double));
    double* _grad = (double*)calloc(3 * dims[0] * dims[1] * dims[2], sizeof(double));
    double* _phi = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    double* _rhs = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    double* _RHS = (double*)calloc((dims[0] - 2) * (dims[1] - 2) * (dims[2] - 2), sizeof(double));
    double* _div = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    for (int i = 0 ; i < flow.size() ; ++i) {
        res[i] = flow[i];
        _grad[i] = grad[i];
    }
    for (int i = 0 ; i < phi.size() ; ++i) {
        _phi[i] = phi[i];
        _rhs[i] = rhs[i];
        _div[i] = div[i];
    }
    
    size_t __size[4] = {3, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]};
    double __spc[4] = {airNaN(), spc[0], spc[1], spc[2]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, res, nrrdTypeDouble, 4, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave(out, nout, NULL);
    
    nrrdNuke(nout);
    nout = nrrdNew();
    nrrdWrap_nva(nout, _grad, nrrdTypeDouble, 4, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave("gradphi.nrrd", nout, NULL);
    
    nrrdNuke(nout);
    nout = nrrdNew();
    nrrdWrap_nva(nout, _phi, nrrdTypeDouble, 3, &__size[1]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    nrrdSave("phi.nrrd", nout, NULL);
    
    nrrdNuke(nout);
    nout = nrrdNew();
    nrrdWrap_nva(nout, _rhs, nrrdTypeDouble, 3, &__size[1]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    nrrdSave("div.nrrd", nout, NULL);
    
    nrrdNuke(nout);
    nout = nrrdNew();
    nrrdWrap_nva(nout, _div, nrrdTypeDouble, 3, &__size[1]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    nrrdSave("finaldiv.nrrd", nout, NULL);
    
    for (int i = 0 ; i < size ; ++i) {
        _RHS[i] = RHS[i];
    }
    nrrdNuke(nout);
    nout = nrrdNew();
    __size[1] -= 2;
    nrrdWrap_nva(nout, _RHS, nrrdTypeDouble, 3, &__size[1]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    nrrdSave("RHS.nrrd", nout, NULL);
    
    return 0;
}
