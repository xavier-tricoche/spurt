#include <math/poisson_gmm.hpp>
#include <data/grid.hpp>
#include <image/nrrd_wrapper.hpp>
#include <teem/nrrd.h>

#include <gmm/gmm_solver_cg.h>
#include <gmm/gmm_precond_ildlt.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_least_squares_cg.h>

char*   in, *out;
double  eps;
int     maxiter, solverid;

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
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Apply divergence cleaning to a 3D vector field defined over a uniform grid",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef spurt::raster_grid<3, double>  grid_type;
typedef grid_type::coord_type           ivec_type;
typedef grid_type::size_type            index_type;
typedef grid_type::vec_type             vec_type;

typedef gmm::wsvector<double>           sparse_vector;
typedef gmm::row_matrix<sparse_vector>  sparse_matrix;
typedef std::vector<double>             dense_vector;

struct index_converter {
    index_converter(const ivec_type& size) : sg(size) {
        si = sg - ivec_type(2, 2, 2);
    }
    
    ivec_type glob2in(const ivec_type& g) const {
        return g - ivec_type(1, 1, 1);
    }
    
    ivec_type in2glob(const ivec_type& i) const {
        return i + ivec_type(1, 1, 1);
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


inline double divergence(const ivec_type& c, const std::vector<double>& in,
                         const grid_type& domain)
{
    double div = 0;
    for (int i = 0 ; i < 3 ; ++i) {
        ivec_type step(0);
        step[i] = 1;
        int hi = 3 * domain.index(c + step) + i;
        int lo = 3 * domain.index(c - step) + i;
        div += 0.5 / domain.spacing()[i] * (in[hi] - in[lo]);
    }
    return div;
}

inline vec_type gradient(const ivec_type& c, const std::vector<double>& in,
                         const grid_type& domain)
{
    vec_type g;
    for (int i = 0 ; i < 3 ; ++i) {
        ivec_type step(0);
        step[i] = 1;
        int hi = domain.index(c + step);
        int lo = domain.index(c - step);
        g[i] = 0.5 / domain.spacing()[i] * (in[hi] - in[lo]);
    }
    return g;
}

inline void divergence(std::vector<double>& out,
                       const std::vector<double>& in,
                       const grid_type& domain)
{
    ivec_type dims = domain.resolution();
    
    // compute divergence for all interior vertices
    std::fill(out.begin(), out.end(), 0.);
    for (int k = 1 ; k < dims[2] - 1 ; ++k) {
        for (int j = 1 ; j < dims[1] - 1 ; ++j) {
            for (int i = 1 ; i < dims[0] - 1 ; ++i) {
                ivec_type c(i, j, k);
                out[domain.index(c)] = divergence(c, in, domain);
            }
        }
    }
}

inline void gradient(std::vector<double>& out,
                     const std::vector<double>& in,
                     const grid_type& domain)
{
    ivec_type dims = domain.resolution();
    // compute gradient of input scalar field at all interior vertices
    std::fill(out.begin(), out.end(), 0.);
    for (int k = 1 ; k < dims[2] - 1 ; ++k) {
        for (int j = 1 ; j < dims[1] - 1 ; ++j) {
            for (int i = 1 ; i < dims[0] - 1 ; ++i) {
                ivec_type c(i, j, k);
                int n = domain.index(i, j, k);
                nvis::vec3 g = gradient(c, in, domain);
                out[3*n  ] = g[0];
                out[3*n+1] = g[1];
                out[3*n+2] = g[2];
            }
        }
    }
}

template<typename PC>
void singular_cg(const double& mu,
                 const sparse_matrix& A,
                 dense_vector& x,
                 const dense_vector& b,
                 const PC& P,
                 gmm::iteration& iter)
{
    const unsigned int N = gmm::vect_size(x);
    
    double rho, rho_1(0), a;
    dense_vector p(N), q(N), r(N), z(N), c(N, 1.0);
    
    iter.set_rhsnorm(gmm::vect_norm2(b));
    
    gmm::clear(x);
    gmm::copy(b, r);
    gmm::mult(P, r, z);
    
    rho = gmm::vect_sp(z, r);
    gmm::copy(z, p);
    
    std::cerr << "in singular_cg\n";
    
    while (!iter.finished_vect(r)) {
        if (!iter.first()) {
            gmm::mult(P, r, z);
            rho = gmm::vect_sp(z, r);
            gmm::add(z, gmm::scaled(p, rho / rho_1), p);
        }
        
        gmm::mult(A, p, q);
        gmm::add(q, gmm::scaled(c, mu*gmm::vect_sp(c, p)), q);
        
        a = rho / gmm::vect_sp(q, p);
        gmm::add(gmm::scaled(p,  a), x);
        gmm::add(gmm::scaled(q, -a), r);
        
        rho_1 = rho;
        ++iter;
    }
}

int main(int argc, const char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin = nrrdNew();
    nin = spurt::readNrrd(in);
    spurt::nrrd_traits traits(nin);
    
    // verify data type
    if (traits.dim() != 4 || traits.sizes()[0] != 3) {
        std::cerr << "invalid input NRRD file.\n";
        return -1;
    }
    
    std::vector<double> flow;
    spurt::to_vector(flow, nin);
    ivec_type dims(traits.sizes()[1], traits.sizes()[2], traits.sizes()[3]);
    nvis::vec3 spc(traits.spacings()[1], traits.spacings()[2], traits.spacings()[3]);
    nvis::vec3 orig(traits.mins()[1], traits.mins()[2], traits.mins()[3]);
    grid_type domain(dims, orig, spc);
    
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
    
    index_converter converter(domain.resolution());
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
    
    nvis::bounding_box<nvis::ivec3> valid_coords(nvis::ivec3(0, 0, 0), nvis::ivec3(dims[0] - 3, dims[1] - 3, dims[2] - 3));
    
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
            if (valid_coords.inside(lowc)) {
                A(row, converter.in(lowc)) = invhsq[dim];
            }
            if (valid_coords.inside(hic)) {
                A(row, converter.in(hic)) = invhsq[dim];
            }
            A(row, row) -= 2 * invhsq[dim];
        }
        RHS[row] = rhs[converter.glob(gc)];
    }
    std::cerr << "matrix set\n";
    if (size < 50) {
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
    // gmm::ilut_precond<sparse_matrix> PC(A);
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
    __size[2] -= 2;
    __size[3] -= 2;
    nrrdWrap_nva(nout, _RHS, nrrdTypeDouble, 3, &__size[1]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    nrrdSave("RHS.nrrd", nout, NULL);
    
    
    
    return 0;
}



















































