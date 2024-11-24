#include <math/poisson_gmm.hpp>
#include <data/image.hpp>
#include <image/nrrd_wrapper.hpp>

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

typedef spurt::raster_grid<long, double, 3, lvec3, vec3> grid_type;
typedef grid_type::coord_type                            ivec_type;
typedef grid_type::size_type                             index_type;
typedef grid_type::point_type                            vec_type;

typedef gmm::wsvector<double>           sparse_vector;
typedef gmm::row_matrix<sparse_vector>  sparse_matrix;
typedef std::vector<double>             dense_vector;


inline double divergence(const ivec_type& c, const std::vector<double>& in,
                         const grid_type& domain)
{
    double div = 0;
    for (int i = 0 ; i < 3 ; ++i) {
        ivec_type step(0, 0, 0);
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
        ivec_type step(0, 0, 0);
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
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
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
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
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
    
    // verify data type
    if (nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        return -1;
    }
    
    std::vector<double> flow;
    spurt::to_vector(flow, nin);
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    
    int size = domain.size();
    
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
              
    // set up matrix
    sparse_matrix A(size, size);
    dense_vector RHS(size, 0);
    double invhsq[] = {1. / (spc[0]*spc[0]), 1 / (spc[1]*spc[1]), 1 / (spc[2]*spc[2])};
    for (int row = 0 ; row < size ; ++row) {
        nvis::ivec3 c = domain.coordinates(row);
        A(row, row) = 0;
        if (c[0] > 0 && c[0] < dims[0] - 1) {
            // discrete Laplace operator
            for (int dim = 0 ; dim < 3 ; ++dim) {
                ivec_type lowc(c), hic(c);
                lowc[dim] -= 1;
                hic[dim] += 1;
                if (dim > 0) {
                    if (lowc[dim] == -1) {
                        lowc[dim] = dims[dim] - 1;
                    }
                    if (hic[dim] == dims[dim]) {
                        hic[dim] = 0;
                    }
                    A(row, domain.index(lowc)) = invhsq[dim];
                    A(row, domain.index(hic)) = invhsq[dim];
                } else {
                    A(row, domain.index(lowc)) = invhsq[dim];
                    A(row, domain.index(hic)) = invhsq[dim];
                }
                A(row, row) -= 2 * invhsq[dim];
            }
            RHS[row] = rhs[row];
        }
        // homogeneous Neumann boundary condition
        else if (c[0] == 0) {
            A(row, row) = 1. / spc[0];
            A(row, row + 1) = -1. / spc[0];
            RHS[row] = 0;
        } else {
            A(row, row) = 1. / spc[0];
            A(row, row - 1) = -1. / spc[0];
            RHS[row] = 0;
        }
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
    
    std::vector<double> phi(domain.size(), 0);
    for (int i = 0 ; i < size ; ++i) {
        phi[i] = X[i];
    }
    
    
    std::cerr << "after divergence cleaning:\n";
    std::vector<double> grad(3*domain.size(), 0);
    gradient(grad, phi, domain);
    for (int i = 0 ; i < flow.size() ; ++i) {
        flow[i] -= grad[i];
    }
    std::vector<double> div(domain.size(), 0);
    divergence(div, flow, domain);
    maxdiv = *std::max_element(div.begin(), div.end());
    mindiv = *std::min_element(div.begin(), div.end());
    sum = 0;
    for (int i = 0 ; i < div.size() ; ++i) {
        sum += fabs(div[i]);
    }
    sum /= (double)rhs.size();
    std::cerr << "OUTPUT: min divergence is " << mindiv << ", max divergence is " << maxdiv
              << ", mean divergence = " << sum << '\n';
              
    double* res = (double*)calloc(3 * dims[0] * dims[1] * dims[2], sizeof(double));
    double* _grad = (double*)calloc(3 * dims[0] * dims[1] * dims[2], sizeof(double));
    double* _phi = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    double* _rhs = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
    double* _RHS = (double*)calloc(dims[0] * dims[1] * dims[2], sizeof(double));
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
    nrrdWrap_nva(nout, _RHS, nrrdTypeDouble, 3, &__size[1]);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, &__spc[1]);
    nrrdSave("RHS.nrrd", nout, NULL);
    
    
    
    return 0;
}

























































