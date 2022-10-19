#include <math/poisson_gmm.hpp>
#include <data/grid.hpp>
#include <image/nrrd_wrapper.hpp>
#include <teem/hest.h>

#include <gmm/gmm_solver_cg.h>
#include <gmm/gmm_precond_ildlt.h>
#include <gmm/gmm_inoutput.h>

char*    in, *out;
double  eps;
int     maxiter;

#if 0
typedef Eigen::SparseMatrix<double> matrix_t;
typedef Eigen::Triplet<double> triplet_t;
typedef Eigen::VectorXd vector_t;
#endif

void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me = argv[0];
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",  "input",    airTypeString,  1,  1,  &in,        NULL,       "input file name");
    hestOptAdd(&hopt, "o",  "output",   airTypeString,  1,  1,  &out,       NULL,       "output file name");
    hestOptAdd(&hopt, "e",  "eps",      airTypeDouble,  0,  1,  &eps,       "1.0e-16",  "integration length for flow map computation");
    hestOptAdd(&hopt, "m",  "max iter", airTypeInt,     0,  1,  &maxiter,   "50",       "max number of solver iterations");
    std::string comment("Apply divergence cleaning to a 3D vector field defined over a uniform grid");
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, comment.c_str(),
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef xavier::raster_grid<3, double>  grid_type;
typedef grid_type::coord_type           ivec_type;
typedef grid_type::size_type            index_type;
typedef grid_type::point_type           vec_type;

typedef gmm::wsvector<double>           sparse_vector;
typedef gmm::row_matrix<sparse_vector>  sparse_matrix;
typedef std::vector<double>             dense_vector;

inline void divergence(std::vector<double>& out,
                       const std::vector<double>& in,
                       const grid_type& domain)
{
    ivec_type dims = domain.resolution();
    
    // compute divergence for all vertices
    // periodicity of the mesh in Y and Z
    std::fill(out.begin(), out.end(), 0.);
    for (int k = 1 ; k < dims[2]-1 ; ++k) {
        for (int j = 1 ; j < dims[1]-1 ; ++j) {
            for (int i = 1 ; i < dims[0] - 1 ; ++i) {
                ivec_type c(i, j, k);
                double div = 0;
                for (int d = 0 ; d < 3 ; ++d) {
                    nvis::ivec3 __hi(c), __lo(c);
                    ++__hi[d];
                    --__lo[d];
                    if (__lo[d] == -1) {
                        __lo[d] = dims[d] - 1;
                    }
                    if (__hi[d] == dims[d]) {
                        __hi[d] = 0;
                    }
                    
                    index_type hi = 3 * domain.index(__hi) + d;
                    index_type lo = 3 * domain.index(__lo) + d;
                    
                    div += (in[hi] - in[lo]) / (2.*domain.spacing()[d]);
                }
                out[domain.index(c)] = div;
            }
        }
    }
}

inline void poisson(sparse_matrix& A, const grid_type& domain)
{
    ivec_type dims = domain.resolution();
    
    size_t ndof = (dims[0] - 2) * dims[1] * dims[2];
    
    A.resize(ndof, ndof);
    
    //
    // Type h[] = { domain.spacing()[0], domain.spacing()[1], domain.spacing()[2] };
    // typename grid_type::vec_type __dims(dims);
    // __dims[0] -= 2;
    // grid_type interior(__dims, domain.spacing());
    //
    // for (int k = 0 ; k < dims[2] ; ++k) {
    //  for (int j = 0 ; j < dims[1] ; ++j) {
    //      for (int i = 1 ; i < dims[0] - 1; ++i) {
    //          ivec_type coord(i, j, k);
    //          index_type row = interior.index(i - 1, j, k);
    //
    //          // diagonal terms
    //          A(row, row) = 2 * (1 / h[0] + 1 / h[1] + 1 / h[2]); // diagonal term
    //          rhs[row] = -grhs[row];
    //          for (int dim = 0 ; dim < 3 ; ++dim) {
    //              for (int s = 0 ; s < 2 ; ++s) {
    //                  nvis::ivec3 c(coord);
    //                  c[dim] += (s ? 1 : -1);
    //                  if (c[dim] == -1) c[dim] = dims[dim] - 1;
    //                  else if (c[dim] == dims[dim]) c[dim] = 0;
    //                  else if (c[0] == 0 || c[0] == dims[0] - 1) continue;
    //                  int col = interior.index(c[0] - 1, c[1], c[2]);
    //                  A(row, col) = -1 / h[dim];
    //              }
    //          }
    //      }
    //  }
    // }
}

inline nvis::vec3 gradient(const std::vector<double>& vals,
                           const grid_type& domain,
                           const ivec_type& c)
{
    nvis::vec3 g(0, 0, 0);
    ivec_type dims = domain.resolution();
    // int id = domain.index(c);
    for (int i = 0 ; i < 3 ; ++i) {
        nvis::ivec3 hi(c), lo(c);
        ++hi[i];
        --lo[i];
        double h = 2.*domain.spacing()[i];
        if (i == 1) {
            if (hi[1] == dims[1]) {
                hi[1] = dims[1] - 1;
                h = domain.spacing()[1];
            } else if (lo[1] == -1) {
                lo[1] = 0;
                h = domain.spacing()[1];
            }
        } else if (hi[i] == dims[i]) {
            hi[i] = 0;
        } else if (lo[i] == -1) {
            lo[i] = dims[i] - 1;
        }
        
        index_type __lo = domain.index(lo);
        index_type __hi = domain.index(hi);
        g[i] = (vals[__hi] - vals[__lo]) / h;
    }
    return g;
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
                int id = domain.index(i, j, k);
                nvis::vec3 g = gradient(in, domain, ivec_type(i, j, k));
                out[3*id  ] = g[0];
                out[3*id+1] = g[1];
                out[3*id+2] = g[2];
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
    
    gmm::clear(x);      // x := 0
    gmm::copy(b, r);    // r := b 
    gmm::mult(P, r, z); // z := Pr
    
    rho = gmm::vect_sp(z, r); // rho := z.r
    gmm::copy(z, p);          // p := z
    
    std::cerr << "in singular_cg\n";
    
    while (!iter.finished_vect(r)) {
        if (!iter.first()) {
            gmm::mult(P, r, z); // z := Pr
            rho = gmm::vect_sp(z, r); // rho := z.r
            gmm::add(z, gmm::scaled(p, rho / rho_1), p); // p = z + rho/rho1*p
        }
        
        gmm::mult(A, p, q); // q := Ap
        // gmm::add(q, gmm::scaled(c, mu*gmm::vect_sp(c, p)), q);
        
        a = rho / gmm::vect_sp(q, p); // a := rho/||p||_A
        gmm::add(gmm::scaled(p,  a), x); // x += a*p
        gmm::add(gmm::scaled(q, -a), r); // r -= a*q
        
        rho_1 = rho;
        ++iter;
    }
}

#if 0
void singular_cg_Eigen(const matrix_t>& A, vector_t& x, const vector_t& b) {                         
    Eigen::ConjugateGradient<matrix_t> cg;
    cg.setMaxIterations(1);
    cg.compute(A);
    
    int n=0;
    for (; n<maxiter ; ++n) {
        x = cg.solve(b);
        std::cout << n << "\testimated error: " << cg.error() << '\n';
        if (cg.error() < eps) return;
    }
    
    std::cout << "Cg failed\n";
}
#endif

int main(int argc, const char* argv[])
{
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
    grid_type domain(dims, vec_type(0), spc);
    
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
              
    std::vector<double> phi(domain.size());
    
    ivec_type __dims(dims);
    __dims -= vec_type(2, 2, 2);
    grid_type interior(__dims, vec_type(0), domain.spacing());
    size_t size = interior.size(); // size of the interior of the domain
    
    // gmm
    gmm::iteration iter(eps, 1, maxiter);
    sparse_matrix A(size, size);
    dense_vector X(size, 0), RHS(size, 0);

#if 0    
    // Eigen
    matrix_t A__(size, size);
    vector_t x__(size), rhs__(size);
#endif
    
    xavier::setup_system_tokamak_dirichlet_interior<grid_type, double, sparse_matrix>(domain, A);
    for (int i = 0 ; i < size ; ++i) {
        ivec_type ic = interior.coordinates(i);
        int gi = domain.index(ic + ivec_type(1, 1, 1));
        RHS[i] = rhs[gi];
    }
    
    std::cerr << "matrix set\n";
    
    gmm::ildlt_precond<sparse_matrix> PC(A);
    std::cout << "preconditioner set" << std::endl;
    
    double hmin = *std::min(spc.begin(), spc.end());
    
    singular_cg(hmin*hmin, A, X, RHS, PC, iter);
    
    for (int i = 0 ; i < size ; ++i) {
        ivec_type gc = domain.coordinates(i);
        if (gc[0] == 0 || gc[0] == dims[0] - 1 ||
                gc[1] == 0 || gc[1] == dims[1] - 1 ||
                gc[2] == 0 || gc[2] == dims[2] - 1 ) {
            phi[i] = 0;
            continue;
        }
        int li = interior.index(gc - ivec_type(1, 1, 1));
        phi[li] = X[i];
    }
    
    std::cerr << "after divergence cleaning:\n";
    std::vector<double> grad(3*domain.size());
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
    double* _RHS = (double*)calloc((dims[0]-2) * (dims[1]-2) * (dims[2]-2), sizeof(double));
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
    
    for (int i=0 ; i<size ; ++i) {
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




























