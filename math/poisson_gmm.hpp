#include <assert.h>

#include <gmm/gmm.h>
#include <gmm/gmm_least_squares_cg.h>

#include <complex>

#include <data/raster.hpp>

// debug
#include <teem/nrrd.h>

namespace {
template<typename Grid>
struct is_boundary_cube {
    typedef Grid                            grid_type;
    typedef typename grid_type::ivec_type   ivec_type;
    typedef typename grid_type::index_type  index_type;
    
    is_boundary_cube(const Grid& domain) : _domain(domain) {}
    
    template<typename Index>
    bool operator()(const Index& id) const {
        const ivec_type& dims = _domain.resolution()();
        for (int i = 0 ; i < 3 ; ++i) {
            if (!id[i] || id[i] == dims[i] - 1) {
                return true;
            }
        }
        return false;
    }
    
    const Grid& _domain;
};

template<typename Grid, int N>
struct is_boundary_torus {
    typedef Grid                            grid_type;
    typedef typename grid_type::coord_type  ivec_type;
    typedef typename grid_type::size_type   index_type;
    
    is_boundary_torus(const Grid& domain) : _domain(domain) {}
    
    template<typename Index>
    bool operator()(const Index& id) const {
        const ivec_type& dims = _domain.resolution();
        return (!id[N] || id[N] == dims[N] - 1);
    }
    
    const Grid& _domain;
};

struct c2i {
    c2i(const nvis::ivec3& size) : _size(size) {}
    
    int operator()(const nvis::ivec3& id) const {
        return id[0] + _size[0]*(id[1] + _size[1]*id[2]);
    }
    
    int operator()(int i, int j, int k) const {
        return i + _size[0]*(j + _size[1]*k);
    }
    
    nvis::ivec3 _size;
};

}

namespace xavier {

template<typename Grid_, typename Type_, typename Matrix_>
void setup_system_tokamak_dirichlet_tokamak(const Grid_& domain, Matrix_& A,
        std::vector<Type_>& rhs, const std::vector<Type_>& grhs)
{
    // boundary conditions:
    // face 0 (x=xmin): homogeneous Dirichlet
    // face 1 (x=xmax): homogeneous Dirichlet
    // face 2 (y=ymin): periodic
    // face 3 (y=ymax): periodic
    // face 4 (z=zmin): periodic
    // face 5 (z=zmax): periodic
    
    typedef Grid_                          grid_type;
    typedef typename grid_type::coord_type ivec_type;
    typedef typename grid_type::size_type  index_type;
    typedef typename grid_type::point_type vec_type;
    typedef Type_                          value_type;
    
    ivec_type dims = domain.resolution();
    value_type h[] = { domain.spacing()[0], domain.spacing()[1], domain.spacing()[2] };
    vec_type __dims(dims);
    __dims[0] -= 2;
    grid_type interior(__dims, domain.spacing());
    
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 1 ; i < dims[0] - 1; ++i) {
                ivec_type coord(i, j, k);
                index_type row = interior.index(i - 1, j, k);
                
                // diagonal terms
                A(row, row) = 2 * (1 / h[0] + 1 / h[1] + 1 / h[2]); // diagonal term
                rhs[row] = -grhs[row];
                for (int dim = 0 ; dim < 3 ; ++dim) {
                    for (int s = 0 ; s < 2 ; ++s) {
                        nvis::ivec3 c(coord);
                        c[dim] += (s ? 1 : -1);
                        if (c[dim] == -1) {
                            c[dim] = dims[dim] - 1;
                        } else if (c[dim] == dims[dim]) {
                            c[dim] = 0;
                        } else if (c[0] == 0 || c[0] == dims[0] - 1) {
                            continue;
                        }
                        int col = interior.index(c[0] - 1, c[1], c[2]);
                        A(row, col) = -1 / h[dim];
                    }
                }
            }
        }
    }
}

template<typename Grid_, typename Type_, typename Matrix_>
void setup_system_tokamak_dirichlet_interior(const Grid_& domain, Matrix_& A)
{
    typedef Grid_                           grid_type;
    typedef typename grid_type::coord_type  ivec_type;
    typedef typename grid_type::size_type   index_type;
    typedef Type_                           value_type;
    
    ivec_type dims = domain.resolution();
    value_type h[] = { domain.spacing()[0], domain.spacing()[1], domain.spacing()[2] };
    
    c2i converter(dims - ivec_type(2, 2, 2));
    
    for (int k = 1 ; k < dims[2] - 1 ; ++k) {
        for (int j = 1 ; j < dims[1] - 1 ; ++j) {
            for (int i = 1 ; i < dims[0] - 1; ++i) {
                ivec_type gc(i, j, k);          // global coordinates
                ivec_type ic(i - 1, j - 1, k - 1);  // interior coordinates
                int row = converter(ic);
                A(row, row) = -2 * (1 / (h[0] * h[0]) + 1 / (h[1] * h[1]) + 1 / (h[2] * h[2]));
                
                for (int d = 0 ; d < 3 ; ++d) {
                    ivec_type step(0, 0, 0);
                    step[d] = 1;
                    if (ic[d] > 0) {
                        int col = converter(ic - step);
                        A(row, col) = 1 / (h[d] * h[d]);
                    }
                    if (ic[d] < dims[d] - 3)  {
                        int col = converter(ic + step);
                        A(row, col) = 1 / (h[d] * h[d]);
                    }
                }
            }
        }
    }
}


template<typename Grid_, typename Type_, typename Matrix_>
void setup_system_tokamak_neumann_tokamak(const Grid_& domain, Matrix_& A,
        std::vector<Type_>& rhs, const std::vector<Type_>& grhs)
{
    // boundary conditions:
    // face 0 (x=xmin): homogeneous neumann
    // face 1 (x=xmax): homogeneous neumann
    // face 2 (y=ymin): periodic
    // face 3 (y=ymax): periodic
    // face 4 (z=zmin): periodic
    // face 5 (z=zmax): periodic
    
    typedef Grid_                           grid_type;
    typedef typename grid_type::coord_type  ivec_type;
    typedef typename grid_type::size_type   index_type;
    typedef Type_                           value_type;
    
    ivec_type dims = domain.resolution();
    value_type h[] = { domain.spacing()[0], domain.spacing()[1], domain.spacing()[2] };
    
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 0 ; i < dims[0] ; ++i) {
                ivec_type coord(i, j, k);
                index_type row = domain.index(i, j, k);
                
                if (i == 0 || i == dims[0] - 1) {
                    // homogeneous Neumann boundary condition
                    // dphi/dx = 0
                    int hi, lo;
                    ivec_type _hi(coord), _lo(coord);
                    if (i == 0) {
                        ++_hi[0];
                    } else {
                        --_lo[0];
                    }
                    hi = domain.index(_hi[0], _hi[1], _hi[2]);
                    lo = domain.index(_lo[0], _lo[1], _lo[2]);
                    A(row, hi) = 1 / h[0];
                    A(row, lo) = -1 / h[0];
                    rhs[row] = 0;
                } else {
                    // diagonal terms
                    A(row, row) = 2 * (1 / h[0] + 1 / h[1] + 1 / h[2]); // diagonal term
                    rhs[row] = -grhs[row];
                    for (int dim = 0 ; dim < 3 ; ++dim) {
                        for (int s = 0 ; s < 2 ; ++s) {
                            nvis::ivec3 c(coord);
                            c[dim] += (s ? 1 : -1);
                            if (c[dim] == -1) {
                                c[dim] = dims[dim] - 1;
                            } else if (c[dim] == dims[dim]) {
                                c[dim] = 0;
                            }
                            int idx = domain.index(c);
                            A(row, idx) = -1 / h[dim];
                        }
                    }
                }
            }
        }
    }
}

template<typename Grid_, typename Type_>
void solve_poisson_tokamak(const Grid_& domain,
                           const std::vector<Type_>& rhs,
                           std::vector<Type_>& x, 
                           int maxiter = 50, Type_ eps = 1.0e-16)
{
    typedef Grid_                           grid_type;
    typedef typename grid_type::coord_type  ivec_type;
    typedef typename grid_type::size_type   index_type;
    typedef typename grid_type::point_type  vec_type;
    typedef Type_                           value_type;
    
    typedef gmm::wsvector<value_type>       sparse_vector;
    typedef gmm::row_matrix<sparse_vector>  sparse_matrix;
    typedef std::vector<value_type>         dense_vector;
    
    ivec_type dims = domain.resolution();
    ivec_type __dims(dims);
    __dims[0] -= 2;
    grid_type interior(__dims, domain.spacing());
    size_t size = interior.size(); // size of the interior of the domain
    
    gmm::iteration iter(eps, 1, maxiter);
    sparse_matrix tmp(size, size);
    dense_vector X(size, 0), RHS(size, 0);
    setup_system_tokamak_dirichlet_tokamak<Type_>(domain, tmp, RHS, rhs);
    gmm::csc_matrix<Type_> A;
    gmm::copy(tmp, A);
    
    std::cerr << "matrix set\n";
    
    // debug
    Nrrd* nout = nrrdNew();
    double* array = (double*)calloc(size, sizeof(double));
    for (int i = 0 ; i < size ; ++i) {
        array[i] = RHS[i];
    }
    size_t __size[] = {dims[0], dims[1], dims[2]};
    double __spc[] = {domain.spacing()[0], domain.spacing()[1], domain.spacing()[2]};
    nrrdWrap_nva(nout, array, nrrdTypeDouble, 3, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave("RHS.nrrd", nout, NULL);
    nrrdNuke(nout);
    
    gmm::ildlt_precond<gmm::csc_matrix<Type_> > PC(A);
    // // incomplete Cholesky preconditioner
    // gmm::mr_approx_inverse_precond<gmm::csc_matrix<Type_> > PC(A, 10, 1.0e-6);
    
    std::cerr << "PC set\n";
    
    gmm::cg(A, X, RHS, PC, iter);
    
    std::cerr << "system solved\n";
    
    // copy solution to output vector
    std::fill(x.begin(), x.end(), 0);
    for (int i = 0 ; i < size ; ++i) {
        ivec_type c = interior.coordinates(i);
        ++c[0];
        index_type j = domain.index(c);
        x[j] = X[i];
    }
}

template<typename Grid_, typename Type_>
void solve_poisson_interior(const Grid_& domain,
                            const std::vector<Type_>& rhs,
                            std::vector<Type_>& x, 
                            int maxiter = 50, Type_ eps = 1.0e-16)
{
    typedef Grid_                           grid_type;
    typedef typename grid_type::coord_type  ivec_type;
    typedef typename grid_type::size_type   index_type;
    typedef typename grid_type::point_type  vec_type;
    typedef Type_                           value_type;
    
    typedef gmm::wsvector<value_type>       sparse_vector;
    typedef gmm::row_matrix<sparse_vector>  sparse_matrix;
    typedef std::vector<value_type>         dense_vector;
    
    ivec_type dims = domain.resolution();
    ivec_type __dims(dims);
    __dims -= ivec_type(2, 2, 2);
    grid_type interior(__dims, vec_type(0), domain.spacing());
    size_t size = interior.size(); // size of the interior of the domain
    
    gmm::iteration iter(eps, 1, maxiter);
    sparse_matrix tmp(size, size);
    dense_vector X(size, 0), RHS(size, 0);
    setup_system_tokamak_dirichlet_interior<Type_>(domain, tmp);
    gmm::csc_matrix<Type_> A;
    gmm::copy(tmp, A);
    
    std::cerr << "matrix set\n";
    
    // debug
    Nrrd* nout = nrrdNew();
    double* array = (double*)calloc(size, sizeof(double));
    for (int i = 0 ; i < size ; ++i) {
        array[i] = RHS[i];
    }
    size_t __size[] = {dims[0], dims[1], dims[2]};
    double __spc[] = {domain.spacing()[0], domain.spacing()[1], domain.spacing()[2]};
    nrrdWrap_nva(nout, array, nrrdTypeDouble, 3, __size);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, __spc);
    nrrdSave("RHS.nrrd", nout, NULL);
    nrrdNuke(nout);
    
    gmm::ildlt_precond<gmm::csc_matrix<Type_> > PC(A);
    // // incomplete Cholesky preconditioner
    // gmm::mr_approx_inverse_precond<gmm::csc_matrix<Type_> > PC(A, 10, 1.0e-6);
    
    std::cerr << "PC set\n";
    
    gmm::cg(A, X, RHS, PC, iter);
    
    std::cerr << "system solved\n";
    
    // copy solution to output vector
    std::fill(x.begin(), x.end(), 0);
    for (int i = 0 ; i < size ; ++i) {
        ivec_type c = interior.coordinates(i);
        ++c[0];
        index_type j = domain.index(c);
        x[j] = X[i];
    }
}

} // xavier































































































