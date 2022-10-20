#ifndef __XAVIER_DIVERGENCE_CLEANING_HPP__
#define __XAVIER_DIVERGENCE_CLEANING_HPP__

#include <set>
#include <iostream>
#include <vector>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

#include <gmm/gmm_solver_cg.h>
#include <gmm/gmm_precond_ildlt.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_least_squares_cg.h>

#include <data/raster.hpp>
#include <data/grid.hpp>
#include <data/metric.hpp>


namespace div_cleaning {

typedef double                                         scalar_type;
typedef size_t                                         size_type;
typedef spurt::grid::uniform_grid<scalar_type, 3>     raster_type;
typedef raster_type::coord_type                        ivec_type;
typedef raster_type::vec_type                          vec_type;
typedef nvis::mat3                                     mat_type;
typedef spurt::raster_data<vec_type,3>                vdataset_type;
typedef spurt::raster_data<scalar_type, 3>            sdataset_type;

typedef gmm::wsvector<scalar_type>                     sparse_vector;
typedef gmm::row_matrix<sparse_vector>                 sparse_matrix;
typedef std::vector<scalar_type>                       dense_vector;

struct tetrahedron {
    vec_type p[4];
    vec_type v[4];
};

std::ostream& operator<<(std::ostream& os, const tetrahedron& T)
{
    os << "Tet:\n";
    for (int i = 0 ; i < 4 ; ++i) {
        os << i << ": p=" << T.p[i] << ", v=" << T.v[i] << '\n';
    }
    return os;
}

inline 
vec_type __interpolate(const vec_type& x, const tetrahedron& T, 
                       bool __verbose=false)
{
    // compute barycentric coordinates
    // b1 (p1-p0) + b2 (p2-p0) + b3 (p3-p0) = x - p0
    vec_type b, rhs;
    nvis::mat3 A;
    for (int i = 0 ; i < 3 ; ++i) {
        A(i, 0) = T.p[1][i] - T.p[0][i];
        A(i, 1) = T.p[2][i] - T.p[0][i];
        A(i, 2) = T.p[3][i] - T.p[0][i];
    }
    rhs = x - T.p[0];
    b = nvis::solve(A, rhs);
    if (__verbose) {
        std::cout << "__interpolate(): x = " << x << ", b = " << b << '\n';
    }
    
    // debug
    if (b[0] < 0 || b[1] < 0 || b[2] < 0 || 1 < b[0] - b[1] - b[2]) {
        throw std::runtime_error("invalid barycentric coordinates");
    }
    
    // interpolate
    vec_type v = (1 - b[0] - b[1] - b[2])*T.v[0] + b[0]*T.v[1] + b[1]*T.v[2] + b[2]*T.v[3];
    if (__verbose) {
        std::cout << "v = " << v << '\n';
    }
    return v;
}
} // div_cleaning


namespace {
    inline int even(int a)
    {
        return (a % 2 ? -1 : +1);
    }
    
    const int tets_id[64] = { -1, 18, -1, -1, 5, 6, 4, 7,
                              -1, -1, -1, -1, -1, -1, 10, -1,
                              13, 19, -1, -1, 14, -1, -1, -1,
                              12, -1, -1, -1, 15, -1, 9, -1,
                              -1, 17, -1, 23, -1, -1, -1, 22,
                              -1, -1, -1, 20, -1, -1, 11, 21,
                              -1, 16, -1, -1, -1, -1, -1, -1,
                              1, 2, 0, 3, -1, -1, 8, -1
                            };
}

namespace spurt {
template<typename Data>
struct divfree_field {
    typedef Data    data_type;
    typedef div_cleaning::raster_type grid_type;
    typedef div_cleaning::vec_type vec_type;
    typedef div_cleaning::ivec_type ivec_type;
    typedef div_cleaning::scalar_type scalar_type;
    typedef div_cleaning::tetrahedron tetrahedron;
    
    divfree_field(const grid_type& domain, const data_type& data)
        : _domain(domain), _data(data) {}
        
        
    int coords2id(const vec_type& c) const {
        scalar_type u = c[0], v = c[1], w = c[2];
        int code = 0;
        if (u <= v) {
            code |= 1;
        }
        if (u <= 1 - v) {
            code |= 2;
        }
        if (v <= w) {
            code |= 4;
        }
        if (v <= 1 - w) {
            code |= 8;
        }
        if (w <= u) {
            code |= 16;
        }
        if (w <= 1 - u) {
            code |= 32;
        }
        return code;
    }
    
    void verbose(bool v) const {
        __verbose = v;
    }
    
    vec_type voxel_value(const ivec_type& vid) const {
        // compute value associated with voxel center
        vec_type newval(0);
        
        for (int i = 0 ; i <= 1 ; ++i) {
            for (int j = 0 ; j <= 1 ; ++j) {
                for (int k = 0 ; k <= 1 ; ++k) {
                    const vec_type& fijk = _data(vid + ivec_type(i, j, k));
                    ivec_type ijk(i, j, k);
                    for (int d = 0 ; d < 3 ; ++d) {
                        newval[d] += fijk[d] +
                                     even(ijk[d] + ijk[(d + 1) % 3]) * _domain.spacing()[d] / _domain.spacing()[(d + 1) % 3] * fijk[(d+1)%3] +
                                     even(ijk[d] + ijk[(d + 2) % 3]) * _domain.spacing()[d] / _domain.spacing()[(d + 2) % 3] * fijk[(d+2)%3];
                    }
                }
            }
        }
        
        newval /= 8.;
        if (newval[2] < 0) {
            for (int u = 0 ; u <= 1 ; ++u) {
                for (int v = 0 ; v <= 1 ; ++v) {
                    for (int w = 0 ; w <= 1 ; ++w) {
                        ivec_type c = vid + ivec_type(u, v, w);
                        std::cerr << "v[" << c << "] = " << _data(c) << '\n';
                    }
                }
            }
            std::cerr << "v = " << newval << '\n';
            
            throw std::runtime_error("negative vz at voxel center");
        }
        return newval;
    }
    
    tetrahedron create_tetrahedron(const ivec_type& vid, int tid) const {
        int fid = tid / 4;
        ivec_type coords[4];
        for (int n = 0 ; n < 4 ; ++n) {
            coords[n] = grid::voxel_vertices[grid::voxel_faces[fid][n]];
        }
        vec_type fx = 0.5 * (_domain(coords[0]) + _domain(coords[2]));
        
        // compute value associated with face center
        vec_type fv(0);
        
        const int face2fixed[][2] = { {2, 0}, {2, 1}, {1, 0}, {0, 1}, {1, 1}, {0, 0} };
        int fixed = face2fixed[fid][0];
        int fixed_val = face2fixed[fid][1];
        int dof[] = {(fixed + 1) % 3, (fixed + 2) % 3 };
        for (int u = 0 ; u <= 1 ; ++u) {
            for (int v = 0 ; v <= 1 ; ++v) {
                ivec_type ijk;
                ijk[fixed] = fixed_val;
                ijk[dof[0]] = u;
                ijk[dof[1]] = v;
                fv[dof[0]] += _data(vid + ijk)[dof[0]] +
                              even(u + v) * _domain.spacing()[dof[0]] / _domain.spacing()[dof[1]] * _data(vid + ijk)[dof[1]];
                fv[dof[1]] += _data(vid + ijk)[dof[1]] +
                              even(u + v) * _domain.spacing()[dof[1]] / _domain.spacing()[dof[0]] * _data(vid + ijk)[dof[0]];
                fv[fixed] += _data(vid + ijk)[fixed];
            }
        }
        
        fv *= 0.25;
        if (fv[2] < 0) {
            ivec_type ijk;
            ijk[fixed] = fixed_val;
            for (int u = 0 ; u <= 1 ; ++u) {
                for (int v = 0 ; v <= 1 ; ++v) {
                    ijk[dof[0]] = u;
                    ijk[dof[1]] = v;
                    std::cerr << "ijk = " << ijk << '\n';
                    std::cerr << "data(" << vid << "+ijk) = " << _data(vid + ijk) << '\n';
                }
            }
            
            throw std::runtime_error("negative vz at face center");
        }
        
        vec_type cx = 0.5 * (_domain(ivec_type(0)) + _domain(ivec_type(1)));
        vec_type cv = voxel_value(vid);
        
        tetrahedron T;
        int ltid = tid % 4;
        T.p[0] = _domain(coords[ltid]);
        T.p[1] = _domain(coords[(ltid+1)%4]);
        T.p[2] = fx;
        T.p[3] = cx;
        T.v[0] = _data(vid + coords[ltid]);
        T.v[1] = _data(vid + coords[(ltid+1)%4]);
        T.v[2] = fv;
        T.v[3] = cv;
        
        return T;
    }
    
    vec_type operator()(const vec_type& x) const {
        vec_type actual_x = _domain.dmodulo(x);
        vec_type y = (actual_x - _domain(0, 0, 0)) / _domain.spacing();
        ivec_type vid(floor(y[0]), floor(y[1]), floor(y[2]));
        vec_type local = (y - vec_type(vid)) * _domain.spacing();
        int tet = tets_id[coords2id(local)];
        tetrahedron T = create_tetrahedron(vid, tet);
        
        if (__verbose) {
            std::cout << "divfree_field::operator(): x = " << x << ", actual_x = " << actual_x
                      << ", y = " << y << ", vid = " << vid << ", local = " << local << ", tet = " << tet
                      << ", T = " << T;
        }
        
        return __interpolate(local, T, __verbose);
    }
    
    vec_type interpolate(const vec_type& x) const {
        return (*this)(x);
    }
    
    const grid_type& get_grid() const {
        return _domain;
    }
    
    const data_type& get_data() const {
        return _data;
    }
    
    grid_type        _domain;
    const data_type& _data;
    mutable bool     __verbose;
};

inline void divergence(div_cleaning::sdataset_type& out,
                       const div_cleaning::vdataset_type& in)
{
    using namespace div_cleaning;
    
    div_cleaning::raster_type domain(out.grid());
    ivec_type dims = domain.resolution();
    
    static const ivec_type dc[] = { ivec_type(1, 0, 0),
                                    ivec_type(0, 1, 0),
                                    ivec_type(0, 0, 1)
                                  };
                                  
    // compute divergence for all interior vertices
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 0 ; i < dims[0] ; ++i) {
                ivec_type c(i, j, k);
                if (domain.on_boundary(c)) {
                    out(c) = 0.;
                } else {
                    scalar_type div = 0;
                    for (int dim = 0 ; dim < 3 ; ++dim) {
                        div += 0.5 * (in(c + dc[dim])[dim] - in(c - dc[dim])[dim]) / domain.spacing()[dim];
                    }
                    out(c) = div;
                }
            }
        }
    }
}

inline void gradient(div_cleaning::vdataset_type& out, const div_cleaning::sdataset_type& in)
{
    using namespace div_cleaning;
    
    div_cleaning::raster_type domain(out.grid());
    ivec_type dims = domain.resolution();
    
    static const ivec_type dc[] = { ivec_type(1, 0, 0),
                                    ivec_type(0, 1, 0),
                                    ivec_type(0, 0, 1)
                                  };
                                  
    // compute gradient for all interior vertices
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 0 ; i < dims[0] ; ++i) {
                ivec_type c(i, j, k);
                if (domain.on_boundary(c)) {
                    out(c) = vec_type(0);
                } else {
                    vec_type g(0);
                    for (int dim = 0 ; dim < 3 ; ++dim) {
                        g[dim] = 0.5 * (in(c + dc[dim]) - in(c - dc[dim])) / domain.spacing()[dim];
                    }
                    out(c) = g;
                }
            }
        }
    }
}

inline void divergence(std::vector<div_cleaning::scalar_type>& out,
                       const std::vector<div_cleaning::scalar_type>& in,
                       const div_cleaning::raster_type& domain)
{
    using namespace div_cleaning;
    
    ivec_type dims = domain.resolution();
    
    // compute divergence for all vertices
    std::fill(out.begin(), out.end(), 0.);
    for (int k = 1 ; k < dims[2] - 1 ; ++k) {
        for (int j = 1 ; j < dims[1] - 1 ; ++j) {
            for (int i = 1 ; i < dims[0] - 1 ; ++i) {
                ivec_type c(i, j, k);
                scalar_type div = 0;
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
                    
                    size_type hi = 3 * domain.index(__hi) + d;
                    size_type lo = 3 * domain.index(__lo) + d;
                    
                    div += (in[hi] - in[lo]) / (2.*domain.spacing()[d]);
                }
                out[domain.index(c)] = div;
            }
        }
    }
}

inline div_cleaning::vec_type gradient(const std::vector<div_cleaning::scalar_type>& vals,
                                       const div_cleaning::raster_type& domain,
                                       const div_cleaning::ivec_type& c)
{
    using namespace div_cleaning;
    
    nvis::vec3 g(0, 0, 0);
    ivec_type dims = domain.resolution();
    // int id = domain.index(c);
    for (int i = 0 ; i < 3 ; ++i) {
        nvis::ivec3 hi(c), lo(c);
        ++hi[i];
        --lo[i];
        scalar_type h = 2.*domain.spacing()[i];
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
        
        size_type __lo = domain.index(lo);
        size_type __hi = domain.index(hi);
        g[i] = (vals[__hi] - vals[__lo]) / h;
    }
    return g;
}

inline void gradient(std::vector<div_cleaning::scalar_type>& out,
                     const std::vector<div_cleaning::scalar_type>& in,
                     const div_cleaning::raster_type& domain)
{
    using namespace div_cleaning;
    
    ivec_type dims = domain.resolution();
    // compute gradient of input scalar field at all interior vertices
    std::fill(out.begin(), out.end(), 0.);
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 0 ; i < dims[0] ; ++i) {
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
                 const div_cleaning::sparse_matrix& A,
                 div_cleaning::dense_vector& x,
                 const div_cleaning::dense_vector& b,
                 const PC& P,
                 gmm::iteration& iter)
{
    using namespace div_cleaning;
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
        // gmm::add(q, gmm::scaled(c, mu*gmm::vect_sp(c, p)), q);
        
        a = rho / gmm::vect_sp(q, p);
        gmm::add(gmm::scaled(p,  a), x);
        gmm::add(gmm::scaled(q, -a), r);
        
        rho_1 = rho;
        ++iter;
    }
}


} // namespace spurt

#endif // __XAVIER_DIVERGENCE_CLEANING_HPP__







































