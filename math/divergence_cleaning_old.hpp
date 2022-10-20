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

#include <math/poisson_gmm.hpp>
#include <data/grid.hpp>

typedef spurt::grid<double, 3>     grid_type;
typedef grid_type::ivec_type        ivec_type;
typedef grid_type::index_type       index_type;
typedef grid_type::vec_type         vec_type;
typedef nvis::mat3                  mat_type;

typedef gmm::wsvector<double>           sparse_vector;
typedef gmm::row_matrix<sparse_vector>  sparse_matrix;
typedef std::vector<double>             dense_vector;

namespace {
struct tetrahedron {
    vec_type p[4];
    vec_type v[4];
};

inline vec_type interpolate(const vec_type& x, const tetrahedron& T)
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
    
    // interpolate
    return (1 - b[0] - b[1] - b[2])*T.v[0] + b[0]*T.v[1] + b[1]*T.v[2] + b[2]*T.v[3];
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
} // anonymous


namespace spurt {
inline int even(int a)
{
    return (a % 2 ? -1 : +1);
}

template<typename Data>
struct divfree_field {
    typedef Data    data_type;
    
    divfree_field(const grid_type& domain, const data_type& data)
        : _domain(domain), _data(data) {}
        
        
    int coords2id(const vec_type& c) const {
        double u = c[0], v = c[1], w = c[2];
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
    
    std::pair<vec_type, vec_type> voxel_center(const ivec_type& vid) const {
        // create voxel vertex
        vec_type vx = 0.5 * (_domain(vid) + _domain(vid + ivec_type(1, 1, 1)));
        // compute value associated with voxel center
        vec_type newval(0);
        for (int u = 0 ; u <= 1 ; ++u) {
            for (int v = 0 ; v <= 1 ; ++v) {
                for (int w = 0 ; w <= 1 ; ++w) {
                    ivec_type c = vid + ivec_type(u, v, w);
                    for (int d = 0 ; d < 3 ; ++d) {
                        newval[d] += _data(c)[d] +
                                     even(d + (d + 1) % 3) * _domain.step(d) / _domain.step((d + 1) % 3) * _data(c)[(d+1)%3] +
                                     even(d + (d + 2) % 3) * _domain.step(d) / _domain.step((d + 2) % 3) * _data(c)[(d+2)%3];
                    }
                }
            }
        }
        newval /= 8.;
        return std::make_pair(vx, newval);
    }
    
    tetrahedron create_tetrahedron(const ivec_type& vid, int tid) const {
        int fid = tid / 4;
        nvis::ivec4 face;
        ivec_type coords[4];
        for (int n = 0 ; n < 4 ; ++n) {
            coords[n] = vid + voxel_vertices[voxel_faces[fid][n]];
            face[n] = _domain.index(coords[n]);
        }
        vec_type fx = 0.5 * (_domain(coords[0]) + _domain(coords[2]));
        
        // compute value associated with face center
        vec_type fv(0);
        const int fid2fixed[] = { 2, 2, 1, 0, 1, 0 };
        const int uv2id[] = { 0, 1, 3, 2 };
        int fixed = fid2fixed[fid];
        int dof[2] = {(fixed + 1) % 3, (fixed + 2) % 3};
        for (int n = 0 ; n < 2 ; ++n) {
            fv[dof[n]] = 0.;
            for (int u = 0; u <= 1 ; ++u) {
                for (int v = 0 ; v <= 1 ; ++v) {
                    int id = uv2id[u+2*v];
                    fv[dof[n]] += _data(coords[id])[dof[n]] +
                                  even(u + v) * _domain.step(dof[n]) / _domain.step(dof[(n+1)%2]) * _data(id)[dof[(n+1)%2]];
                }
            }
        }
        fv[fixed] = 0.;
        for (int n = 0 ; n < 4 ; ++n) {
            fv[fixed] += _data(face[n])[fixed];
        }
        fv *= 0.25;
        
        std::pair<vec_type, vec_type> c = voxel_center(vid);
        vec_type cx = c.first();
        vec_type cv = c.second();
        
        tetrahedron T;
        int ltid = tid % 4;
        T.p[0] = _domain(coords[ltid]);
        T.p[1] = _domain(coords[(ltid+1)%4]);
        T.p[2] = fx;
        T.p[3] = cx;
        T.v[0] = _data(coords[ltid]);
        T.v[1] = _data(coords[(ltid+1)%4]);
        T.v[2] = fv;
        T.v[3] = cv;
        return T;
    }
    
    vec_type operator()(const vec_type& x) const {
        vec_type y = (x - _domain(0, 0, 0)) / _domain.spacing();
        ivec_type vid(floor(y[0]), floor(y[1]), floor(y[2]));
        vec_type local = (y - vec_type(vid));
        int tet = tets_id[coords2id(local)];
        tetrahedron T = create_tetrahedron(vid, tet);
        return interpolate(x, T);
    }
    
    grid_type           _domain;
    const data_type&    _data;
};

// implementation of Peikert & Sadlo's divergence free interpolation method
// given a piecewise divergence free rectilinear mesh
void divfree(std::vector<nvis::vec3>& vertices, std::vector<nvis::vec3>& values,
             std::vector< nvis::fixed_vector<int, 4> >& tets,
             const std::vector<double>& flow, const grid_type& domain)
{
    // export divergence free piecewise linear field
    const nvis::ivec3 dims = domain.dimensions();
    
    int nbvertices = dims[0] * dims[1] * dims[2];
    int nbvoxels = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
    int nbfaces = (dims[0] - 1) * (dims[1] - 1) * dims[2] +
                  (dims[0] - 1) * dims[1] * (dims[2] - 1) +
                  dims[0] * (dims[1] - 1) * (dims[2] - 1);
                  
    vertices.reserve(nbvertices + nbfaces + nbvoxels);
    values.reserve(nbvertices + nbfaces + nbvoxels);
    vertices.clear();
    values.clear();
    for (int k = 0 ; k < dims[2] ; ++k) {
        for (int j = 0 ; j < dims[1] ; ++j) {
            for (int i = 0 ; i < dims[0] ; ++i) {
                int id = domain.index(i, j, k);
                vertices.push_back(domain(i, j, k));
                values.push_back(flow[id]);
            }
        }
    }
    
    std::map < id4, int, Lt_id_set > added;
    tets.reserve(24*nbvoxels);
    for (int k = 0 ; k < dims[2] - 1 ; ++k) {
        for (int j = 0 ; j < dims[1] - 1 ; ++j) {
            for (int i = 0 ; i < dims[0] - 1 ; ++i) {
            
                nvis::ivec3 base_coord(i, j, k);
                
                // create new voxel vertex
                nvis::vec3 vx = 0.5 * (domain(i, j, k) + domain(i + 1, j + 1, k + 1));
                vertices.push_back(vx);
                int vid = vertices.size() - 1;
                
                // compute value associated with voxel center
                nvis::vec3 newval(0);
                for (int u = 0 ; u <= 1 ; ++u) {
                    for (int v = 0 ; v <= 1 ; ++v) {
                        for (int w = 0 ; w <= 1 ; ++w) {
                            int ijk = domain.index(base_coord + nvis::ivec3(u, v, w));
                            for (int d = 0 ; d < 3 ; ++d) {
                                newval[d] += flow[3*ijk+d] +
                                             even(d + (d + 1) % 3) * domain.step(d) / domain.step((d + 1) % 3) * flow[3*ijk+(d+1)%3] +
                                             even(d + (d + 2) % 3) * domain.step(d) / domain.step((d + 2) % 3) * flow[3*ijk+(d+2)%3];
                            }
                        }
                    }
                }
                newval /= 8.;
                values.push_back(newval);
                
                // create face vertices
                for (int f = 0 ; f < 6 ; ++f) {
                    int ids[4];
                    nvis::ivec3 coords[4];
                    for (int n = 0 ; n < 4 ; ++n) {
                        coords[n] = base_coord + voxel_vertices[voxel_faces[f][n]];
                        ids[n] = domain.index(coords[n]);
                    }
                    id4 face(ids[0], ids[1], ids[2], ids[3]);
                    int fid;
                    if (added.find(face) == added.end()) {
                        nvis::vec3 fx = 0.5 * (domain(coords[0]) + domain(coords[2]));
                        vertices.push_back(fx);
                        added[face] = vertices.size() - 1;
                        fid = vertices.size() - 1;
                        
                        // compute value associated with face center
                        nvis::vec3 newval;
                        int fixed = f / 2;
                        int dof[2] = {(fixed + 1) % 3, (fixed + 2) % 3};
                        for (int n = 0 ; n < 2 ; ++n) {
                            newval[dof[n]] = 0.;
                            for (int u = 0; u <= 1 ; ++u) {
                                for (int v = 0 ; v <= 1 ; ++v) {
                                    int ijk = ids[u+2*v];
                                    newval[dof[n]] += flow[3*ijk+dof[n]] +
                                                      even(u + v) * domain.step(dof[n]) / domain.step(dof[(n+1)%2]) * flow[3*ijk+dof[(n+1)%2]];
                                }
                            }
                        }
                        newval[fixed] = 0.;
                        for (int n = 0 ; n < 4 ; ++n) {
                            newval[fixed] += flow[3*ids[n] + fixed];
                        }
                        newval *= 0.25;
                        values.push_back(newval);
                    } else {
                        fid = added[face];
                    }
                    
                    // add four tetrahedra to the mesh
                    for (int e = 0 ; e < 4 ; ++e) {
                        int id0 = ids[e];
                        int id1 = ids[(e+1)%4];
                        tets.push_back(nvis::fixed_vector<int, 4>(id0, id1, fid, vid));
                    }
                }
            }
        }
    }
}

inline double divergence(const ivec_type& c, const std::vector<double>& in,
                         const grid_type& domain)
{
    double div = 0;
    for (int i = 0 ; i < 3 ; ++i) {
        ivec_type step(0, 0, 0);
        step[i] = 1;
        int hi = 3 * domain.index(c + step) + i;
        int lo = 3 * domain.index(c - step) + i;
        div += 0.5 / domain.step(i) * (in[hi] - in[lo]);
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
        g[i] = 0.5 / domain.step(i) * (in[hi] - in[lo]);
    }
    return g;
}

inline void divergence(std::vector<double>& out,
                       const std::vector<double>& in,
                       const grid_type& domain)
{
    ivec_type dims = domain.dimensions();
    
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
    ivec_type dims = domain.dimensions();
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

} // namespace spurt

#endif // __XAVIER_DIVERGENCE_CLEANING_HPP__






