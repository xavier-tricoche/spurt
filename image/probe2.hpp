#ifndef __PROBE_HPP__
#define __PROBE_HPP__

#include <teem/gage.h>
#include <teem/nrrd.h>
#include <teem/ten.h>
#include <vector>
#include <iostream>
#include <assert.h>

#include <Eigen/Eigen>

typedef Eigen::Vector<double, 6> vec6;
typedef Eigen::Vector<double, 4> vec4;
typedef Eigen::Vector<double, 3> vec3;
typedef Eigen::Vector<double, 2> vec2;
typedef Eigen::Matrix<double, 3, 3> mat3;
typedef Eigen::Matrix<double, 2, 2> mat2;

namespace spurt {
struct coord_mapping {
	
	coord_mapping() : origin(0,0,0), step(0,0,0) {}

    coord_mapping(const Nrrd* nrrd, int N) : origin(0, 0, 0), step(0, 0, 0) {
        for (int i = 0 ; i < N ; ++i) {
            origin[i] = nrrd->axis[nrrd->dim-N+i].min;
            if (std::isnan(origin[i])) origin[i] = 0.;
            step[i] = nrrd->axis[nrrd->dim-N+i].spacing;
            if (std::isnan(step[i])) step[i] = 1.;
        }
    }

    vec3 to_world(const vec3& idx) const {
        return origin + idx*step;
    }

    vec3 to_index(const vec3& wc) const {
        return (wc - origin).array() / step.array();
    }

    vec3 origin, step;
};


// general remark about interpolation in gage_interface:
// point coordinates MUST be specified with respect to raster coordinates,
// i.e. [0, N-1]
namespace gage_interface {
    
// supporting kernels
typedef enum {
    TENT,
    BC_BLUR,
    BC_INTERP,
    BSPL3,
    BSPL3_INTERP,
    BSPL5,
    BSPL5_INTERP,
    BSPL7,
    BSPL7_INTERP,
    C3_QUINTIC,
    C4_HEXIC,
    GAUSS,
    HANN,
    BLACK
} kernel_idx_t; 

typedef double gage_t;

struct index_grid {
    index_grid(const Nrrd* nrrd, unsigned int upsampling = 1);

    unsigned int id(unsigned int i, unsigned int j, unsigned int k) const;
    vec3 operator()(unsigned int i, unsigned int j, unsigned int k) const;

    double min[3], max[3];
    unsigned int size[3];
    double d[3];
};

int set_kernel(gageContext* ctx, kernel_idx_t);
bool gage_probe(gageContext *ctx, double u, double v, double w, 
                const coord_mapping& wcidx, bool use_wc);

// gage interface for scalar volumes
class scalar_wrapper
{
public:

    scalar_wrapper(const Nrrd* nrrd, kernel_idx_t kidx=BC_BLUR,
                   bool use_grad=true, bool use_hess=true, 
                   bool use_heval=true, bool use_hevec=true);
    ~scalar_wrapper();

    void use_world() {
        _use_wc = true;
    }

    bool value(const vec2& u, double& v) const;   // 2D
    bool value(const vec3& u, double& v) const;   // 3D

    bool gradient(const vec2& u, vec2& g) const;   // 2D
    bool gradient(const vec3& u, vec3& g) const;   // 3D

    bool hessian(const vec2& u, vec3& h) const;   // 2D
    bool hessian(const vec3& u, vec6& h) const;   // 3D
    bool hessian(const vec2& u, mat2& h) const;   // 2D
    bool hessian(const vec3& u, mat3& h) const;   // 3D

    bool hess_evals(const vec2 &u, vec2 &vals) const;   // 3D
    bool hess_evecs(const vec2 &u, mat2 &vecs) const;   // 3D
    bool hess_evals(const vec3& u, vec3& vals) const;   // 3D
    bool hess_evecs(const vec3 &u, mat3 &vecs) const;   // 3D
    
    unsigned int support_radius() const { return ctx->radius; }

    const gageContext *get_ctx() const { return ctx; }

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_v, *_g, *_h, *_heval, *_hevec;
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    bool _use_grad, _use_hess, _use_heval, _use_hevec; 
    spurt::coord_mapping wcidx;
};

// gage interface for vector volumes
class vector_wrapper
{
public:
	vector_wrapper() : ctx(NULL), pv(NULL), _v(NULL), _j(NULL) {}
    vector_wrapper(const Nrrd* nrrd, kernel_idx_t kidx=BC_BLUR,
                   bool _use_jac=true);
	vector_wrapper(const vector_wrapper& other);
    ~vector_wrapper();

    void use_world() {
        _use_wc = true;
    }
    bool value(const vec2& u, vec2& v) const;   // 2D
    bool value(const vec3& u, vec3& v) const;   // 3D

    bool jacobian(const vec2& u, mat2& g) const;   // 2D
    bool jacobian(const vec3& u, mat3& g) const;   // 3D

    unsigned int support_radius() const { return ctx->radius; }

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_v, *_j;
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    bool _use_jac;
    spurt::coord_mapping wcidx;
	kernel_idx_t m_kernel_idx;
};

// gage interface for diffusion tensor volumes
class diffusion_tensor_wrapper {
public:
    diffusion_tensor_wrapper(const Nrrd *nrrd, kernel_idx_t kidx = BC_BLUR,
                   bool use_evals = true, bool use_evecs = true,
                   bool use_aniso = true, bool use_fa = true,
                   bool use_mode = true,
                   bool use_fa_ridge = true, 
                   bool use_mode_ridge = true);
    ~diffusion_tensor_wrapper();

    void use_world() {
        _use_wc = true;
    }

    bool value(const vec3 &p, mat3 &T) const;
    bool confidence(const vec3 &p, double &value) const;
    bool aniso(const vec3 &p, int anisokind, double &val) const;

    bool evals(const vec3 &p, vec3 &val) const;
    bool evecs(const vec3 &p, mat3 &val) const;
    
    bool fa(const vec3 &p, double &val) const;
    bool fa_grad(const vec3 &p, vec3 &grad) const;
    bool fa_hess(const vec3 &p, mat3 &hess) const;
    bool fa_hess_evals(const vec3 &p, vec3 &eval) const;
    bool fa_hess_evecs(const vec3 &p, mat3 &evec) const;
    bool fa_ridge_surf_strength(const vec3 &p, double &val) const;
    bool fa_valley_surf_strength(const vec3 &p, double &val) const;

    bool mode(const vec3& u, double& v) const;
    bool mode_grad(const vec3& u, vec3& g) const;
    bool mode_hess(const vec3& u, mat3& h) const;
    bool mode_hess_evals(const vec3& u, vec3& vals) const;
    bool mode_hess_evecs(const vec3 &u, mat3 &vecs) const;

    unsigned int support_radius() const { return ctx->radius; }

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_t, *_c, *_aniso, *_eval, *_evec;
    const gage_t *_fa, *_fa_grad, *_fa_hess, *_fa_hess_eval, *_fa_hess_evec,
    *_fa_ridge_surf_strength, *_fa_valley_surf_strength;
    const gage_t *_mode, *_mode_grad, *_mode_hess, *_mode_hess_eval, *_mode_hess_evec;
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    bool _use_evals;
    bool _use_evecs;
    bool _use_aniso;
    bool _use_fa;
    bool _use_mode;
    bool _use_fa_ridge;
    bool _use_mode_ridge;
    spurt::coord_mapping wcidx;
};

} // gage_interface

} // spurt

inline
spurt::gage_interface::index_grid::
index_grid(const Nrrd* nrrd, unsigned int upsample)
{
//    NrrdAxisInfo *info;
    assert(nrrd->dim >= 3);
    for (unsigned int i = 0 ; i < 3 ; i++) {
        const NrrdAxisInfo& axis = nrrd->axis[nrrd->dim-3+i];

        // raster coordinates
        min[i] = 0.;
        max[i] = (double)axis.size - 1.0;
        if (axis.center == nrrdCenterCell) {
            min[i] -= 0.499;
            max[i] += 0.499;
        }
        size[i] = axis.size * upsample;
        d[i] = (max[i] - min[i]) / (size[i] - 1);
    }
}

inline
unsigned int
spurt::gage_interface::index_grid::
id(unsigned int i, unsigned int j, unsigned int k) const
{
    return i + size[0]*(j + size[1]*k);
}

inline
vec3
spurt::gage_interface::index_grid::
operator()(unsigned int i, unsigned int j, unsigned int k) const
{
    return vec3(min[0] + i*d[0],
                min[1] + j*d[1],
                min[2] + k*d[2]);
}

inline bool
spurt::gage_interface::
gage_probe(gageContext *ctx, double u, double v, double w, 
           const coord_mapping &wcidx, bool use_wc) {
    if (!use_wc) return gageProbe(ctx, u, v, w);
    vec3 idx = wcidx.to_index(vec3(u, v, w));
    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}

#endif












