#ifndef __PROBE_HPP__
#define __PROBE_HPP__

#include <teem/gage.h>
#include <teem/nrrd.h>
#include <teem/ten.h>
#include <vector>
#include <iostream>
#include <math/fixed_vector.hpp>
#include <assert.h>

#include <Eigen/Core>

namespace nvis {
typedef fixed_vector< double, 6 > vec6;
typedef fixed_vector< double, 4 > vec4;
}

namespace xavier {
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

    nvis::vec3 to_world(const nvis::vec3& idx) const {
        return origin + idx*step;
    }

    nvis::vec3 to_index(const nvis::vec3& wc) const {
        return (wc - origin) / step;
    }

    nvis::vec3 origin, step;
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
    nvis::vec3 operator()(unsigned int i, unsigned int j, unsigned int k) const;

    double min[3], max[3];
    unsigned int size[3];
    double d[3];
};

int set_kernel(gageContext* ctx, kernel_idx_t);

// gage interface for scalar volumes
class scalar_wrapper {
public:
    scalar_wrapper(const Nrrd* nrrd, kernel_idx_t kidx=BC_BLUR,
                   bool use_grad=true, bool use_hess=true, 
                   bool use_heval=true, bool use_hevec=true);
    ~scalar_wrapper();

    void use_world() {
        _use_wc = true;
    }

    bool value(double* u, double& v) const;   // 2D
    bool value(const nvis::vec1& u, double& v) const;   // 1D
    bool value(const nvis::vec2& u, double& v) const;   // 2D
    bool value(const nvis::vec3& u, double& v) const;   // 3D

    template<typename Position_, size_t N, typename Value_=double>
    bool value(const Position_& u, Value_& v) const;

    bool gradient(double* u, double* g) const;   // 2D
    bool gradient(const nvis::vec2& u, nvis::vec2& g) const;   // 2D
    bool gradient(const nvis::vec3& u, nvis::vec3& g) const;   // 3D
    
    template<typename Position_, size_t N, typename Vector_ = Position_>
    bool gradient(const Position_& u, Vector_& g) const; 

    bool hessian(double* u, double* h) const;   // 2D
    bool hessian(const nvis::vec2& u, nvis::vec3& h) const;   // 2D
    bool hessian(const nvis::vec3& u, nvis::vec6& h) const;   // 3D
    
    template<typename Position_, size_t N, typename Matrix_>
    bool hessian(const Position_& u, Matrix_& h) const;
    
    bool hess_evals(const nvis::vec3& u, nvis::vec3& vals) const;   // 3D
    bool hess_evecs(const nvis::vec3& u, std::vector< nvis::vec3 >& vecs) const;   // 3D

    template<typename Position_, size_t N, typename Vector_=Position_>
    bool hess_evals(const Position_& u, Vector_& vals) const;
    
    template<typename Position_, size_t N, typename Vector_=Position_>
    bool hess_evecs(const Position_& u, std::vector<Vector_>& evecs) const;
    
    unsigned int support_radius() const { return ctx->radius; }

    bool __gageProbe(gageContext* ctx, double u, double v, double w) const;

    const gageContext* get_ctx() const { return ctx; }

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_v, *_g, *_h, *_heval, *_hevec;
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    bool _use_grad, _use_hess, _use_heval, _use_hevec; 
    xavier::coord_mapping wcidx;
};

// gage interface for scalar volumes
class vector_wrapper
{
public:
    typedef nvis::fixed_vector<double, 2> point2_t;
    typedef nvis::fixed_vector<double, 2> value2_t;
    typedef nvis::fixed_vector<double, 3> point3_t;
    typedef nvis::fixed_vector<double, 3> value3_t;
    typedef Eigen::Matrix<double, 2, 2> deriv2_t;
    typedef Eigen::Matrix<double, 3, 3> deriv3_t;
    
	vector_wrapper() : ctx(NULL), pv(NULL), _v(NULL), _j(NULL) {}
    vector_wrapper(const Nrrd* nrrd, kernel_idx_t kidx=BC_BLUR,
                   bool _use_jac=true);
	vector_wrapper(const vector_wrapper& other);
    ~vector_wrapper();

    void use_world() {
        _use_wc = true;
    }

    bool value(const point2_t& u, value2_t& v) const;   // 2D
    bool value(const point3_t& u, value3_t& v) const;   // 3D

    bool jacobian(const point2_t& u, deriv2_t& g) const;   // 2D
    bool jacobian(const point3_t& u, deriv3_t& g) const;   // 3D
    
    unsigned int support_radius() const { return ctx->radius; }

    bool __gageProbe(gageContext* ctx, double u, double v, double w) const;

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_v, *_j;
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    bool _use_jac;
    xavier::coord_mapping wcidx;
	kernel_idx_t m_kernel_idx;
};

// gage interface for DTI volumes
class tensor_wrapper {
public:
    tensor_wrapper(const Nrrd* nrrd, kernel_idx_t kidx);
    ~tensor_wrapper();

    void use_world() {
        _use_wc = true;
    }

    bool tensor(double *u, double* t) const;
    bool confidence(double *u, double& c) const;
    bool aniso(double *u, int anisokind, double& fa) const;
    bool evals(double* u, double* eval) const;
    bool evecs(double* u, double* evec) const;

    bool confidence(const nvis::vec3& u, double& v) const;

    bool fa(double *u, double *fa) const;
    bool fa_grad(double *u, double *grad) const;
    bool fa_hess(double *u, double *hess) const;
    bool fa_hess_evals(double *u, double *eval) const;
    bool fa_hess_evecs(double *u, double *evec) const;
    bool fa_ridge_surf_strength(double *u, double *val) const;
    bool fa_valley_surf_strength(double *u, double *val) const;

    bool fa(const nvis::vec3& u, double& v) const;
    bool fa_grad(const nvis::vec3& u, nvis::vec3& g) const;
    bool fa_hess(const nvis::vec3& u, nvis::vec6& h) const;
    bool fa_hess_evals(const nvis::vec3& u, nvis::vec3& vals) const;
    bool fa_hess_evecs(const nvis::vec3& u, std::vector< nvis::vec3 >& vecs) const;

    bool mode(double *u, double *mode) const;
    bool mode_grad(double *u, double *grad) const;
    bool mode_hess(double *u, double *hess) const;
    bool mode_hess_evals(double *u, double *eval) const;
    bool mode_hess_evecs(double *u, double *evec) const;

    bool mode(const nvis::vec3& u, double& v) const;
    bool mode_grad(const nvis::vec3& u, nvis::vec3& g) const;
    bool mode_hess(const nvis::vec3& u, nvis::vec6& h) const;
    bool mode_hess_evals(const nvis::vec3& u, nvis::vec3& vals) const;
    bool mode_hess_evecs(const nvis::vec3& u, std::vector< nvis::vec3 >& vecs) const;
    
    unsigned int support_radius() const { return ctx->radius; }

    bool __gageProbe(gageContext* ctx, double u, double v, double w) const;

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_t, *_c, *_aniso, *_eval, *_evec;
    const gage_t *_fa, *_fa_grad, *_fa_hess, *_fa_hess_eval, *_fa_hess_evec,
    *_fa_ridge_surf_strength, *_fa_valley_surf_strength;
    const gage_t *_mode, *_mode_grad, *_mode_hess, *_mode_hess_eval, *_mode_hess_evec;
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    xavier::coord_mapping wcidx;
};


#if 0
template<typename Scalar_, unsigned int Dim_, unsigned int Order_, 
         typename Value_=std::array<Scalar_,Dim_>, 
         typename Deriv_=std::array<Scalar_,void, typename SecondDeriv_=void>
class smooth_wrapper {
    typedef Scalar_ scalar_type;
    constexpr unsigned int dimension = Dim_;
    constexpr unsigned int order = Order_;
    
    smooth_wrapper(const Nrrd* nrrd, kernel_idx_t kidx=BC_BLUR,
                   bool use_1st_deriv=true, bool use_2nd_deriv=false);
    ~smooth_wrapper();
    
    /*
    order O
    dim D
    D: 1 s: 1, v: 1, t: 1
    D: 2 s: 1, v: 2, t: 4 / 3 / 4
    D: 3 s: 1, v: 3. t: 9 / 6 / 7
    
    */

    void use_world() {
        _use_wc = true;
    }

    bool value(double* u, double* v) const;

    bool derivative(double* u, double* d) const;
    
    bool second_derivative(double* u, double* dd) const;
    
    unsigned int support_radius() const { return ctx->radius; }

    bool __gageProbe(gageContext* ctx, double u, double v, double w) const;

    const gageContext* get_ctx() const { return ctx; }

protected:
    gageContext *ctx;
    gagePerVolume *pv;
    const gage_t *_v, *_g, *_h;
    const gage_t *_t, *_c, *_aniso, *_eval, *_evec;
    const gage_t *_fa, *_fa_grad, *_fa_hess, *_fa_hess_eval, *_fa_hess_evec,
    *_fa_ridge_surf_strength, *_fa_valley_surf_strength;
    const gage_t *_mode, *_mode_grad, *_mode_hess, *_mode_hess_eval, *_mode_hess_evec;
    const gage_t *_v, *_j;
    
    double kparm[NRRD_KERNEL_PARMS_NUM];
    bool _use_wc;
    bool _use_grad, _use_hess, _use_heval, _use_hevec; 
    xavier::coord_mapping wcidx;
    
};
#endif

};

};

inline
xavier::gage_interface::index_grid::
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
xavier::gage_interface::index_grid::
id(unsigned int i, unsigned int j, unsigned int k) const
{
    return i + size[0]*(j + size[1]*k);
}

inline
nvis::vec3
xavier::gage_interface::index_grid::
operator()(unsigned int i, unsigned int j, unsigned int k) const
{
    return nvis::vec3(min[0] + i*d[0],
                      min[1] + j*d[1],
                      min[2] + k*d[2]);
}

inline bool xavier::gage_interface::scalar_wrapper::
__gageProbe(gageContext* ctx, double u, double v, double w) const
{
    if (!_use_wc) return gageProbe(ctx, u, v, w);
    nvis::vec3 idx = wcidx.to_index(nvis::vec3(u, v, w));
    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}


template<typename Position_, size_t N, typename Value_>
bool xavier::gage_interface::scalar_wrapper::
value(const Position_& u, Value_& v) const
{
    bool asw;
    if (N == 2) asw = __gageProbe(ctx, u[0], u[1], 0);
    else if (N==3) asw = __gageProbe(ctx, u[0], u[1], u[2]);
    else throw std::runtime_error("Invalid spatial dimension: " + std::to_string(N));
    if (asw) return false;
    v = *_v;
    return true;
}

template<typename Position_, size_t N, typename Vector_>
bool xavier::gage_interface::scalar_wrapper::
gradient(const Position_& u, Vector_& g) const
{   
    if (!_use_grad) {
        throw std::runtime_error("Gradient computation deactivated");
    }

    bool asw;
    if (N == 2) asw = __gageProbe(ctx, u[0], u[1], 0);
    else if (N==3) asw = __gageProbe(ctx, u[0], u[1], u[2]);
    else throw std::runtime_error("Invalid spatial dimension: " + std::to_string(N));
    if (asw) return false;
    g[0] = _g[0];
    g[1] = _g[1];
    if (N==3) g[2] = _g[2];
    return true;
}

template<typename Position_, size_t N, typename Matrix_>
bool xavier::gage_interface::scalar_wrapper::
hessian(const Position_& u, Matrix_& h) const
{   
    if (!_use_hess) {
        throw std::runtime_error("Hessian computation deactivated");
    }

    bool asw;
    if (N == 2) asw = __gageProbe(ctx, u[0], u[1], 0);
    else if (N==3) asw = __gageProbe(ctx, u[0], u[1], u[2]);
    else throw std::runtime_error("Invalid spatial dimension: " + std::to_string(N));
    if (asw) return false;
    /*
    0 1 2
    3 4 5
    6 7 8
    */
    h(0,0) = _h[0];
    h(0,1) = h(1,0) = _h[1];
    h(1,1) = _h[4];

    if (N==3) {
        h(0,2) = h(2,0) = _h[2];
        h(1,2) = h(2,1) = _h[5];
        h(2,2) = _h[8];
    }
    return true;
}

template<typename Position_, size_t N, typename Vector_>
bool xavier::gage_interface::scalar_wrapper::
hess_evals(const Position_& u, Vector_& vals) const {
    if (!_use_heval) {
        throw std::runtime_error("Hessian eigenvalues computation deactivated");
    }

    bool asw;
    if (N==2) asw = __gageProbe(ctx, u[0], u[1], 0);
    else asw =  __gageProbe(ctx, u[0], u[1], u[2]);
    if (asw) return false;
    
    vals[0] = _heval[0];
    vals[1] = _heval[1];
    if (N==3) vals[2] = _heval[2];
    return true;
}

template<typename Position_, size_t N, typename Vector_>
bool xavier::gage_interface::scalar_wrapper::
hess_evecs(const Position_& u, std::vector<Vector_>& evecs) const {
    if (!_use_hevec) {
        throw std::runtime_error("Hessian eigenvectors computation deactivated");
    }
    evecs.resize(N);
    bool asw;
    if (N==2) asw = __gageProbe(ctx, u[0], u[1], 0);
    else asw = __gageProbe(ctx, u[0], u[1], u[2]);
    
    if (asw) return false;
    for (unsigned int i = 0 ; i < N ; i++) {
        evecs[0][i] = _hevec[  i];
        evecs[1][i] = _hevec[3+i];
        if (N==3) evecs[2][i] = _hevec[6+i];
    }
    return true;
}

#endif












