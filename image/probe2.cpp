#include "probe2.hpp"
#include <stdexcept>

using namespace spurt;

/*
inline bool gage_interface::scalar_wrapper::
__gageProbe(gageContext* ctx, double u, double v, double w) const
{
    if (!_use_wc) return gageProbe(ctx, u, v, w);
    vec3 idx = wcidx.to_index(vec3(u, v, w));
    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}

inline bool gage_interface::vector_wrapper::
__gageProbe(gageContext* ctx, double u, double v, double w) const
{
    if (!_use_wc) return gageProbe(ctx, u, v, w);
    vec3 idx = wcidx.to_index(vec3(u, v, w));
    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}

inline bool gage_interface::tensor_wrapper::
__gageProbe(gageContext* ctx, double u, double v, double w) const
{
    if (!_use_wc) return gageProbe(ctx, u, v, w);
    vec3 idx = wcidx.to_index(vec3(u, v, w));
    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}
*/

int gage_interface::set_kernel(gageContext* ctx, kernel_idx_t kidx) {
    if (kidx<0 || kidx>C4_HEXIC) {
        throw std::runtime_error("Unrecognized kernel type" + std::to_string(kidx));
    }

    double kparm[NRRD_KERNEL_PARMS_NUM];

    if (kidx == TENT) {
        kparm[0] = 1; // scale
        kparm[1] = kparm[2] = 0;
    }
    else if (kidx == BC_INTERP) {
        // parameter setting for Cubic BC kernel
        kparm[0] = 1.0; /* scale parameter, in units of samples */
        kparm[1] = 0.0; /* B */
        kparm[2] = 0.5; /* C */
    }
    else if (kidx == BC_BLUR) {
        kparm[0] = 1.0; // scale
        kparm[1] = 1.0; // B
        kparm[2] = 0.;  // C
    }
    else if (kidx == C3_QUINTIC || kidx == C4_HEXIC) {
        // parameter setting for Quintic C3 / C4 Hexic kernel
        kparm[0] = 1.0;
        kparm[1] = kparm[2] = 0;
    }
    else if (kidx == HANN || kidx == BLACK) {
        kparm[0] = 1; // scale
        kparm[1] = 3; // cut-off
        kparm[2] = 0;
    }
    else if (kidx == GAUSS) {
        kparm[0] = 2; // sigma
        kparm[1] = 3; // cut-off
        kparm[2] = 0;
    }

    int E = 0;
    if (kidx == TENT) {
        // linear
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelTent, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelCentDiff, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBCCubicDD, kparm);
    }
    else if (kidx == BC_BLUR || kidx == BC_INTERP) {
        // Cubic BC
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBCCubic, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBCCubicD, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBCCubicDD, kparm);
    }
    else if (kidx == BSPL3) {
        //
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBSpline3, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBSpline3D, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBSpline3DD, kparm);
    }
    else if (kidx == BSPL3_INTERP) {
        //
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBSpline3ApproxInverse, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBSpline3D, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBSpline3DD, kparm);
    }
    else if (kidx == BSPL5) {
        //
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBSpline5, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBSpline5D, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBSpline5DD, kparm);
    }
    else if (kidx == BSPL5_INTERP) {
        //
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBSpline5ApproxInverse, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBSpline5D, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBSpline5DD, kparm);
    }
    else if (kidx == BSPL7) {
        //
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBSpline7, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBSpline7D, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBSpline7DD, kparm);
    }
    else if (kidx == BSPL7_INTERP) {
        //
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBSpline7ApproxInverse, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBSpline7D, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBSpline7DD, kparm);
    }
    else if (kidx == C3_QUINTIC) {
        // Quintic C3
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelC3Quintic, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelC3QuinticD, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelC3QuinticDD, kparm);
    }
    else if (kidx == C4_HEXIC) {
        // Hexic C4 kernel
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelC4Hexic, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelC4HexicD, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelC4HexicDD, kparm);
    }
    else if (kidx == HANN) {
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelHann, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelHannD, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelHannDD, kparm);
    }
    else if (kidx == BLACK) {
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBlackman, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBlackmanD, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBlackmanDD, kparm);
    }
    else if (kidx == GAUSS) {
        if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelGaussian, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelGaussianD, kparm);
        if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelGaussianDD, kparm);
    }

    return E;
}

gage_interface::
scalar_wrapper::scalar_wrapper(const Nrrd* nrrd, kernel_idx_t kidx,
                               bool use_grad, bool use_hess,
                               bool use_heval, bool use_hevec)
        : wcidx(nrrd, nrrd->dim), _use_grad(use_grad),
          _use_hess(use_hess), _use_heval(use_heval),
          _use_hevec(use_hevec), _use_wc(false)
{
    if (gageKindVolumeCheck(gageKindScl, nrrd)) {
        std::cout << "ERROR: gage_interface:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }

    ctx = gageContextNew();

    int E = 0;
    if (!E) E |= !(pv = gagePerVolumeNew(ctx, nrrd, gageKindScl));
    if (!E) E |= gagePerVolumeAttach(ctx, pv);
    if (!E) E |= set_kernel(ctx, kidx);
    if (!E) E |= gageQueryItemOn(ctx, pv, gageSclValue);
    if (!E && _use_grad) E |= gageQueryItemOn(ctx, pv, gageSclGradVec);
    if (!E && _use_hess) E |= gageQueryItemOn(ctx, pv, gageSclHessian);
    if (!E && _use_heval) E |= gageQueryItemOn(ctx, pv, gageSclHessEval);
    if (!E && _use_hevec) E |= gageQueryItemOn(ctx, pv, gageSclHessEvec);
    if (!E) E |= gageUpdate(ctx);
    if (E) {
        std::cout << "ERROR: gage_interface:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }
    _v = gageAnswerPointer(ctx, pv, gageSclValue);
    if (_use_grad) _g = gageAnswerPointer(ctx, pv, gageSclGradVec);
    if (_use_hess) _h = gageAnswerPointer(ctx, pv, gageSclHessian);
    if (_use_heval) _heval = gageAnswerPointer(ctx, pv, gageSclHessEval);
    if (_use_hevec) _hevec = gageAnswerPointer(ctx, pv, gageSclHessEvec);
}

gage_interface::
scalar_wrapper::~scalar_wrapper()
{
    ctx = gageContextNix(ctx);
    pv = NULL;
}

bool gage_interface::
scalar_wrapper::value(const vec2& u, double& v) const
{
    if (gage_probe(ctx, u[0], u[1], 0, wcidx, _use_wc)) {
        v = 0;
        return false;
    }
    v = *_v;
    return true;
}

bool gage_interface::
scalar_wrapper::value(const vec3& u, double& v) const
{
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        v = 0;
        return false;
    }
    v = *_v;
    return true;
}

bool gage_interface::
scalar_wrapper::gradient(const vec2& u, vec2& g) const
{
    if (!_use_grad) {
        throw std::runtime_error("Gradient computation deactivated");
    }

    if (gage_probe(ctx, u[0], u[1], 0, wcidx, _use_wc)) {
        g[0] = g[1] = 0;
        return false;
    }
    g[0] = _g[0];
    g[1] = _g[1];
    return true;
}

bool gage_interface::
scalar_wrapper::gradient(const vec3& u, vec3& g) const
{
    if (!_use_grad) {
        throw std::runtime_error("Gradient computation deactivated");
    }

    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        g[0] = g[1] = g[2] = 0;
        return false;
    }
    g[0] = _g[0];
    g[1] = _g[1];
    g[2] = _g[2];
    return true;
}

bool gage_interface::
scalar_wrapper::hessian(const vec2& u, mat2& h) const
{
    if (!_use_hess) {
        throw std::runtime_error("Hessian computation deactivated");
    }

    if (gage_probe(ctx, u[0], u[1], 0.5, wcidx, _use_wc)) {
        h = mat2::Zero();
        return false;
    }
    h(0,0) = _h[0];
    h(0,1) = h(0,1) = _h[1];
    h(1,1) = _h[2];
    return true;
}

bool gage_interface::
scalar_wrapper::hessian(const vec3& u, mat3& h) const
{
    if (!_use_hess) {
        throw std::runtime_error("Hessian computation deactivated");
    }

    h = mat3::Zero();
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        return false;
    }
    h(0,0) = _h[0];
    h(0,1) = h(1,0) = _h[1];
    h(0,2) = h(2,0) = _h[2];
    h(1,1) = _h[4];
    h(1,2) = _h[5];
    h(2,2) = _h[8];
    return true;
}

bool gage_interface::
scalar_wrapper::hess_evals(const vec3& u, vec3& vals) const
{
    if (!_use_heval) {
        throw std::runtime_error("Hessian eigenvalues computation deactivated");
    }

    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        vals[0] = vals[1] = vals[2] = 0;
        return false;
    }
    vals[0] = _heval[0];
    vals[1] = _heval[1];
    vals[2] = _heval[2];
    return true;
}

// -----------------------------------------------------------------------------

gage_interface::
vector_wrapper::vector_wrapper(const Nrrd* nrrd, kernel_idx_t kidx,
                               bool use_jac)
        : wcidx(nrrd, nrrd->dim-1), _use_jac(use_jac), m_kernel_idx(kidx), _use_wc(false)
{
    if (gageKindVolumeCheck(gageKindVec, nrrd)) {
        std::cout << "gage_interface:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }

    ctx = gageContextNew();

    int E = 0;
    if (!E) E |= !(pv = gagePerVolumeNew(ctx, nrrd, gageKindVec));
    if (!E) E |= gagePerVolumeAttach(ctx, pv);
    if (!E) E |= set_kernel(ctx, kidx);
    if (!E) E |= gageQueryItemOn(ctx, pv, gageVecVector);
    if (!E && _use_jac) E |= gageQueryItemOn(ctx, pv, gageVecJacobian);
    if (!E) E |= gageUpdate(ctx);
    if (E) {
        std::cout << "gage_interface:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }
    _v = gageAnswerPointer(ctx, pv, gageVecVector);
    if (_use_jac) _j = gageAnswerPointer(ctx, pv, gageVecJacobian);
}

gage_interface::
vector_wrapper::vector_wrapper(const vector_wrapper& other)
	: vector_wrapper(other.pv->nin, other.m_kernel_idx, other._use_jac) {
        _use_wc = other._use_wc;
    }

gage_interface::
vector_wrapper::~vector_wrapper()
{
    ctx = gageContextNix(ctx);
    pv = NULL;
}

bool gage_interface::
vector_wrapper::value(const vec2& u, vec2& v) const
{
    if (gage_probe(ctx, u[0], u[1], 0, wcidx, _use_wc)) {
        v = vec2::Zero();
        return false;
    }
    v[0] = _v[0];
    v[1] = _v[1];
    return true;
}

bool gage_interface::
vector_wrapper::value(const vec3& u, vec3& v) const
{
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        v = vec3::Zero();
        return false;
    }
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];
    return true;
}

bool gage_interface::
vector_wrapper::jacobian(const vec2 &u, mat2 &j) const
{
    if (!_use_jac)
    {
        throw std::runtime_error("Jacobian computation deactivated");
    }

    if (gage_probe(ctx, u[0], u[1], 0, wcidx, _use_wc))
    {
        j = mat2::Zero();
        return false;
    }
    for (int i = 0; i < 4; ++i)
        j.data()[i] = _j[i];

    return true;
}

bool gage_interface::
vector_wrapper::jacobian(const vec3& u, mat3& j) const
{
    if (!_use_jac) {
        throw std::runtime_error("Jacobian computation deactivated");
    }

    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        j=mat3::Zero();
        return false;
    }
    for (int i=0; i<9 ; ++i) j.data()[i]=_j[i];

    return true;
}

// -----------------------------------------------------------------------------

gage_interface::
diffusion_tensor_wrapper::
diffusion_tensor_wrapper(const Nrrd* nrrd, kernel_idx_t kidx,
                         bool use_evals, bool use_evecs,
                         bool use_aniso, bool use_fa,
                         bool use_mode, bool use_fa_ridge, 
                         bool use_mode_ridge)
        : wcidx(nrrd, 3), _use_wc(false), _use_evals(use_evals), 
         _use_evecs(use_evecs), _use_fa(use_fa), _use_aniso(use_aniso),
         _use_fa_ridge(use_fa_ridge), _use_mode_ridge(use_mode_ridge)
{
    if (gageKindVolumeCheck(tenGageKind, nrrd)) {
        std::cout << "gage_interface #0:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }

    ctx = gageContextNew();
    int E = 0;
    if (!E) E |= !(pv = gagePerVolumeNew(ctx, nrrd, tenGageKind));
    if (!E) E |= gagePerVolumeAttach(ctx, pv);
    if (!E) E |= set_kernel(ctx, kidx);

    // 1st order
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageTensor);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageConfidence);
    if (_use_aniso && !E)
        E |= gageQueryItemOn(ctx, pv, tenGageAniso);
    if (_use_evals && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageEval);
    if (_use_evecs && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageEvec);

    std::cout << "1st order OK" << std::endl;

    // FA Hessian
    if (_use_fa && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageFA);
    if (_use_fa_ridge && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageFAGradVec);
    if (_use_fa_ridge && !E)
        E |= gageQueryItemOn(ctx, pv, tenGageFAHessian);
    if (_use_fa_ridge && !E)
        E |= gageQueryItemOn(ctx, pv, tenGageFAHessianEval);
    if (_use_fa_ridge && !E)
        E |= gageQueryItemOn(ctx, pv, tenGageFAHessianEvec);
    if (_use_fa_ridge && !E)
        E |= gageQueryItemOn(ctx, pv, tenGageFARidgeSurfaceStrength);
    if (_use_fa_ridge && !E)
        E |= gageQueryItemOn(ctx, pv, tenGageFAValleySurfaceStrength);

    std::cout << "FA items OK" << std::endl;

    // Mode Hessian
    if (_use_mode_ridge && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageModeHessian);
    if (_use_mode_ridge && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageModeGradVec);
    if (_use_mode_ridge && !E) 
        E |= gageQueryItemOn(ctx, pv, tenGageModeHessianEval);
    if (_use_mode_ridge && !E) 
    E |= gageQueryItemOn(ctx, pv, tenGageModeHessianEvec);

    std::cout << "mode items OK" << std::endl;

    if (!E) E |= gageUpdate(ctx);

    if (E) {
        std::cout << "gage_interface #1:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }
    _t = gageAnswerPointer(ctx, pv, tenGageTensor);
    _c = gageAnswerPointer(ctx, pv, tenGageConfidence);
    if (_use_aniso)
        _aniso = gageAnswerPointer(ctx, pv, tenGageAniso);
    if (_use_evecs)
        _evec = gageAnswerPointer(ctx, pv, tenGageEvec);
    if (_use_evals)
        _eval = gageAnswerPointer(ctx, pv, tenGageEval);

    if (_use_fa)
        _fa = gageAnswerPointer(ctx, pv, tenGageFA);
    if (_use_fa_ridge)
    {
        _fa_grad = gageAnswerPointer(ctx, pv, tenGageFAGradVec);
        _fa_hess = gageAnswerPointer(ctx, pv, tenGageFAHessian);
        _fa_hess_eval = gageAnswerPointer(ctx, pv, tenGageFAHessianEval);
        _fa_hess_evec = gageAnswerPointer(ctx, pv, tenGageFAHessianEvec);
        _fa_ridge_surf_strength = gageAnswerPointer(ctx, pv, tenGageFARidgeSurfaceStrength);
        _fa_valley_surf_strength = gageAnswerPointer(ctx, pv, tenGageFAValleySurfaceStrength);
    }

    if (_use_mode)
        _mode = gageAnswerPointer(ctx, pv, tenGageMode);
    if (_use_mode_ridge)
    {
        _mode_grad = gageAnswerPointer(ctx, pv, tenGageModeGradVec);
        _mode_hess = gageAnswerPointer(ctx, pv, tenGageModeHessian);
        _mode_hess_eval = gageAnswerPointer(ctx, pv,      tenGageModeHessianEval);
        _mode_hess_evec = gageAnswerPointer(ctx, pv,      tenGageModeHessianEvec);
    }
}

// --------------------------------------------------------------------------
// fixed_vector based API

bool gage_interface::
diffusion_tensor_wrapper::value(const vec3 &u, mat3 &v) const
{
    v = mat3::Zero();
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc))
    {
        return false;
    }
    v(0,0) = _t[1];
    v(1,0) = v(0,1) = _t[2];
    v(2,0) = v(0,2) = _t[3];
    v(1,1) = _t[4];
    v(2,1) = v(1,2) = _t[5];
    v(2,2) = _t[6];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::fa(const vec3& u, double& v) const
{
    if (!_use_fa)
    {
        throw std::runtime_error("FA computation deactivated");
    }
    v = 0;
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        return false;
    }
    v = *_fa;
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::fa_grad(const vec3& u, vec3& g) const
{
    if (!_use_fa_ridge)
    {
        throw std::runtime_error("FA gradient computation deactivated");
    }
    g[0] = g[1] = g[2] = 0;
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        return false;
    }
    g[0] = _fa_grad[0];
    g[1] = _fa_grad[1];
    g[2] = _fa_grad[2];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::fa_hess(const vec3& u, mat3& h) const
{
    if (!_use_fa_ridge)
    {
        throw std::runtime_error("FA hessian computation deactivated");
    }
    h = mat3::Zero();
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        return false;
    }
    h(0,0) = _fa_hess[0];
    h(1,0) = h(0,1) = _fa_hess[1];
    h(2,0) = h(0,2) = _fa_hess[2];
    h(1,1) = _fa_hess[4];
    h(2,1) = h(1,2) = _fa_hess[5];
    h(2,2) = _fa_hess[8];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::fa_hess_evals(const vec3& u, vec3& vals) const
{
    if (!_use_fa_ridge)
    {
        throw std::runtime_error("FA hessian eigenvalues computation deactivated");
    }
    vals = vec3::Zero();
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        return false;
    }
    vals[0] = _fa_hess_eval[0];
    vals[1] = _fa_hess_eval[1];
    vals[2] = _fa_hess_eval[2];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::fa_hess_evecs(const vec3& u,
                                        mat3& vecs) const
{
    if (!_use_fa_ridge)
    {
        throw std::runtime_error("FA hessian eigenvectors computation deactivated");
    }
    vecs = mat3::Zero();
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        return false;
    }
    for (unsigned int row= 0 ; row < 3 ; row++) {
        vecs(row, 0) = _fa_hess_evec[  row];
        vecs(row, 1) = _fa_hess_evec[3+row];
        vecs(row, 2) = _fa_hess_evec[6+row];
    }
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::mode(const vec3& u, double& v) const
{
    if (!_use_mode)
    {
        throw std::runtime_error("mode computation deactivated");
    }
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        v = 0;
        return false;
    }
    v = *_mode;
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::mode_grad(const vec3& u, vec3& g) const
{
    if (!_use_mode_ridge)
    {
        throw std::runtime_error("mode gradient computation deactivated");
    }
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        g[0] = g[1] = g[2] = 0;
        return false;
    }
    g[0] = _mode_grad[0];
    g[1] = _mode_grad[1];
    g[2] = _mode_grad[2];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::mode_hess(const vec3& u, mat3& h) const
{
    if (!_use_mode_ridge)
    {
        throw std::runtime_error("mode hessian computation deactivated");
    }
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        h = mat3::Zero();
        return false;
    }
    h(0,0) = _mode_hess[0];
    h(1,0) = h(0,1) = _mode_hess[1];
    h(2,0) = h(0,2) = _mode_hess[2];
    h(1,1) = _mode_hess[4];
    h(2,1) = h(1,2) = _mode_hess[5];
    h(2,2) = _mode_hess[8];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::mode_hess_evals(const vec3& u, vec3& vals) const
{
    if (!_use_mode_ridge)
    {
        throw std::runtime_error("mode hessian eigenvalues computation deactivated");
    }
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        vals[0] = vals[1] = vals[2] = 0;
        return false;
    }
    vals[0] = _mode_hess_eval[0];
    vals[1] = _mode_hess_eval[1];
    vals[2] = _mode_hess_eval[2];
    return true;
}

bool gage_interface::
diffusion_tensor_wrapper::mode_hess_evecs(const vec3& u, mat3& vecs) const
{
    if (!_use_mode_ridge)
    {
        throw std::runtime_error("mode hessian eigenvectors computation deactivated");
    }
    if (gage_probe(ctx, u[0], u[1], u[2], wcidx, _use_wc)) {
        vecs = mat3::Zero();
        return false;
    }
    for (unsigned int row = 0; row < 3; row++)
    {
        vecs(row, 0) = _mode_hess_evec[row];
        vecs(row, 1) = _mode_hess_evec[3 + row];
        vecs(row, 2) = _mode_hess_evec[6 + row];
    }
    return true;
}