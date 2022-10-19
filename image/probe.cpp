#include "probe.hpp"
#include <stdexcept>

using namespace xavier;

// inline bool gage_interface::scalar_wrapper::
// __gageProbe(gageContext* ctx, double u, double v, double w) const
// {
//     if (!_use_wc) return gageProbe(ctx, u, v, w);
//     nvis::vec3 idx = wcidx.to_index(nvis::vec3(u, v, w));
//     return gageProbe(ctx, idx[0], idx[1], idx[2]);
// }

inline bool gage_interface::vector_wrapper::
__gageProbe(gageContext* ctx, double u, double v, double w) const
{
    if (!_use_wc) return gageProbe(ctx, u, v, w);
    nvis::vec3 idx = wcidx.to_index(nvis::vec3(u, v, w));

    // std::cout << "world position (" << u << ", " << v << ", " << w << ") "
    //     << "mapped to " << idx << "\n";

    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}

inline bool gage_interface::tensor_wrapper::
__gageProbe(gageContext* ctx, double u, double v, double w) const
{
    if (!_use_wc) return gageProbe(ctx, u, v, w);
    nvis::vec3 idx = wcidx.to_index(nvis::vec3(u, v, w));
    return gageProbe(ctx, idx[0], idx[1], idx[2]);
}

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
scalar_wrapper::value(double* u, double& v) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        v = 0;
        return false;
    }
    v = *_v;
    return true;
}

bool gage_interface::
scalar_wrapper::value(const nvis::vec1& u, double& v) const
{
    if (__gageProbe(ctx, u[0], 0, 0)) {
        v = 0;
        return false;
    }
    v = *_v;
    return true;
}

bool gage_interface::
scalar_wrapper::value(const nvis::vec2& u, double& v) const
{
    if (__gageProbe(ctx, u[0], u[1], 0)) {
        v = 0;
        return false;
    }
    v = *_v;
    return true;
}

bool gage_interface::
scalar_wrapper::value(const nvis::vec3& u, double& v) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        v = 0;
        return false;
    }
    v = *_v;
    return true;
}

bool gage_interface::
scalar_wrapper::gradient(double* u, double* g) const
{
    if (!_use_grad) {
        throw std::runtime_error("Gradient computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        g[0] = g[1] = 0;
        return false;
    }
    g[0] = _g[0];
    g[1] = _g[1];
    return true;
}

bool gage_interface::
scalar_wrapper::gradient(const nvis::vec2& u, nvis::vec2& g) const
{
    if (!_use_grad) {
        throw std::runtime_error("Gradient computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], 0)) {
        g[0] = g[1] = 0;
        return false;
    }
    g[0] = _g[0];
    g[1] = _g[1];
    return true;
}

bool gage_interface::
scalar_wrapper::gradient(const nvis::vec3& u, nvis::vec3& g) const
{
    if (!_use_grad) {
        throw std::runtime_error("Gradient computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        g[0] = g[1] = g[2] = 0;
        return false;
    }
    g[0] = _g[0];
    g[1] = _g[1];
    g[2] = _g[2];
    return true;
}

bool gage_interface::
scalar_wrapper::hessian(double* u, double* h) const
{
    if (!_use_hess) {
        throw std::runtime_error("Hessian computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        h[0] = h[1] = h[2] = 0;
        return false;
    }
    h[0] = _h[0];
    h[1] = _h[1];
    h[2] = _h[4];
    return true;
}

bool gage_interface::
scalar_wrapper::hessian(const nvis::vec2& u, nvis::vec3& h) const
{
    if (!_use_hess) {
        throw std::runtime_error("Hessian computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], 0.5)) {
        h[0] = h[1] = h[2] = 0;
        return false;
    }
    h[0] = _h[0];
    h[1] = _h[1];
    h[2] = _h[4];
    return true;
}

bool gage_interface::
scalar_wrapper::hessian(const nvis::vec3& u, nvis::vec6& h) const
{
    if (!_use_hess) {
        throw std::runtime_error("Hessian computation deactivated");
    }

    h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    h[0] = _h[0];
    h[1] = _h[1];
    h[2] = _h[2];
    h[3] = _h[4];
    h[4] = _h[5];
    h[5] = _h[8];
    return true;
}

bool gage_interface::
scalar_wrapper::hess_evals(const nvis::vec3& u, nvis::vec3& vals) const
{
    if (!_use_heval) {
        throw std::runtime_error("Hessian eigenvalues computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        vals[0] = vals[1] = vals[2] = 0;
        return false;
    }
    vals[0] = _heval[0];
    vals[1] = _heval[1];
    vals[2] = _heval[2];
    return true;
}

bool gage_interface::
scalar_wrapper::hess_evecs(const nvis::vec3& u,
                           std::vector< nvis::vec3 >& vecs) const
{
    if (!_use_hevec) {
        throw std::runtime_error("Hessian eigenvectors computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        vecs[0][i] = _hevec[  i];
        vecs[1][i] = _hevec[3+i];
        vecs[2][i] = _hevec[6+i];
    }
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
vector_wrapper::value(const nvis::vec2& u, nvis::vec2& v) const
{
    if (__gageProbe(ctx, u[0], u[1], 0)) {
        v = 0;
        return false;
    }
    v[0] = _v[0];
    v[1] = _v[1];
    return true;
}

bool gage_interface::
vector_wrapper::value(const nvis::vec3& u, nvis::vec3& v) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        v = 0;
        return false;
    }
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];
    return true;
}

bool gage_interface::
vector_wrapper::jacobian(const nvis::vec2& u, deriv2_t& j) const
{
    if (!_use_jac) {
        throw std::runtime_error("Jacobian computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], 0)) {
        std::cout << "gage return an error\n";
        j=deriv2_t::Zero();
        return false;
    }
    for (int i=0; i<4 ; ++i) j.data()[i]=_j[i];
    return true;
}

bool gage_interface::
vector_wrapper::jacobian(const nvis::vec3& u, deriv3_t& j) const
{
    if (!_use_jac) {
        throw std::runtime_error("Jacobian computation deactivated");
    }

    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        j=deriv3_t::Zero();
        return false;
    }
    for (int i=0; i<9 ; ++i) j.data()[i]=_j[i];

    return true;
}

// -----------------------------------------------------------------------------

gage_interface::
tensor_wrapper::tensor_wrapper(const Nrrd* nrrd, kernel_idx_t kidx)
        : wcidx(nrrd, 3), _use_wc(false)
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
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageAniso);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageEval);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageEvec);

    std::cout << "1st order OK" << std::endl;

    // FA Hessian
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFA);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFAGradVec);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFAHessian);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFAHessianEval);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFAHessianEvec);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFARidgeSurfaceStrength);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageFAValleySurfaceStrength);

    std::cout << "FA items OK" << std::endl;

    // Mode Hessian
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageModeHessian);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageModeGradVec);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageModeHessianEval);
    if (!E) E |= gageQueryItemOn(ctx, pv, tenGageModeHessianEvec);

    std::cout << "mode items OK" << std::endl;

    if (!E) E |= gageUpdate(ctx);

    if (E) {
        std::cout << "gage_interface #1:" << std::endl
        << biffGetDone(GAGE)
        << std::endl;
        return;
    }
    _t =                        gageAnswerPointer(ctx, pv, tenGageTensor);
    _c =                        gageAnswerPointer(ctx, pv, tenGageConfidence);
    _aniso =                    gageAnswerPointer(ctx, pv, tenGageAniso);
    _evec =                     gageAnswerPointer(ctx, pv, tenGageEvec);
    _eval =                     gageAnswerPointer(ctx, pv, tenGageEval);

    _fa =                       gageAnswerPointer(ctx, pv, tenGageFA);
    _fa_grad =                  gageAnswerPointer(ctx, pv, tenGageFAGradVec);
    _fa_hess =                  gageAnswerPointer(ctx, pv, tenGageFAHessian);
    _fa_hess_eval =             gageAnswerPointer(ctx, pv, tenGageFAHessianEval);
    _fa_hess_evec =             gageAnswerPointer(ctx, pv, tenGageFAHessianEvec);
    _fa_ridge_surf_strength =   gageAnswerPointer(ctx, pv, tenGageFARidgeSurfaceStrength);
    _fa_valley_surf_strength =  gageAnswerPointer(ctx, pv, tenGageFAValleySurfaceStrength);

    _mode =                     gageAnswerPointer(ctx, pv, tenGageMode);
    _mode_grad =                gageAnswerPointer(ctx, pv, tenGageModeGradVec);
    _mode_hess =                gageAnswerPointer(ctx, pv, tenGageModeHessian);
    _mode_hess_eval =           gageAnswerPointer(ctx, pv, tenGageModeHessianEval);
    _mode_hess_evec =           gageAnswerPointer(ctx, pv, tenGageModeHessianEvec);
}

gage_interface::
tensor_wrapper::~tensor_wrapper()
{
    ctx = gageContextNix(ctx);
    pv = NULL;
}

bool gage_interface::
tensor_wrapper::tensor(double* u, double* t) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 7 ; i++) {
            t[i] = 0;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < 7 ; i++) {
        t[i] = _t[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::confidence(double* u, double& c) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        c = 0;
        return false;
    }
    c = *_c;
    return true;
}

bool gage_interface::
tensor_wrapper::confidence(const nvis::vec3& u, double& val) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        val = 0;
        return false;
    }
    val = *_c;
    return true;
}

bool gage_interface::
tensor_wrapper::aniso(double* u, int anisokind, double& aniso) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        aniso = 0;
        return false;
    }
    aniso = _aniso[anisokind];
    return true;
}

bool gage_interface::
tensor_wrapper::evals(double* u, double* eval) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        eval[0] = eval[1] = eval[2] = 0;
        return false;
    }
    eval[0] = _eval[0];
    eval[1] = _eval[1];
    eval[2] = _eval[2];
    return true;
}

bool gage_interface::
tensor_wrapper::evecs(double* u, double* evec) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 9 ; i++) {
            evec[i] = 0;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < 9 ; i++) {
        evec[i] = _evec[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::fa(double* u, double* fa) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        *fa = 0;
        return false;
    }
    *fa = *_fa;
    return true;
}

bool gage_interface::
tensor_wrapper::fa_grad(double* u, double* grad) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 3 ; i++) {
            grad[0] = 0;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        grad[i] = _fa_grad[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::fa_hess(double* u, double* hess) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 6 ; i++) {
            hess[i] = 0;
        }
        return false;
    }
    hess[0] = _fa_hess[0];
    hess[1] = _fa_hess[1];
    hess[2] = _fa_hess[2];
    hess[3] = _fa_hess[4];
    hess[4] = _fa_hess[5];
    hess[5] = _fa_hess[8];
    return true;
}

bool gage_interface::
tensor_wrapper::fa_hess_evals(double* u, double* eval) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 3 ; i++) {
            eval[i] = 0;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        eval[i] = _fa_hess_eval[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::fa_hess_evecs(double* u, double* evec) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 9 ; i++) {
            evec[i] = 0;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < 9 ; i++) {
        evec[i] = _fa_hess_evec[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::fa_ridge_surf_strength(double* u, double* val) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        *val = 0;
        return false;
    }
    *val = *_fa_ridge_surf_strength;
    return true;
}

bool gage_interface::
tensor_wrapper::fa_valley_surf_strength(double* u, double* val) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        *val = 0;
        return false;
    }
    *val = *_fa_valley_surf_strength;
    return true;
}

bool gage_interface::
tensor_wrapper::mode(double* u, double* mode) const
{
    *mode = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    *mode = *_mode;
    return true;
}

bool gage_interface::
tensor_wrapper::mode_grad(double* u, double* grad) const
{
    grad[0] = grad[1] = grad[2] = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        grad[i] = _mode_grad[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::mode_hess(double* u, double* hess) const
{
    for (unsigned int i = 0 ; i < 9 ; i++) hess[i] = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    for (unsigned int i = 0 ; i < 9 ; i++) {
        hess[i] = _mode_hess[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::mode_hess_evals(double* u, double* eval) const
{
    eval[0] = eval[1] = eval[2];
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        eval[i] = _mode_hess_eval[i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::mode_hess_evecs(double* u, double* evec) const
{
    for (unsigned int i = 0 ; i < 9 ; i++) {
        evec[i] = 0;
    }
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    for (unsigned int i = 0 ; i < 9 ; i++) {
        evec[i] = _mode_hess_evec[i];
    }
    return true;
}

// --------------------------------------------------------------------------
// nvis::fixed_vector based API

bool gage_interface::
tensor_wrapper::fa(const nvis::vec3& u, double& v) const
{
    v = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    v = *_fa;
    return true;
}

bool gage_interface::
tensor_wrapper::fa_grad(const nvis::vec3& u, nvis::vec3& g) const
{
    g[0] = g[1] = g[2] = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    g[0] = _fa_grad[0];
    g[1] = _fa_grad[1];
    g[2] = _fa_grad[2];
    return true;
}

bool gage_interface::
tensor_wrapper::fa_hess(const nvis::vec3& u, nvis::vec6& h) const
{
    h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    h[0] = _fa_hess[0];
    h[1] = _fa_hess[1];
    h[2] = _fa_hess[2];
    h[3] = _fa_hess[4];
    h[4] = _fa_hess[5];
    h[5] = _fa_hess[8];
    return true;
}

bool gage_interface::
tensor_wrapper::fa_hess_evals(const nvis::vec3& u, nvis::vec3& vals) const
{
    vals[0] = vals[1] = vals[2] = 0;
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    vals[0] = _fa_hess_eval[0];
    vals[1] = _fa_hess_eval[1];
    vals[2] = _fa_hess_eval[2];
    return true;
}

bool gage_interface::
tensor_wrapper::fa_hess_evecs(const nvis::vec3& u,
                              std::vector< nvis::vec3 >& vecs) const
{
    for (unsigned int i = 0 ; i < 3 ; i++) {
        vecs[i][0] = vecs[i][1] = vecs[i][2] = 0;
    }
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        vecs[0][i] = _fa_hess_evec[  i];
        vecs[1][i] = _fa_hess_evec[3+i];
        vecs[2][i] = _fa_hess_evec[6+i];
    }
    return true;
}

bool gage_interface::
tensor_wrapper::mode(const nvis::vec3& u, double& v) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        v = 0;
        return false;
    }
    v = *_mode;
    return true;
}

bool gage_interface::
tensor_wrapper::mode_grad(const nvis::vec3& u, nvis::vec3& g) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        g[0] = g[1] = g[2] = 0;
        return false;
    }
    g[0] = _mode_grad[0];
    g[1] = _mode_grad[1];
    g[2] = _mode_grad[2];
    return true;
}

bool gage_interface::
tensor_wrapper::mode_hess(const nvis::vec3& u, nvis::vec6& h) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        h[0] = h[1] = h[2] = h[3] = h[4] = h[5] = 0;
        return false;
    }
    h[0] = _mode_hess[0];
    h[1] = _mode_hess[1];
    h[2] = _mode_hess[2];
    h[3] = _mode_hess[4];
    h[4] = _mode_hess[5];
    h[5] = _mode_hess[8];
    return true;
}

bool gage_interface::
tensor_wrapper::mode_hess_evals(const nvis::vec3& u, nvis::vec3& vals) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        vals[0] = vals[1] = vals[2] = 0;
        return false;
    }
    vals[0] = _mode_hess_eval[0];
    vals[1] = _mode_hess_eval[1];
    vals[2] = _mode_hess_eval[2];
    return true;
}

bool gage_interface::
tensor_wrapper::mode_hess_evecs(const nvis::vec3& u,
                                std::vector< nvis::vec3 >& vecs) const
{
    if (__gageProbe(ctx, u[0], u[1], u[2])) {
        for (unsigned int i = 0 ; i < 3 ; i++) {
            vecs[i][0] = vecs[i][1] = vecs[i][2] = 0;
        }
        return false;
    }
    for (unsigned int i = 0 ; i < 3 ; i++) {
        vecs[0][i] = _mode_hess_evec[  i];
        vecs[1][i] = _mode_hess_evec[3+i];
        vecs[2][i] = _mode_hess_evec[6+i];
    }
    return true;
}
