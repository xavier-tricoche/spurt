#include <iostream>

#include <image/probe.hpp>
#include <crease/measure_wrapper.hpp>

using namespace spurt;
using namespace std;

bool wrong;



double spurt::MeasureWrapper::
my_directional_derivative(const vec3& p, const unsigned int dim, double eps) const
{
    vec3 q(p);
    q[dim] += eps;
    double vplus = this->value(q);
    q[dim] -= 2 * eps;
    double vminus = this->value(q);
    return (vplus -vminus) / (2*eps);
}

spurt::vec3 spurt::MeasureWrapper::my_gradient(const vec3& p, double eps) const
{
    vec3 g;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        g[i] = my_directional_derivative(p, i, eps);
    }
    return g;
}

spurt::vec6 spurt::MeasureWrapper::my_hessian(const vec3& p, double eps) const
{
    vec6 h;
    vec3 g[6];
    vec3 q(p);
    for (unsigned int dim = 0 ; dim < 3 ; ++dim) {
        q[dim] += eps;
        g[2*dim] = my_gradient(q, eps);
        q[dim] -= 2 * eps;
        g[2*dim+1] = my_gradient(q, eps);
        q[dim] += eps;
    }
    unsigned int n = 0;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        for (unsigned int j = i ; j < 3 ; ++j) {
            h[n++] = (g[2*j][i] - g[2*j+1][i]) / (2 * eps);
        }
    }
    return h;
}

spurt::vec3 spurt::MeasureWrapper::
my_eigenvector(const vec3& p, unsigned int idx, double eps) const
{
    vec6 h = my_hessian(p, eps);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    return vec3(evec[3*idx], evec[3*idx+1], evec[3*idx+2]);
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

int spurt::MeasureWrapper::what_measure() const
{
    return _measure;
}

// --------------------------------------------------------------------------

void spurt::MeasureWrapper::turn_buffer_on() const
{
    active_buffer = true;
}

void spurt::MeasureWrapper::turn_buffer_off() const
{
    buffer.clear();
    active_buffer = false;
}

bool spurt::MeasureWrapper::is_stored(measured_values& values, const vec3& p) const
{
    std::map< vec3, measured_values, eps_order >::const_iterator it =
        buffer.find(p);
    if (it != buffer.end()) {
        values = it->second;
        return true;
    }

    return false;
}

void spurt::MeasureWrapper::store(measured_values& mv, const vec3& p) const
{
    mv.conf = this->_confidence(p);
    mv.val = this->_value(p);
    mv.g = this->_gradient(p);
    mv.Hg = to_mat3(this->hessian(p))*mv.g;
    for (unsigned int i = 0 ; i < 3 ; ++i) {
        mv.e[i] = this->_eigenvector(p, i);
        mv.eval[i] = this->_eigenvalue(p, i);
    }

    buffer[p] = mv;
}

// --------------------------------------------------------------------------

double spurt::MeasureWrapper::confidence(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);
        return mv.conf;
    }

    return this->_confidence(p);
}

// --------------------------------------------------------------------------

double spurt::MeasureWrapper::_confidence(const vec3& p) const
{
    double v;
    switch (_measure) {
    case 0: {
        v = 1;
        break;
    }
    case 1: {
        _tens_wrap->confidence(p, v);
        break;
    }
    case 2: {
        _tens_wrap->confidence(p, v);
        break;
    }
    default: {
        assert(false);
    }
    }

    return v;
}

// --------------------------------------------------------------------------

double spurt::MeasureWrapper::value(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);
        return mv.val;
    }

    return this->_value(p);
}

// --------------------------------------------------------------------------

double spurt::MeasureWrapper::_value(const vec3& p) const
{
    double v;

    switch (_measure) {
    case 0: {
        _scal_wrap->value(p, v);
        break;
    }
    case 1: {
        _tens_wrap->fa(p, v);
        break;
    }
    case 2: {
        _tens_wrap->mode(p, v);
        break;
    }
    default: {
        assert(false);
    }
    }

    return v;
}

// --------------------------------------------------------------------------

spurt::vec3 spurt::MeasureWrapper::eigenvector(const vec3& p, unsigned int idx) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.e[idx];
    }

    return this->_eigenvector(p, idx);
}

// --------------------------------------------------------------------------

spurt::vec3
spurt::MeasureWrapper::_eigenvector(const vec3& p, unsigned int idx) const
{
#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->hess_evecs(p, _evecs);
        break;
    }
    case 1: {
        _tens_wrap->fa_hess_evecs(p, _evecs);
        break;
    }
    case 2: {
        _tens_wrap->mode_hess_evecs(p, _evecs);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec6 h = my_hessian(p);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    vec3 ev(evec[3*idx], evec[3*idx+1], evec[3*idx+2]);
    // std::cout << "measured eigenvector = " << _evecs[idx] << ", approx eigenvector = " << ev << std::endl;

    return ev;
#endif

    return _evecs[idx];
}

// --------------------------------------------------------------------------

double spurt::MeasureWrapper::eigenvalue(const vec3& p, unsigned int idx) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.eval[idx];
    }

    return this->_eigenvalue(p, idx);
}

// --------------------------------------------------------------------------

double
spurt::MeasureWrapper::_eigenvalue(const vec3& p, unsigned int idx) const
{
    vec3 _tmp;

#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->hess_evals(p, _tmp);
        break;
    }
    case 1: {
        _tens_wrap->fa_hess_evals(p, _tmp);
        break;
    }
    case 2: {
        _tens_wrap->mode_hess_evals(p, _tmp);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec6 h = my_hessian(p);
    double mat[9] = { h[0], h[1], h[2],
                      h[1], h[3], h[4],
                      h[2], h[4], h[5]
                    };
    double eval[3], evec[9];
    ell_3m_eigensolve_d(eval, evec, mat, 1);
    // std::cout << "measured eigenvalue = " << _tmp[idx] << ", approx eigenvalue = " << eval[idx] << std::endl;
    return eval[idx];
#endif

    return _tmp[idx];
}

// --------------------------------------------------------------------------

spurt::vec3 spurt::MeasureWrapper::gradient(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.g;
    }

    return this->_gradient(p);
}

// --------------------------------------------------------------------------

spurt::vec3 spurt::MeasureWrapper::Hgradient(const vec3& p) const
{
    if (active_buffer) {
        measured_values mv;
        if (!is_stored(mv, p))
            store(mv, p);

        return mv.Hg;
    }

    return to_mat3(this->hessian(p)) * this->_gradient(p);
}

// --------------------------------------------------------------------------

spurt::vec3 spurt::MeasureWrapper::_gradient(const vec3& p) const
{
#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->gradient(p, _grad);
        break;
    }
    case 1: {
        _tens_wrap->fa_grad(p, _grad);
        break;
    }
    case 2: {
        _tens_wrap->mode_grad(p, _grad);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec3 g = my_gradient(p);
    // std::cout << "measured gradient = " << _grad << ", approx gradient = " << g << std::endl;
    return g;
#endif

    return _grad;
}

// --------------------------------------------------------------------------

spurt::vec6 spurt::MeasureWrapper::hessian(const vec3& p) const
{
    vec6 _tmp;

#ifndef DEBUG
    switch (_measure) {
    case 0: {
        _scal_wrap->hessian(p, _tmp);
        break;
    }
    case 1: {
        _tens_wrap->fa_hess(p, _tmp);
        break;
    }
    case 2: {
        _tens_wrap->mode_hess(p, _tmp);
        break;
    }
    default: {
        assert(false);
    }
    }

#else
    vec6 h = my_hessian(p);
    // std::cout << "measured hessian = " << h << ", approx hessian = " << h << std::endl;
    return h;
#endif

    return _tmp;
}

// ------------------------------------------------------------------------

MeasureWrapper::MeasureWrapper(const Nrrd* nrrd, int aniso, kernel_idx_t kidx)
        : _scal_wrap(0), _tens_wrap(0), _evecs(3), _evals(3), _measure(aniso)
{
    // check if this is a tensor field
    bool is_tensor = (nrrd->dim == 4 && nrrd->axis[0].size == 7);

    active_buffer = false;

    std::cout << "measure is " << _measure << std::endl;

    if (is_tensor) {
        _tens_wrap = new gage_interface::tensor_wrapper(nrrd, kidx);
        _tens_wrap->use_world();
    }
    else {
        _scal_wrap = new gage_interface::scalar_wrapper(nrrd, kidx);
        _scal_wrap->use_world();
    }
}

// --------------------------------------------------------------------------

unsigned int N;
vec3 _gradient(const MeasureWrapper& wrap, const vec3& p)
{
    double v[2], dv_prev, dv_cur;
    vec3 q, grad;

    wrong = false;

    for (unsigned int i = 0 ; i < 3 ; i++) {
        double h = 0.01;
        N = 0;
        for (unsigned int n = 0 ; true ; n++, N++) {
            q = p;
            q[i] -= h;
            v[0] = wrap.value(q);
            q[i] += 2 * h;
            v[1] = wrap.value(q);

            dv_cur = 0.5 * (v[1] - v[0]) / h;
            std::cout << dv_cur << "... " << std::flush;
            if (n) {
                double delta = (dv_cur - dv_prev);
                double denom = (dv_prev != 0. ? dv_prev : 1.0);
                double err = fabs(delta / denom);
                if (err < 1.0e-6) {
                    grad[i] = dv_cur;
                    break;
                }
            }

            dv_prev = dv_cur;
            h *= 0.5;
            if (h < 1.0e-5) {
                std::cout << "unable to compute derivative" << std::endl;
                wrong = true;
                return grad; // meaningless value
            }
        }
        std::cout << std::endl;
    }

    return grad;
}

// --------------------------------------------------------------------------

vec6 _hessian(const MeasureWrapper& wrap, const vec3& p)
{
    vec3 grad[2], dg_cur, dg_prev;
    vec3 q;
    vec6 hess;

    wrong = false;

    for (unsigned int i = 0 ; i < 3 ; i++) {
        double h = 0.01;
        N = 0;
        for (unsigned int n = 0 ; true ; n++, N++) {
            q = p;
            q[i] -= h;
            grad[0] = wrap.gradient(q);
            q[i] += 2 * h;
            grad[1] = wrap.gradient(q);

            dg_cur = 0.5 / h * (grad[1] - grad[0]);
            std::cout << dg_cur << "... " << std::flush;
            if (n) {
                double delta = spurt::norm(vec3(dg_cur - dg_prev));
                double denom = spurt::norm(dg_prev);
                double err = fabs(delta / denom);
                if (err < 1.0e-6) {
                    switch (i) {
                    case 0:
                        hess[0] = dg_cur[0];
                        hess[1] = dg_cur[1];
                        hess[2] = dg_cur[2];
                        break;
                    case 1:
                        hess[3] = dg_cur[1];
                        hess[4] = dg_cur[2];
                        break;
                    case 2:
                        hess[5] = dg_cur[2];
                        break;
                    default:
                        break;
                    }
                    break;
                }
            }

            dg_prev = dg_cur;
            h *= 0.5;
            if (h < 1.0e-5) {
                std::cout << "unable to compute derivative" << std::endl;
                wrong = true;
                return hess; // meaningless value
            }
        }
        std::cout << std::endl;
    }

    return hess;
}

// --------------------------------------------------------------------------

void spurt::test_derivatives(const Nrrd* nrrd, int aniso)
{
    srand48(time(0));
    MeasureWrapper wrapper(nrrd, aniso);

    bbox3 box = spurt::nrrd_utils::get_bounds<3>(nrrd);

    // check gradient
    double avg = 0, total_max = 0;
    unsigned int k = 0;
    for (unsigned int n = 0 ; n < 100 ; k++) {
        std::cout << "computing gradient at point #" << k
                  << " (" << n << ")" << std::endl;
        vec3 ran(drand48(), drand48(), drand48());
        vec3 p = box.min() + ran * box.size();

        // make sure we are at a valid position within the volume
        if (wrapper.confidence(p) < 0.5) continue;

        // compute gradient
        vec3 grad_ref = wrapper.gradient(p);
        vec3 grad_app = _gradient(wrapper, p);
        if (wrong) continue;

        double scale = grad_ref[0] / grad_app[0];
        grad_app *= scale;

        std::cout << "gage: " << grad_ref << std::endl;
        std::cout << "approx: " << grad_app << " (" << N << ")" << std::endl;
        std::cout << "scale: " << scale << std::endl;

        double max = 0;
        for (unsigned int i = 0 ; i < 3 ; i++) {
            double diff = fabs((grad_ref[i] - grad_app[i]) /
                               (grad_ref[i] ? grad_ref[i] : 1.));
            if (diff > max) max = diff;
        }

        avg += max;
        if (max > total_max) total_max = max;
        n++;
    }

    avg /= (double)100;

    // check hessian
    avg = 0;
    total_max = 0;
    k = 0;
    for (unsigned int n = 0 ; n < 100 ; k++) {
        std::cout << "computing hessian at point #" << k
                  << " (" << n << ")" << std::endl;
        vec3 random(drand48(), drand48(), drand48());
        vec3 p = box.min() + random * box.size();

        // make sure we are at a valid position within the volume
        if (wrapper.confidence(p) < 0.5) continue;

        // compute gradient
        vec6 hess_ref = wrapper.hessian(p);
        vec6 hess_app = _hessian(wrapper, p);
        if (wrong) continue;

        double scale = hess_ref[0] / hess_app[0];
        hess_app *= scale;

        std::cout << "gage: " << hess_ref << std::endl;
        std::cout << "approx: " << hess_app << " (" << N << ")" << std::endl;
        std::cout << "scale: " << scale << std::endl;

        double max = 0;
        for (unsigned int i = 0 ; i < 3 ; i++) {
            double delta = norm(hess_ref - hess_app);
            double denom = norm(hess_ref);
            if (denom == 0) denom = 1.;
            double diff = fabs(delta / denom);
            if (diff > max) max = diff;
        }

        avg += max;
        if (max > total_max) total_max = max;
        n++;
    }

    avg /= (double)100;

    std::cout << "average relative error = " << avg << std::endl
              << "maximum relative error = " << total_max << std::endl;
}



