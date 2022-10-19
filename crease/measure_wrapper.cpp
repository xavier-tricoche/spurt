#include <iostream>

#include <math/bounding_box.hpp>

#include <crease/measure_wrapper.hpp>
#include <image/nrrd_wrapper.hpp>

using namespace spurt;
using namespace nvis;
using namespace std;

bool wrong;

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
                double delta = spurt::norm(dg_cur - dg_prev);
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

    bbox3 box = spurt::bounds<3>(nrrd);

    // check gradient
    double avg = 0, total_max = 0;
    unsigned int k = 0;
    for (unsigned int n = 0 ; n < 100 ; k++) {
        std::cout << "computing gradient at point #" << k
                  << " (" << n << ")" << std::endl;
        vec3 random(drand48(), drand48(), drand48());
        vec3 p = box.min() + random * box.size();

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
            double delta = spurt::norm(hess_ref - hess_app);
            double denom = spurt::norm(hess_ref);
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



