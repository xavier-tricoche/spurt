#ifndef __dopri5_hpp
#define __dopri5_hpp

#include <cmath>
#include <iostream>
#include <stdexcept>

// beware of commenting out the QUIET define in multithreaded implementations
#define QUIET


namespace
{

const double safe = 0.9;
const double epsilon = std::numeric_limits<double>::epsilon();
const double facl = 0.2;
const double facr = 10.0;
const double beta = 0.04;
const unsigned int nstiff = 100;

const double
    a21 = 0.2,
    a31 = 3.0 / 40.0, a32 = 9.0 / 40.0,
    a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0,
    a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0,
    a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
    a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;

const double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0;

const double d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
             d4 = -10690763975.0 / 1880347072.0, d5 = 701980252875.0 / 199316789632.0,
             d6 = -1453857185.0 / 822651844.0, d7 = 69997945.0 / 29380423.0;

const double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
             e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

double sign(double a, double b)
{
    return (b > 0.0) ? fabs(a) : -fabs(a);
}
}

namespace nvis
{

struct invalid_position_exception: public std::runtime_error {
    typedef std::runtime_error base_type;
    explicit invalid_position_exception(const std::string& what_arg) : base_type(what_arg), on_start(false) {}
    explicit invalid_position_exception(const char* what_arg ) : base_type(what_arg), on_start(false) {}
    explicit invalid_position_exception() : base_type("Invalid position queried for interpolation"), on_start(false) {}

    bool on_start;
};

template<typename V> struct dopri5 {
    typedef V value_type;

    enum result {
        OK,
        T_MAX_REACHED,
        STEPSIZE_UNDERFLOW,
        STIFFNESS_DETECTED,
    };

    double reltol;
    double abstol;

    double h;
    double h_max;

    unsigned int n_accepted;
    unsigned int n_rejected;
    unsigned int n_steps;
    unsigned int n_eval;

    double     t, t_max;
    value_type y, k1;

    // stepsize control stabilization
    double     facold;

    // stiffness detection
    double     hlamb;
    int        iasti;
    int        nonsti;

    // error monitoring
    value_type total_ee;

public:
    struct step {
        value_type _p[5];
        value_type _y0, _y1;
        double     _t0, _t1;

        value_type y(const double& t) const {
            const double a = (t - _t0) / (_t1 - _t0);
            const double b = 1.0 - a;

            return _p[0] + a*(_p[1] + b*(_p[2] + a*(_p[3] + b*_p[4])));
        }

        double y(const double& t, unsigned int n) const {
            const double a = (t - _t0) / (_t1 - _t0);
            const double b = 1.0 - a;

            return _p[0][n] + a*(_p[1][n] + b*(_p[2][n] + a*(_p[3][n] + b*_p[4][n])));
        }

        value_type dy(const double& t) const {
            double a = (t - _t0) / (_t1 - _t0);
            double b = 1.0 - a;

            value_type pc = _p[3] + b * _p[4];

            return (_p[1] + b*(_p[2] + a*pc) -
                    a*(_p[2] + a*pc -  b*(_p[3] + (b - a)*_p[4]))) / (_t1 - _t0);
        }

        double dy(const double& t, unsigned int n) const {
            double a = (t - _t0) / (_t1 - _t0);
            double b = 1.0 - a;

            double pc = _p[3][n] + b * _p[4][n];

            return (_p[1][n] + b*(_p[2][n] + a*pc) -
                    a*(_p[2][n] + a*pc -  b*(_p[3][n] + (b - a)*_p[4][n]))) / (_t1 - _t0);
        }

        const value_type& y0() const {
            return _y0;
        }
        const value_type& y1() const {
            return _y1;
        }

        const double& t0() const {
            return _t0;
        }
        const double& t1() const {
            return _t1;
        }
    };

    dopri5() :
            reltol(1e-10), abstol(1e-10), h(0.0), h_max(0.0),
            n_accepted(0), n_rejected(0), n_steps(0), n_eval(0),
            facold(1e-4), hlamb(0.0), iasti(0), total_ee(0) {

    }

    value_type error_estimate() const {
        return total_ee;
    }

    template<typename RHS> double h_init(const RHS& rhs,
                                         const double& h_max) {
        double direction = sign(1.0, t_max - t);

        double sk, sqr;
        double dnf = 0.0;
        double dny = 0.0;

        double h;

        for (int i = 0; i < V::size(); i++) {
            sk = abstol + reltol * fabs(y[i]);
            sqr = k1[i] / sk;
            dnf += sqr * sqr;
            sqr = y[i] / sk;
            dny += sqr * sqr;
        }

        if ((dnf <= 1.0e-10) || (dny <= 1.0e-10))
            h = 1.0e-6;
        else
            h = sqrt(dny / dnf) * 0.01;

        h = std::min(h, h_max);
        h = sign(h, direction);

        // perform an explicit Euler step
        value_type k3 = y + h * k1;

        value_type k2 = rhs(t + h, k3);
        n_eval++;

        // estimate the second derivative of the solution
        double der2 = 0.0;

        for (int i = 0; i < V::size(); i++) {
            sk = abstol + reltol * fabs(y[i]);
            sqr = (k2[i] - k1[i]) / sk;
            der2 += sqr * sqr;
        }

        der2 = sqrt(der2) / h;

        // step size is computed such that
        // h**(1/5) * max( norm(k1), norm(der2) ) = 0.01
        double der12 = std::max(fabs(der2), sqrt(dnf));
        double h1;

        if (der12 <= 1.0e-15)
            h1 = std::max(1.0e-6, fabs(h) * 1.0e-3);
        else
            h1 = pow(0.01 / der12, 0.2);

        h = std::min(100.0 * h, std::min(h1, h_max));

        return sign(h, direction);
    }

    template<typename RHS> result do_step(const RHS& rhs,
                                          step& _step, bool verbose=false) {
        const double direction = sign(1.0, t_max - t);

        value_type k2, k3, k4, k5, k6, k7;
        bool reject = false;

        // compute maximum stepsize
        double local_h_max = h_max;

        if (local_h_max == 0.0)
            local_h_max = fabs(t_max - t);

        // compute k1 to ensure first-same-as-last principle,
        // maybe also needed for hinit())
        if (n_steps == 0) {
            try {
                k1 = rhs(t, y);
            }
            catch(invalid_position_exception& e) {
                e.on_start = true;
                throw;
            }
            n_eval++;
        }

        // determine stepsize (user-specified or educated guess)
        if (h == 0.0)
            try {
                h = h_init(rhs, local_h_max);
            }
            catch(invalid_position_exception& e) {
                e.on_start = true;
                throw;
            }
        else
            h = sign(h, direction);

        // integration step loop
        while (true) {
            bool last = false;
            value_type y_new, y_stiff;

            // stepsize underflow?
            if (0.1*fabs(h) <= fabs(t)*epsilon || fabs(h) < 0.1*abstol) {
#ifndef QUIET
                if (verbose) {
                    std::cout << "dopri5(): exiting at t = " << t
                    << ", step size too small (h = "
                    << h << ")\n";
                }
#endif

                return STEPSIZE_UNDERFLOW;
            }

            // do not run past integration end
            if ((t + 1.01*h - t_max) * direction > 0.0) {
                last = true;
                h = t_max - t;
            }

            n_steps++;
            // std::cout << "dopri5 step: t = " << t
            //       << ", y = " << y
            //       << ", h = " << h
            //       << ", t+h = " << t+h << '\n';

            // perform stages
            bool invalid_query = false;
            try {
                y_new = y + h * a21 * k1;
                k2 = rhs(t + c2 * h, y_new);

                y_new = y + h * (a31 * k1 + a32 * k2);
                k3 = rhs(t + c3 * h, y_new);

                y_new = y + h * (a41 * k1 + a42 * k2 + a43 * k3);
                k4 = rhs(t + c4 * h, y_new);

                y_new = y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4);
                k5 = rhs(t + c5 * h, y_new);

                y_stiff = y_new = y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5);
                k6 = rhs(t + h, y_new);

                y_new = y + h * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6);
                k7 = rhs(t + h, y_new);
            }
            catch(invalid_position_exception& e) {
                invalid_query = true;
#ifndef QUIET
                if (verbose) {
                    std::cout << "dopri5::do_step: invalid query\n";
                }
#endif
            }

            n_eval += 6;

            double err = 0.0, h_new, fac11;

            // error estimation
            value_type ee = h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7);
            double sk, sqr;

            for (int i = 0; i < V::size(); i++) {
                sk = abstol + reltol * std::max(fabs(y[i]), fabs(y_new[i]));
                sqr = ee[i] / sk;
                err += sqr * sqr;
            }

            err = sqrt(err / (double)V::size());

            // compute next (potential) stepsize
            fac11 = pow(err, 0.2 - beta * 0.75);

            // Lund-stabilization
            double fac = fac11 / pow(facold, beta);

            // we require facl <= h_new/h <= facr
            fac = std::max(1.0 / facr, std::min(1.0 / facl, fac / safe));

            h_new = h / fac;

            if (err <= 1.0 && !invalid_query) {
                // step accepted
                facold = std::max(err, 1.0e-4);

                n_accepted++;
                total_ee += nvis::abs(ee);

                // stiffness detection
                if (!(n_accepted % nstiff) || (iasti > 0)) {
                    double stnum = 0.0, stden = 0.0, sqr;

                    for (int i = 0; i < V::size(); i++) {
                        sqr = k7[i] - k6[i];
                        stnum += sqr * sqr;
                        sqr = y_new[i] - y_stiff[i];
                        stden += sqr * sqr;
                    }

                    if (stden > 0.0)
                        hlamb = h * sqrt(stnum / stden);

                    if (hlamb > 3.25) {
                        nonsti = 0;
                        iasti++;

                        if (iasti == 15) {
#ifndef QUIET
                        if (verbose) {
                            std::cout << "dopri5(): exiting at t = " << t
                            << ", problem seems stiff (y = "
                            << y << ")\n";
                        }
#endif

                            return STIFFNESS_DETECTED;
                        }
                    }
                    else {
                        nonsti++;
                        if (nonsti == 6)
                            iasti = 0;
                    }
                }

                // --- step looks ok - prepare for return

                if (reject)
                    h_new = direction * std::min(fabs(h_new), fabs(h));

                // fill in step
                // make interpolation polynomial
                _step._p[0] = y;
                _step._p[1] = y_new - y;
                _step._p[2] = h * k1 - _step._p[1];
                _step._p[3] = -h * k7 + _step._p[1] - _step._p[2];
                _step._p[4] = h * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6 + d7 * k7);

                _step._t0 = t;
                _step._y0 = y;

                _step._t1 = t = t + h;
                _step._y1 = y = y_new;

                h = h_new;

                // update internal state
                // first-same-as-last for k1
                k1 = k7;

                // normal exit
                return last ? T_MAX_REACHED : OK;
            }
            else {
                // step rejected
                if (invalid_query) {
                    h_new = 0.2*h;
#ifndef QUIET
                    if (verbose) {
                        std::cout << "dopri5: after invalid query, step size reduced fron "
                        << std::setprecision(std::numeric_limits<double>::digits10+1) << h << " to " << h_new << '\n';
                    }
#endif
                }
                else h_new = h / std::min(1.0 / facl, fac11 / safe);
                reject = true;

                if (n_accepted >= 1)
                    n_rejected++;

                h = h_new;
            }
        }
    }
};


} // namespace nvis

#endif // __dopri5_hpp
