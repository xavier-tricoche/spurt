#ifndef __integral_curve_hpp
#define __integral_curve_hpp

#include "nvis-math/fixed_vector.hpp"
#include "nvis-math/bezier_spline.hpp"
#include "nvis-math/dopri5.hpp"

#include <iosfwd>
#include <vector>

namespace nvis
{
template<int N>
class integral_curve: public bezier_spline<fixed_vector<double, N+1>, 4>
{
public:

    typedef nvis::fixed_vector<double, N>		pos_type;
    typedef nvis::fixed_vector<double, N+1 >	vec_type;
    typedef nvis::fixed_vector<double, N+2 >  	int_var;
    typedef nvis::dopri5< int_var >      		int_type;

    typedef typename int_type::step             int_step;
    typedef bezier_spline<vec_type, 4>        	base;
    typedef integral_curve<N>					self_type;

    enum {
        ARC_INDEX    = N,
        DOMAIN_INDEX = N+1,
    };

    // integral_curve integration result
    enum state {
        OK = 0,
        LEFT_DOMAIN,
        ERROR,
        CUSTOM,
    };

    static int_var to_int_var(const pos_type& x) {
        int_var y(0);
        for (int i=0 ; i<N ; ++i) y[i] = x[i];
		return y;
    }

    static int_var to_int_var(const vec_type& x) {
        int_var y(0);
        for (int i=0 ; i<N+1 ; ++i) y[i] = x[i];
		return y;
    }

    integral_curve(const pos_type& p0 = pos_type(0),
                   const double& t0 = 0.0, const double& a0 = 0.0) :
        _state(OK) {
        bwd.t = fwd.t = start_t = t0;
        int_var tmp = to_int_var(p0);
        tmp[N] = a0;
        bwd.y = fwd.y = tmp;

        abstol = 1e-7;
        reltol = 1e-7;
        stepsz = 0.0;
        record = true;

        this->_knots.push_back(t0);
        this->_ctrlp.push_back(subv<0, N+1>(fwd.y));
    }

    integral_curve(std::istream& in) {
        in.read((char*)&start_t, sizeof(start_t));
        base::read(in);

        if (!in)
            return;

        bwd.t = this->_knots.front();
        fwd.t = this->_knots.back();

        const vec_type& p0 = this->_ctrlp.front();
        const vec_type& p1 = this->_ctrlp.back();

        bwd.y = int_var(to_int_var(p0));
        fwd.y = int_var(to_int_var(p1));

        abstol = 1e-7;
        reltol = 1e-7;
        stepsz = 0.0;
        record = true;
    }

    ~integral_curve() {
    }

    // ---

    struct no_stop {
        bool operator()(const int_step&) {
            return false;
        }
    };

    template<typename RHS>
    state advance(const RHS& rhs, double t) {
        no_stop ns;
        return advance(rhs, t, ns);
    }

    template<typename RHS, typename STOP>
    state advance(const RHS& rhs, double t, STOP& stop) {
        if (t < this->_knots.front()) {
            do_advance(bwd, rhs, t, stop);

            this->_ctrlp.insert(this->_ctrlp.begin(), lctrlp.rbegin(), lctrlp.rend());
            this->_knots.insert(this->_knots.begin(), lknots.rbegin(), lknots.rend());
        }
        else if (t > this->_knots.back()) {
            do_advance(fwd, rhs, t, stop);

            this->_ctrlp.insert(this->_ctrlp.end(), lctrlp.begin(), lctrlp.end());
            this->_knots.insert(this->_knots.end(), lknots.begin(), lknots.end());
        }

        std::vector<vec_type>().swap(lctrlp);
        std::vector<double>().swap(lknots);

        return _state;
    }

    template<typename RHS, typename STOP>
    state do_advance(int_type& integrator, const RHS& rhs, double t, STOP& stop) {	
        rhs_wrapper<RHS> wrapper(rhs);

        integrator.abstol = abstol;
        integrator.reltol = reltol;

        int_step step;

        for (;;) {
            int_type int_save = integrator;
            integrator.t_max = t;
            integrator.h     = stepsz;

            typename int_type::result res = integrator.do_step(wrapper, step);

            // integration error?
            if (res != int_type::OK && res != int_type::T_MAX_REACHED)
                return _state = ERROR;

            if (stop(step)) {
                append(step); // final step
                return _state = CUSTOM;
            }

            // domain left?
            if (integrator.y[DOMAIN_INDEX] > 0.0) {
                double likely_end = integrator.t - integrator.y[DOMAIN_INDEX];

                if (likely_end > 10.0*abstol) {
                    int_save.t_max = likely_end;
                    int_save.do_step(wrapper, step);
                    append(step); // final step
                }

                return _state = LEFT_DOMAIN;
            }

            // record step
            if (record) append(step);

            // t reached?
            if (res == int_type::T_MAX_REACHED)
                return _state = OK;
        }
    }

    void append(const int_step& s) {
        vec_type p[4] = {
            subv<0, N+1>(s._p[0] + (s._p[1] + s._p[2]) / 4),
            subv<0, N+1>(s._p[0] + s._p[1] / 2 + s._p[2] / 3 + (s._p[3] + s._p[4]) / 6),
            subv<0, N+1>(s._p[0] + (3*s._p[1] + s._p[2] + s._p[3]) / 4),
            subv<0, N+1>(s._p[0] + s._p[1]),
        };

        lknots.insert(lknots.end(), s._t1);
        lctrlp.insert(lctrlp.end(), p, p + 4);
    }

    const vec_type& front() const {

    }

    const double& a_min() const {
        return this->_ctrlp.front()[N];
    }

    const double& a_max() const {
        return this->_ctrlp.back()[N];
    }

    const double& t_start() const {
        return start_t;
    }

    pos_type operator()(const double& t) const {
        return subv<0, N>(base::operator()(t));
    }

    double t(const double& a) const {
        vec_type result;

        if (a > a_max())
            return this->t_min();

        if (a < a_min())
            return this->t_max();

        double tmin_ = this->t_min(), tmax_ = this->t_max(), tmid_ = (tmin_ + tmax_) / 2;
        double amin_ = this->a_min(), amax_ = this->a_max(), amid_ = base::operator()(tmid_)[N];

        while (fabs(amid_ - a) / (amax_ - amin_) > 1e-5) {
            if (a < amid_)
                tmax_ = tmid_;
            else
                tmin_ = tmid_;

            tmid_ = (tmin_ + tmax_) / 2;
            amid_ = base::operator()(tmid_)[N];
        }

        return tmid_;
    }

    double a(const double& t) const {
        return base::operator()(t)[N];
    }

    pos_type operator[](const double& a) const {
        return subv<0, N>(base::operator()(t(a)));
    }

private:

    void write(std::ostream& out) const {
        out.write((const char*)&start_t, sizeof(start_t));
        base::write(out);
    }

    // --- integrand wrapper helper class ---

    template<typename I> struct rhs_wrapper {
        typedef I int_type;

        const int_type& i;

        rhs_wrapper(const int_type& i_) : i(i_) {
        }

        int_var operator()(const double& t, const int_var& p) const {
            pos_type f, v(subv<0, N>(p));

            bool success = i(t, v, f);

            if (success)
            {
                int_var iv = to_int_var(f);
                iv[N] = nvis::norm(f);
                return iv;
            }
            else {
                int_var iv(0);
                iv[N+1] = 1.0;
                return iv;
            }
        }


    };

protected:

    // --- variables ---

    double      start_t;
    int_type    fwd, bwd;
    state       _state;

    std::vector<double> lknots;
    std::vector<vec_type>   lctrlp;

    // --- integral_curve parameters ---

public:

    double abstol;
    double reltol;
    double stepsz;
    bool   record;

    template<int K>
    friend std::ostream& operator<<(std::ostream& out, const integral_curve<K>& s);
};

// -----------------------------------------

template<int N>
inline std::ostream& operator<<(std::ostream& out, const integral_curve<N>& s)
{
    s.write(out);
    return out;
}

// --------------------------------------------------------------------------

} // namespace nvis

#endif // __integral_curve_hpp








