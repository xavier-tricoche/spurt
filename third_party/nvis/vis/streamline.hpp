#ifndef __streamline_hpp
#define __streamline_hpp

#include "math/fixed_vector.hpp"
#include "math/bezier_spline.hpp"
#include "math/dopri5.hpp"

#include <iosfwd>
#include <vector>

namespace nvis
{

class streamline: public bezier_spline<vec4, 4>
{
public:

	typedef nvis::fixed_vector<double, 5>  int_var;
	typedef nvis::dopri5<int_var>         int_type;
	typedef int_type::step                int_step;

	typedef bezier_spline<vec4, 4>         base;

	enum {
		ARC_INDEX    = 3,
		DOMAIN_INDEX = 4,
	};

	// streamline integration result
	enum state {
		OK = 0,
		LEFT_DOMAIN,
		ERROR,
		CUSTOM,
	};

	streamline(const nvis::vec3& p0 = nvis::vec3(0, 0, 0),
	           const double& t0 = 0.0, const double& a0 = 0.0) :
			_state(OK) {
		bwd.t = fwd.t = start_t = t0;
		bwd.y = fwd.y = int_var(p0[0], p0[1], p0[2], a0, 0.0);

		abstol = 1e-7;
		reltol = 1e-7;
		stepsz = 0.0;
		record = true;

		_knots.push_back(t0);
		_ctrlp.push_back(subv<0, 4>(fwd.y));
	}

	streamline(std::istream& in) {
		in.read((char*)&start_t, sizeof(start_t));
		base::read(in);

		if (!in)
			return;

		bwd.t = _knots.front();
		fwd.t = _knots.back();

		const nvis::vec4& p0 = _ctrlp.front();
		const nvis::vec4& p1 = _ctrlp.back();

		bwd.y = int_var(p0[0], p0[1], p0[2], p0[3], 0.0);
		fwd.y = int_var(p1[0], p1[1], p1[2], p1[3], 0.0);

		abstol = 1e-7;
		reltol = 1e-7;
		stepsz = 0.0;
		record = true;
	}

	~streamline() {
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
		if (t < _knots.front()) {
			do_advance(bwd, rhs, t, stop);

			_ctrlp.insert(_ctrlp.begin(), lctrlp.rbegin(), lctrlp.rend());
			_knots.insert(_knots.begin(), lknots.rbegin(), lknots.rend());
		}
		else if (t > _knots.back()) {
			do_advance(fwd, rhs, t, stop);

			_ctrlp.insert(_ctrlp.end(), lctrlp.begin(), lctrlp.end());
			_knots.insert(_knots.end(), lknots.begin(), lknots.end());
		}

		std::vector<vec4>()  .swap(lctrlp);
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

			int_type::result result = integrator.do_step(wrapper, step);

			// integration error?
			if (result != int_type::OK && result != int_type::T_MAX_REACHED)
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
			if (result == int_type::T_MAX_REACHED)
				return _state = OK;
		}
	}

	void append(const int_step& s) {
		vec4 p[4] = {
			subv<0, 4>(s._p[0] + (s._p[1] + s._p[2]) / 4),
			subv<0, 4>(s._p[0] + s._p[1] / 2 + s._p[2] / 3 + (s._p[3] + s._p[4]) / 6),
			subv<0, 4>(s._p[0] + (3*s._p[1] + s._p[2] + s._p[3]) / 4),
			subv<0, 4>(s._p[0] + s._p[1]),
		};

		lknots.insert(lknots.end(), s._t1);
		lctrlp.insert(lctrlp.end(), p, p + 4);
	}

	const double& a_min() const {
		return _ctrlp.front()[3];
	}

	const double& a_max() const {
		return _ctrlp.back()[3];
	}

	const double& t_start() const {
		return start_t;
	}

	nvis::vec3 operator()(const double& t) const {
		return subv<0, 3>(base::operator()(t));
	}

	double t(const double& a) const {
		vec4 result;

		if (a > a_max())
			return t_min();

		if (a < a_min())
			return t_max();

		double tmin = t_min(), tmax = t_max(), tmid = (tmin + tmax) / 2;
		double amin = a_min(), amax = a_max(), amid = base::operator()(tmid)[3];

		while (fabs(amid - a) / (amax - amin) > 1e-5) {
			if (a < amid)
				tmax = tmid;
			else
				tmin = tmid;

			tmid = (tmin + tmax) / 2;
			amid = base::operator()(tmid)[3];
		}

		return tmid;
	}

	double a(const double& t) const {
		return base::operator()(t)[3];
	}

	nvis::vec3 operator[](const double& a) const {
		return subv<0, 3>(base::operator()(t(a)));
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
			nvis::vec3 f, v(p[0], p[1], p[2]);

			bool success = i(t, v, f);

			if (success)
				return int_var(f[0], f[1], f[2], norm(f), 0.0);
			else
				return int_var(0, 0, 0, 0, 1.0);
		}


	};

protected:

	// --- variables ---

	double      start_t;
	int_type    fwd, bwd;
	state       _state;

	std::vector<double> lknots;
	std::vector<vec4>   lctrlp;

	// --- streamline parameters ---

public:

	double abstol;
	double reltol;
	double stepsz;
	bool   record;

	friend std::ostream& operator<<(std::ostream& out, const streamline& s);
};

// -----------------------------------------

inline std::ostream& operator<<(std::ostream& out, const streamline& s)
{
	s.write(out);
	return out;
}

// --------------------------------------------------------------------------

} // namespace nvis

#endif // __streamline_hpp








