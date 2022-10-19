#include "poincare_map.hpp"

#include <math/dopri5.hpp>
#include <assert.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <fstream>

using namespace nvis;

// --------------------------------------------------------------------------

void poincare_map::map_RK56(const vec2& in, std::vector<value_type>& out, int niter) const
{
	out.clear();
	out.reserve(std::abs(niter));

	dopri5<vec3> intg;

	intg.t = 0;
	intg.y = _field->unproject(in);
	intg.t_max = niter < 0 ? -1e9 : 1e9;
	intg.h = 1e-1;
	intg.abstol = _prec;
	intg.reltol = _prec;
	// 
	// std::cerr << "RK56 starts at " << in << '\n';

	dopri5<vec3>::result res;
	dopri5<vec3>::step   step;

	unsigned int nsteps = 0;
	while (out.size() < std::abs(niter)) {
		++nsteps;
		// 
		// std::cerr << "\t\tRK56 at " << step.y0() << '\n';

		try {
			res = intg.do_step(*_field, step);

			if (res != dopri5<vec3>::OK) {				// 
							// std::cerr << "dopri returned an error code: " << res << '\n';
				break;
			}

			if (_field->intersect(step.y0(), step.y1())) {
				double t0 = step.t0();
				double t1 = step.t1();

				vec3 y0 = step.y0();
				vec3 y1 = step.y1();

				double tmid;
				vec3   ymid;

				for (unsigned int i = 0; i < 10; ++i) {
					tmid = (t0 + t1) / 2;
					ymid = step.y(tmid);

					if (_field->intersect(y0, ymid)) {
						t1 = tmid;
						y1 = ymid;
					}
					else if (_field->intersect(ymid, y1)) {
						t0 = tmid;
						y0 = ymid;
					}
					else
						assert(false);
				}

				value_type v;
				v.x = _field->project(ymid);
				v.t = tmid;
				v.err = _field->project(intg.error_estimate());
				v.nsteps = intg.n_accepted;
				v.total_nsteps = intg.n_steps;
				out.push_back(v);
			}
		}
		catch (...) {
			// std::cout << "exception caught in map_RK56\n";
			break;
		}
	}
}

nvis::mat3 integrate_jacobian(const nvis::mat3& w, const tokamak_field* field,
                              const dopri5<vec3>::step& step, double tmax,
                              bool verbose = false)
{
	nvis::mat3 _w = w;
	double t0 = step.t0();
	double t1 = tmax;
	double dt = (t1 - t0) / 10.;
	nvis::mat3 J;
	for (int i = 0 ; i < 10 ; ++i) {
		double t = t0 + (double)i * dt;
		nvis::vec3 x = step.y(t);
		J = field->derivative(0, x);
		_w += dt * (J * _w);
	}
	if (verbose) std::cerr << "integrating Jacobian from " << w << " to " << _w << '\n';
	return _w;
}

void poincare_map::map_and_jac_RK56(const vec2& in, std::vector<value_type>& out, int niter) const
{
	out.clear();
	out.reserve(std::abs(niter));

	dopri5<vec3> intg;

	intg.t = 0;
	intg.y = _field->unproject(in);
	intg.t_max = niter < 0 ? -1e9 : 1e9;
	intg.h = 1e-1;
	intg.abstol = _prec;
	intg.reltol = _prec;

	dopri5<vec3>::result res;
	dopri5<vec3>::step   step;

	nvis::mat3 w, dvdx, aux;
	w = nvis::mat3::identity(); // derivative of flow map with respect to initial condition at time 0
	unsigned int nsteps = 0;
	while (out.size() < std::abs(niter)) {
		++nsteps;

		try {
			res = intg.do_step(*_field, step);

			if (res != dopri5<vec3>::OK)
				break;

			// std::cout << "step #" << nsteps << ": " << step.y0() << " -- " << step.y1() << std::endl;

			if (_field->intersect(step.y0(), step.y1())) {
				double t0 = step.t0();
				double t1 = step.t1();

				vec3 y0 = step.y0();
				vec3 y1 = step.y1();

				double tmid;
				vec3   ymid;

				for (unsigned int i = 0; i < 10; ++i) {
					tmid = (t0 + t1) / 2;
					ymid = step.y(tmid);

					if (_field->intersect(y0, ymid)) {
						t1 = tmid;
						y1 = ymid;
					}
					else if (_field->intersect(ymid, y1)) {
						t0 = tmid;
						y0 = ymid;
					}
					else
						assert(false);
				}

				value_type v;
				v.x = _field->project(ymid);
				v.t = tmid;
				v.err = _field->project(intg.error_estimate());
				aux = integrate_jacobian(w, _field, step, tmid, _verbose_jacobian);
				v.J = _field->project(aux);
				if (_verbose_jacobian) {
					std::cerr << "projecting integrated Jacobian " << aux
					          << " to " << v.J << 'n';
				}
				v.nsteps = intg.n_accepted;
				v.total_nsteps = intg.n_steps;
				out.push_back(v);
			}

			// // Euler step on w
			// double dt = step.t1() - step.t0();
			// w += dt * (dvdx * w);
			// // update Jacobian
			// dvdx = _field->derivative(0, intg.y);
			w = integrate_jacobian(w, _field, step, step.t1(), _verbose_jacobian);

			// std::cout << "w = " << w << ", dvdx = " << dvdx << ", dt = " << dt << '\n';
			// std::cout << "det w = " << nvis::det(w) << '\n';
		}
		catch (...) {
			break;
		}
	}
}