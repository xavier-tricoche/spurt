#ifndef __poincare_map_hpp
#define __poincare_map_hpp

#include <vector>
#include <exception>
#include <iostream>

#include <nvis-math/fixed_vector.hpp>
#include <nvis-math/bounding_box.hpp>

#include "tokamak_field.hpp"
#include "map2d.hpp"

class poincare_map: public map2d
{
public:

	enum Mode {
		RK56,
		ADAMS
	};

	poincare_map(const tokamak_field* field) :
			_field(field), _prec(1e-6), _mode(RK56), _verbose_jacobian(false) {
	}

	poincare_map(const poincare_map& other) {
		_prec  = other._prec;
		_mode  = other._mode;
		_field = other._field->clone();
		_verbose_jacobian = other._verbose_jacobian;
	}

	void mode(Mode mode) {
		_mode = mode;
	}

	void precision(const double& prec) {
		_prec = prec;
	}
	
	double get_precision() const {
		return _prec;
	}

	void map_complete(const nvis::vec2& in, std::vector<value_type>& out, int niter) const {
		if (_mode == RK56)
			map_and_jac_RK56(in, out, niter);
		else
			map_ADAMS(in, out, niter);
	}

	void map(const nvis::vec2& in, std::vector<nvis::vec2>& out, int niter) const {
		std::vector<value_type> __out;
		map_complete(in, __out, niter);
		out.resize(__out.size());
		for (unsigned int i = 0 ; i < __out.size() ; ++i) {
			out[i] = __out[i].x;
		}
	}

	void jmap_complete(const nvis::vec2& in, std::vector<value_type>& out, int niter) const {
		if (_mode != RK56) {
			std::cerr << "jacobian integration undefined for Adams-Bashforth\n";
			throw map_undefined();
		}
		map_and_jac_RK56(in, out, niter);
	}

	void jmap(const nvis::vec2& in, std::vector< nvis::vec2 >& hits,
	          std::vector< nvis::mat2 >& jacobians, int niter) const {
		std::vector<value_type> __out;
		jmap_complete(in, __out, niter);
		hits.resize(__out.size());
		jacobians.resize(__out.size());
		for (unsigned int i = 0 ; i < __out.size() ; ++i) {
			hits[i] = __out[i].x;
			jacobians[i] = __out[i].J;
		}
	}

	value_type map_complete(const nvis::vec2& in, int niter) const {
		std::vector<value_type> out;
		map_complete(in, out, niter);
		return out.back();
	}

	nvis::vec2 map(const nvis::vec2& in, int niter) const {
		std::vector<value_type> out;
		map_complete(in, out, niter);
		if (out.size() != std::abs(niter)) throw map_undefined();
		return out.back().x;
	}

	value_type jmap_complete(const nvis::vec2& in, int niter) const {
		std::vector<value_type> out;
		jmap_complete(in, out, niter);
		return out.back();
	}

	void jmap_complete(const nvis::vec2& in, value_type& hit, int niter) const {
		std::vector<value_type> out;
		jmap_complete(in, out, niter);
		hit = out.back();
	}

	void jmap(const nvis::vec2& in, nvis::vec2& hit, nvis::mat2& jacobian, int niter) const {
		std::vector<value_type> out;
		jmap_complete(in, out, niter);
		if (out.size() != std::abs(niter)) throw map_undefined();
		hit = out.back().x;
		jacobian = out.back().J;
	}

	nvis::bbox2 bounds() const {
		return _field->bounds();
	}

	nvis::bvec2 periodic() const {
		return nvis::bvec2(false, false);
	}

	poincare_map* clone() const {
		return new poincare_map(*this);
	}

	const tokamak_field* field() const {
		return _field;
	}

	void set_verbose(bool verbose = true) {
		_verbose_jacobian = verbose;
	}

private:

	void map_RK56(const nvis::vec2& in, std::vector<value_type>& out, int niter) const;
	void map_and_jac_RK56(const nvis::vec2& in, std::vector<value_type>& out, int niter) const;
	void map_ADAMS(const nvis::vec2& in, std::vector<value_type>& out, int niter) const {};

	const tokamak_field* _field;
	double               _prec;
	Mode                 _mode;
	bool 				 _verbose_jacobian;
};

#endif // __poincare_map_hpp
















