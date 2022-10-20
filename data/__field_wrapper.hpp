#ifndef __FIELD_WRAPPER_HPP__
#define __FIELD_WRAPPER_HPP__

#include <stdexcept>
#include <iostream>
#include <vector>
#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>

namespace spurt
{
template<typename _Format> class data_traits {};

template<>
class data_traits<Nrrd*>
{	
public:
	typedef double 								scalar_type;
	typedef nvis::fixed_vector<scalar_type, 3>	vector_type;
	typedef nvis::vec3							point_type;

	data_traits(const Nrrd* nin) {
		if (nin->dim < 3)
			throw std::runtime_error("data_traits<Nrrd*>: invalid dimension");
		for (int i = 0 ; i < 3 ; ++i) __size[i] = nin->axis[nin->dim-3+i].size;
		__bounds = spurt::bounds<3>(nin);
		std::cerr << "bounds are " << __bounds << "\n";
		__step = __bounds.size() / point_type(__size - nvis::ivec3(1, 1, 1));
		__scalar_values = spurt::to_array<scalar_type>(nin);
		__offset[0] = 0;
		__offset[1] = 1;
		__offset[2] = 1+__size[0];
		__offset[3] = __size[0];
		for (int i=0 ; i<4 ; ++i) {
			__offset[4+i] = __offset[i] + __size[1];
		}
	}
	
	~data_traits() {
		delete[] __scalar_values;
	}

	bool get_value(const point_type& x, scalar_type& f) const {
		point_type y;
		if (!g2l(y, x)) return false;
		nvis::ivec3 id(floor(y[0]), floor(y[1]), floor(y[2]));
		point_type z = y - point_type(id);
		scalar_type u = z[0], v = z[1], w = z[2];
		scalar_type U = 1. - u, V = 1. - v, W = 1. - w;
		scalar_type weight[8] = {
			U*V*W, u*V*W, u*v*W, U*v*W,
			U*V*w, u*V*w, u*v*w, U*v*w
		};
		int n = index(id);
		f = 0;
		for (int i=0 ; i<8 ; ++i) f += weight[i] * __scalar_values[n+__offset[i]];

		return true;
	}

	bool get_value(const point_type& x, vector_type& f) const {
		point_type y;
		if (!g2l(y, x)) return false;
		nvis::ivec3 id(floor(y[0]), floor(y[1]), floor(y[2]));
		point_type z = y - point_type(id);
		scalar_type u = z[0], v = z[1], w = z[2];
		scalar_type U = 1. - u, V = 1. - v, W = 1. - w;
		scalar_type weight[8] = {
			U*V*W, u*V*W, u*v*W, U*v*W,
			U*V*w, u*V*w, u*v*w, U*v*w
		};
		int n = index(id);
		f = vector_type(0);
		for (int i=0 ; i<8 ; ++i) {
			f += weight[i] * nvis::suba<scalar_type, 3>(__scalar_values, n+__offset[i]);
		}
		return true;
	}

	const nvis::bbox3& bounds() const {
		return __bounds;
	}

private:

	bool g2l(point_type& l, const point_type& g) const {
		if (!__bounds.inside(g)) return false;
		l = (g - __bounds.min()) / __step;
		return true;
	}

	int index(const nvis::ivec3& id) const {
		return id[0] + __size[0]*(id[1] + __size[1]*id[2]);
	}

	nvis::ivec3 							__size;
	nvis::bbox3 							__bounds;
	point_type 								__step;
	scalar_type*							__scalar_values;
	int										__offset[8];
};

template<typename _Data_Traits>
struct rhs_wrapper {
	typedef _Data_Traits						data_traits;
	typedef typename data_traits::vector_type	vector_type;
	typedef typename data_traits::point_type	point_type;

	rhs_wrapper(const data_traits& traits) : __traits(traits) {}

	bool operator()(double t, const point_type& x, vector_type& f) const {
		return __traits.get_value(x, f);
	}

	const data_traits& __traits;
};

}



#endif








