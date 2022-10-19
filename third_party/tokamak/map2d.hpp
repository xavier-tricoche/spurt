#ifndef __map2d_hpp
#define __map2d_hpp

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <vector>

namespace nvis
{
typedef fixed_vector< bool, 2 > bvec2;

struct return_map_info {
	nvis::vec2 x, err;
	nvis::mat2 J;
	double t;
	unsigned int nsteps;
	unsigned int total_nsteps;
};

}

struct map2d {
	
	typedef nvis::return_map_info	value_type;
	
	struct map_undefined {
	};

	virtual nvis::vec2 map(const nvis::vec2& in, int niter) const = 0;
	virtual void map(const nvis::vec2& in, std::vector<nvis::vec2>& out, int niter) const = 0;
	
	virtual value_type map_complete(const nvis::vec2& in, int niter) const = 0;
	virtual void map_complete(const nvis::vec2& in, std::vector<value_type>& out, int niter) const = 0;

	virtual nvis::bbox2 bounds() const = 0;
	virtual nvis::bvec2 periodic() const = 0;

	virtual map2d* clone() const = 0;

	static map2d* load(const std::string& desc);
};

#endif // __map2d_hpp




