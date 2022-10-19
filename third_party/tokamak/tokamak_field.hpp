#ifndef __tokamak_field_hpp
#define __tokamak_field_hpp

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>

#include <stdexcept>

namespace nvis
{
typedef fixed_vector<double, 9> 	vec9;
typedef fixed_vector<float, 9>		fvec9;
}

class tokamak_field
{
public:

	struct undefined_point : public std::runtime_error {
		undefined_point() : std::runtime_error("undefined point") {
		}
	};

	virtual nvis::vec3 operator()(const double& t, const nvis::vec3& y) const = 0;
	virtual nvis::mat3 derivative(const double& t, const nvis::vec3& y) const = 0;

	virtual nvis::vec4 plane() const = 0;
	virtual nvis::vec2 project(const nvis::vec3& pos) const = 0;
	virtual nvis::mat2 project(const nvis::mat3& pos) const = 0;
	virtual nvis::vec3 unproject(const nvis::vec2& pos) const = 0;

	virtual bool intersect(const nvis::vec3& p0, const nvis::vec3& p1) const {
		return false;
	};

	virtual tokamak_field* clone() const = 0;

	virtual nvis::bbox2 bounds() const = 0;

	virtual void periodic_coordinates(bool periodic = false) const {}
};

#endif // __tokamak_base_hpp

