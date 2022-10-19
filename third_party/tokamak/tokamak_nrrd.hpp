#ifndef __tokamak_nimrod_nrrd_hpp
#define __tokamak_nimrod_nrrd_hpp

#include <cmath>
#include <vector>
#include <teem/nrrd.h>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>

#include <tokamak/tokamak_field.hpp>

namespace nvis
{
typedef fixed_matrix<float, 3> 	fmat3;
}

class tokamak_nrrd: public tokamak_field
{
public:

	tokamak_nrrd(const std::string& file);

	nvis::vec3 operator()(const double& t, const nvis::vec3& y) const;
	nvis::mat3 derivative(const double& t, const nvis::vec3& y) const;

	nvis::vec2 project(const nvis::vec3& v) const;
	nvis::mat2 project(const nvis::mat3& v) const;
	nvis::vec3 unproject(const nvis::vec2& v) const;

	nvis::vec4 plane() const {
		return nvis::vec4(0, -1, 0, 0);
	}

	tokamak_nrrd* clone() const {
		return new tokamak_nrrd(*this);
	}

	nvis::bbox2 bounds() const {
		return nvis::bbox2(nvis::vec2(0, 0), nvis::vec2(size[1] - 1, size[0] - 1));
	}

	void get_slice(std::vector<nvis::vec2>& parm, std::vector<nvis::vec2>& phys,
	               unsigned int& nx, unsigned int& ny) const;

	// ---

	bool intersect(const nvis::vec3& p0, const nvis::vec3& p1) const;

	static void load(const std::string& file,
	                 nvis::fvec3*& points, nvis::fvec3*& data, nvis::uvec3& size);

	void periodic_coordinates(bool periodic = true) const {
		periodic_coords = periodic;
	}

	const nvis::fvec3* get_data(unsigned int& nx, unsigned int& ny, unsigned int& nz) const {
		nx = size[0];
		ny = size[1];
		nz = size[2];
		return data;
	}

private:

	nvis::vec2 poloidal_project(const nvis::vec3& v) const;
	nvis::mat2 poloidal_project(const nvis::mat3& v) const;
	nvis::vec3 poloidal_unproject(const nvis::vec2& v) const;

	double poloidal_angle(const nvis::vec2&) const;
	double poloidal_radius(const nvis::vec2&) const;

	bool interpolate(const nvis::vec3& p, nvis::vec3& r, const nvis::fvec3* variable = 0) const;
	bool jacobian(const nvis::vec3& p, nvis::mat3& r, const nvis::fmat3* variable = 0) const;

	void precompute_jacobian() ;

	nvis::fvec3*             data;
	nvis::fmat3*			 jdata;
	nvis::fvec3*             points;
	nvis::uvec3              size;
	mutable bool			 periodic_coords;
};

#endif // __tokamak_nrrd_hpp



