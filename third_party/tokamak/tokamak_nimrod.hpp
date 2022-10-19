#ifndef __tokamak_nimrod_hpp
#define __tokamak_nimrod_hpp

#include <cmath>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/bezier_spline.hpp>

#include "tokamak_field.hpp"

class tokamak_nimrod: public tokamak_field
{
public:

	tokamak_nimrod(const std::string& h5file,
	               const std::string& step,
	               bool square = false);

	// tokamak_nimrod( const tokamak_nimrod& other );

	// ~tokamak_nimrod();

	nvis::vec3 operator()(const double& t, const nvis::vec3& y) const;
	nvis::mat3 derivative(const double& t, const nvis::vec3& y) const;

	nvis::vec2 project(const nvis::vec3& v) const;
	nvis::mat2 project(const nvis::mat3& pos) const;
	nvis::vec3 unproject(const nvis::vec2& v) const;

	nvis::vec4 plane() const {
		return nvis::vec4(0, -1, 0, 0);
	}

	tokamak_nimrod* clone() const {
		return new tokamak_nimrod(*this);
	}

	nvis::bbox2 bounds() const {
		return bbox;
	}

	// ---

	bool intersect(const nvis::vec3& p0, const nvis::vec3& p1) const;

private:

	nvis::vec2 poloidal_project(const nvis::vec3& v) const;
	nvis::vec3 poloidal_unproject(const nvis::vec2& v) const;

	double poloidal_angle(const nvis::vec2&) const;
	double poloidal_radius(const nvis::vec2&) const;

	bool interpolate(const nvis::vec3& p, nvis::vec3& r) const;

	nvis::bezier_spline<double, 3> rspline;

	nvis::fvec3*             data;
	std::vector<double>      ang;
	std::vector<double>      rad;
	std::vector<nvis::fvec2> slice;
	std::vector<nvis::bbox2> cellbb;

	nvis::vec2               center;
	nvis::bbox2              bbox;

	bool         square;
	unsigned int Nx, Ny, Nz;
};

#endif // __tokamak_nimrod_hpp


