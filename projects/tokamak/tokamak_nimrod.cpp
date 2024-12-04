#include <teem/nrrd.h>

#include <stdexcept>
#include <iostream>
#include <cassert>

#include "tokamak_nimrod.hpp"

using namespace nvis;

#include <hdf5.h>

// #define TWO_ARGS

bool
H5NIMROD_read_dims(hid_t parent_id,
                   const char *dataset_name, int *ndims, hsize_t * grid_dims)
{
	hsize_t dims[16];
	int i, j;

#ifdef TWO_ARGS
	hid_t dataset_id = H5Dopen(parent_id, dataset_name);
#else
	hid_t dataset_id = H5Dopen(parent_id, dataset_name, NULL);
#endif
	if (dataset_id < 0)
		return false;
	hid_t dataspace_id = H5Dget_space(dataset_id);
	if (dataspace_id < 0)
		return false;

	*ndims = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

	for (i = 0, j = *ndims - 1; i < *ndims; i++, j--)
		grid_dims[i] = dims[j];

	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
	return true;
}

bool H5NIMROD_read_float32_array(hid_t parent_id,
                                 const char *dataset_name,
                                 hsize_t * offset,
                                 int ndims, hsize_t * dims, float *array)
{
	hid_t dataspace, dataset, memspace;
#ifdef TWO_ARGS
	dataset = H5Dopen(parent_id, dataset_name);
#else
	dataset = H5Dopen(parent_id, dataset_name, NULL);
#endif
	if (dataset < 0) {
		printf("could not open dataset %s\n", dataset_name);
		return false;
	}
	if (offset == NULL)
		dataspace = H5S_ALL;
	else
		dataspace = H5Dget_space(dataset);	/* dataspace identifier */

	if (dims == NULL) {
		memspace = H5S_ALL;
	}
	else {
		memspace = H5Screate_simple(ndims, dims, NULL);
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
		                    dims, NULL);
	}

	/*  hid_t datatype=H5Dget_type(dataset); */
	H5Dread(dataset,		/* handle for the dataset */
	        H5T_NATIVE_FLOAT,	/* the datatype we use in memory
				   you can change it to FLOAT if you want */
	        memspace,		/* shape/size of data
				   in memory (the complement to disk hyperslab) */
	        dataspace,		/* shape/size of data on disk  i
				   (get hyperslab if needed) */
	        H5P_DEFAULT,		/* ignore... its for parallel reads */
	        array);		/* the data array we are reading into */

	if (memspace != H5S_ALL)
		H5Sclose(memspace);
	if (dataspace != H5S_ALL)
		H5Sclose(dataspace);

	H5Dclose(dataset);		/* release the dataset handle */
	return true;
}

// ---------------------------------------------------------

tokamak_nimrod::tokamak_nimrod(const std::string& h5file,
                               const std::string& step,
                               bool _square) :
		square(_square)
{
	// load geometry
	hid_t file_id, root_id, group_id;
	char* string_attrib;

	file_id = H5Fopen(h5file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	int ndims;
	hsize_t dims[4];
	H5NIMROD_read_dims(file_id, "cartGrid", &ndims, dims);

	if (ndims != 4)
		throw std::runtime_error("NIMROD: can only handle three-dimensional data");

	std::cout << "grid dimensions: " << dims[0] << 'x' << dims[1] << 'x' << dims[2] << 'x' << dims[3] << '\n';

	hsize_t npoints = dims[1] * dims[2] * dims[3];

	fvec3* points = new fvec3[npoints];
	H5NIMROD_read_float32_array(file_id, "cartGrid", NULL, ndims, NULL, (float*)points);

	char stepname[256];
	snprintf(stepname, 256, "step_%07u/B", atoi(step.c_str()));

	data = new fvec3[npoints];
	H5NIMROD_read_float32_array(file_id, stepname, NULL, ndims, NULL, (float*)data);

	// compute the poloidal slice

	Nx = dims[1];
	Ny = dims[2];
	Nz = dims[3];

	// compute poloidal slice (z=0)

	slice.resize(Nx*Ny);

	for (unsigned int i = 0; i < Ny*Nx; ++i)
		slice[i] = poloidal_project(points[i]);

	center = slice[0];

	// compute outer angles and fit radius spline

	ang.resize(Ny);
	rad.resize(Ny);

	for (unsigned int i = 0, j = Nx - 1; i < Ny; ++i, j += Nx) {
		rad[i] = poloidal_radius(slice[j]);
		ang[i] = poloidal_angle(slice[j]);
	}

	ang[0] = ang[Ny-1] - 2 * M_PI;
	rspline = make_catmull_rom_spline(std::vector<double>(ang.begin(), ang.end()), rad);

	// compute slice cell bounding boxes
	cellbb.resize((Nx - 1)*(Ny - 1));

	const double eps = sqrt(std::numeric_limits<double>::epsilon());

	for (unsigned int y = 0; y < Ny - 1; ++y) {
		for (unsigned int x = 0; x < Nx - 1; ++x) {
			bbox2& bb = cellbb[y*(Nx-1)+x];
			int b = y * Nx + x;

			bb.add(slice[b]);
			bb.add(slice[b+1]);
			bb.add(slice[b+Nx]);
			bb.add(slice[b+Nx+1]);

			bb.min() -= vec2(eps);
			bb.max() += vec2(eps);
		}
	}

	delete[] points;

	// compute bounding box
	if (square)
		bbox = bbox2(vec2(0, 0), vec2(1, 1));
	else
		bbox = bbox2(slice.begin(), slice.end());
}

// --------------------------------------------------------------------------

vec2 tokamak_nimrod::project(const vec3& p) const
{
	vec2 pp = poloidal_project(p);

	if (square) {
		double r = poloidal_radius(pp);
		double a = poloidal_angle(pp);

		pp[0] = a / (2 * M_PI);
		pp[1] = r / rspline(a);
	}

	return pp;
}

mat2 tokamak_nimrod::project(const mat3& p) const
{
	std::cerr << "tokamak_nimrod::project is not implemented for jacobian matrices\n";
	throw std::exception();

}

vec3 tokamak_nimrod::unproject(const vec2& p) const
{
	vec2 pp = p;

	if (square) {
		double a = 2 * M_PI * p[0];
		double r = p[1] * rspline(a);

		pp = center + r * vec2(cos(a), sin(a));
	}

	return poloidal_unproject(pp);
}

double tokamak_nimrod::poloidal_angle(const vec2& p) const
{
	double dy = center[1] - p[1];
	double dx = center[0] - p[0];

	return atan2(dy, dx) + M_PI;
}

double tokamak_nimrod::poloidal_radius(const vec2& p) const
{
	return norm(center - p);
}

vec2 tokamak_nimrod::poloidal_project(const nvis::vec3& v) const
{
	return nvis::vec2(sqrt(v[0]*v[0] + v[1]*v[1]), v[2]);
}

vec3 tokamak_nimrod::poloidal_unproject(const nvis::vec2& v) const
{
	return nvis::vec3(v[0], 0.0, v[1]);
}

// --------------------------------------------------------------------------

bool tokamak_nimrod::intersect(const nvis::vec3& p0, const nvis::vec3& p1) const
{
	double f0 = p0[1];
	double f1 = p1[1];

	if (p0[0] < 0)
		return false;

	return f0 < f1 ? (f0 < 0 && f1 > 0) : (f0 > 0 && f1 < 0);
}

// --------------------------------------------------------------------------

bool tokamak_nimrod::interpolate(const vec3& p, vec3& r) const
{
	const double lepsilon = 0.0 - sqrt(std::numeric_limits<double>::epsilon());
	const double uepsilon = 1.0 + sqrt(std::numeric_limits<double>::epsilon());

	const vec2 pp = poloidal_project(p);

	// determine y
	double y = poloidal_angle(pp);
	int yi = std::distance(ang.begin(), std::lower_bound(ang.begin(), ang.end(), y));

	if (yi >= Ny)
		yi = 0;
	else if (yi > 0)
		--yi;

	// std::cout << "  pp = " << pp << '\n';
	// std::cout << "  y = " << y << ' ' << ang[yi] << ' ' << ang[yi+1] << ", ";
	// std::cout << "  i = " << yi << '\n';

	for (unsigned int i = 0, c_id = yi * (Nx - 1); i < Nx - 1; ++i, ++c_id) {
		if (!cellbb[c_id].inside(pp))
			continue;

		int base = c_id % (Nx - 1) + c_id / (Nx - 1) * Nx;

		const vec2& pos0 = slice[base];
		const vec2& pos1 = slice[base+1];
		const vec2& pos2 = slice[base+Nx+1];
		const vec2& pos3 = slice[base+Nx];

		vec2 p1 = pos1 - pos0;
		vec2 p2 = pos2 - pos0;
		vec2 p3 = pos3 - pos0;
		vec2 p4 = p2 - p1 - p3;

		double det1 = cross(p3, p4);
		double det2 = cross(p3, p1);
		double det3 = cross(p1, p4);

		vec2 rp = pp - pos0;

		double a = det1;
		double b = det2 + cross(p4, rp);
		double c = cross(p1, rp);

		// std::cout << "OPD: " << det1 << ' ' << det2 << ' ' << det3 << '\n';
		// std::cout << "OPL: " << a << ' ' << b << ' ' << c << '\n';
		// solve quadratic equation

		double v = -1.0;

		if (std::abs(a) < 1e-9) {
			// std::cout << "linear case\n";

			if (std::abs(b) > 1e-9)
				v = -c / b;
		}
		else {
			double d = b * b - 4 * a * c;

			if (d > 0) {
				double q = (b >= 0.0 ? -0.5 * (b + sqrt(d)) : -0.5 * (b - sqrt(d)));

				v = q / a;

				if (v < lepsilon || v > uepsilon)
					v = c / q;
			}
			else if (d == 0.0) {
				if (a != 0.0)
					v = -0.5 * b / a;
				else {
					std::cout << "degenerate case\n";
					continue;
				}
			}
			else {
				std::cout << "imaginary solutions\n";
				continue;
			}
		}

		if (v < lepsilon || v > uepsilon) {
			// std::cout << "not inside\n";
			continue;
		}

		vec2 d = p1 + v * p4;

		double u = std::abs(d[0]) > std::abs(d[1]) ?
		           (rp[0] - v * p3[0]) / d[0] :
		           (rp[1] - v * p3[1]) / d[1];

		if (u < lepsilon || u > uepsilon) {
			// std::cout << "not inside\n";
			continue;
		}

		// std::cout << "inside, uv = " << u << ' ' << v << '\n';

		double z = 0.5 * (Nz - 1) * (atan2(p[1], -p[0]) / M_PI + 1);
		int zi = floor(z);

		// std::cout << "z = " << z << ", zi = " << zi << '\n';

		vec3 r0, r1;

		// std::cout << ((zi+0)%Nz)*(Nx*Ny) + base << '\n'
		//           << ((zi+0)%Nz)*(Nx*Ny) + base+1 << '\n'
		//           << ((zi+0)%Nz)*(Nx*Ny) + base+Nx+1 << '\n'
		//           << ((zi+0)%Nz)*(Nx*Ny) + base+Nx << '\n';
		// std::cout << ((zi+1)%Nz)*(Nx*Ny) + base << '\n'
		//           << ((zi+1)%Nz)*(Nx*Ny) + base+1 << '\n'
		//           << ((zi+1)%Nz)*(Nx*Ny) + base+Nx+1 << '\n'
		//           << ((zi+1)%Nz)*(Nx*Ny) + base+Nx << '\n';
		//
		{
			const fvec3* slice = data + ((zi + 0) % Nz) * Nx * Ny;

			fvec3 t0 = slice[base];
			fvec3 t1 = slice[base+1];
			fvec3 t2 = slice[base+Nx+1];
			fvec3 t3 = slice[base+Nx];

			r0 = (1 - v) * ((1 - u) * t0 + u * t1) + v * ((1 - u) * t3 + u * t2);
		}

		{
			const fvec3* slice = data + ((zi + 1) % Nz) * Nx * Ny;

			fvec3 t0 = slice[base];
			fvec3 t1 = slice[base+1];
			fvec3 t2 = slice[base+Nx+1];
			fvec3 t3 = slice[base+Nx];

			r1 = (1 - v) * ((1 - u) * t0 + u * t1) + v * ((1 - u) * t3 + u * t2);
		}

		r = (z - zi) * r1 + (zi + 1 - z) * r0;
		return true;
	}

	return false;
}

vec3 tokamak_nimrod::operator()(const double&,
                                const vec3& p) const
{
	vec3 r;

	if (!interpolate(p, r))
		throw undefined_point();

	return r;
}

mat3 tokamak_nimrod::derivative(const double&,
                                const vec3& p) const
{
	throw undefined_point();
	return mat3();
}







