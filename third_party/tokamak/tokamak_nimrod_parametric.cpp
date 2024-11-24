#include <teem/gage.h>

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <fstream>

#include <nvis-math/fixed_matrix.hpp>
#include "tokamak_nimrod_parametric.hpp"

using namespace nvis;

#include <hdf5.h>

#define USE_TRILINEAR
#define THETA_NORMALIZED

// #define TWO_ARGS

// ---------------------------------------------------------

#define _L(x,y,z) ((x)+(y)*size[0]+(z)*(size[0]*size[1]))

typedef float gage_t;

static inline double det3(const vec3& v0, const vec3& v1, const vec3& v2)
{
	return v0[0]*(v1[1]*v2[2] - v2[1]*v1[2]) +
	       v1[0]*(v2[1]*v0[2] - v0[1]*v2[2]) +
	       v2[0]*(v0[1]*v1[2] - v1[1]*v0[2]);
}

tokamak_nimrod_parametric::tokamak_nimrod_parametric(const std::string& h5file,
        const std::string& step) : periodic_coords(true)
{
	load(h5file, step, points, data, size);
	precompute_jacobian();
}

void tokamak_nimrod_parametric::precompute_jacobian()
{
#ifndef USE_TRILINEAR
	const int offset = 10;

	// create block of floats with periodic boundaries in y and z
	int nx = size[0], ny = size[1], nz = size[2];
	int Nx = size[0], Ny = size[1] + 2 * offset, Nz = size[2] + 2 * offset;
	float *tmp_data = (float*)calloc(3 * Nx * Ny * Nz, sizeof(float));
	for (int k = 0 ; k < nz ; ++k) {
		for (int j = 0 ; j < ny ; ++j) {
			for (int i = 0 ; i < nx ; ++i) {
				// add offset in y and z
				int j__ = j + offset;
				int k__ = k + offset;
				const nvis::fvec3& v = data[_L(i, j, k)];
				int id = i + Nx * (j__ + Ny * k__);
				tmp_data[3*id+0] = v[0];
				tmp_data[3*id+1] = v[1];
				tmp_data[3*id+2] = v[2];
				if (j < offset) {
					int _id = i + Nx * ((j__ + ny) + Ny * k__);
					tmp_data[3*_id+0] = v[0];
					tmp_data[3*_id+1] = v[1];
					tmp_data[3*_id+2] = v[2];
				}
				else if (j > ny - 1 - offset) {
					int _id = i + Nx * ((j__ - ny) + Ny * k__);
					tmp_data[3*_id+0] = v[0];
					tmp_data[3*_id+1] = v[1];
					tmp_data[3*_id+2] = v[2];
				}
				if (k < offset) {
					int _id = i + Nx * (j__ + Ny * (k__ + nz));
					tmp_data[3*_id+0] = v[0];
					tmp_data[3*_id+1] = v[1];
					tmp_data[3*_id+2] = v[2];
				}
				else if (k > nz - 1 - offset) {
					int _id = i + Nx * (j__ + Ny * (k__ - nz));
					tmp_data[3*_id+0] = v[0];
					tmp_data[3*_id+1] = v[1];
					tmp_data[3*_id+2] = v[2];
				}
			}
		}
	}

	// wrap it in a Nrrd
	Nrrd *nrrd = nrrdNew();
	if (nrrdWrap_va(nrrd, tmp_data, nrrdTypeFloat, 4, 3, Nx, Ny, Nz)) {
		std::cerr << biffGetDone(NRRD) << std::endl;
		exit(-1);
	}
	nrrdAxisInfoSet_va(nrrd, nrrdAxisInfoSpacing, AIR_NAN, 1.0, 1.0, 1.0);
	nrrdAxisInfoSet_va(nrrd, nrrdAxisInfoCenter, nrrdCenterUnknown,
	                   nrrdCenterNode, nrrdCenterNode, nrrdCenterNode);

	// use gage to measure Jacobian matrix
	gageContext *ctx;
	gagePerVolume *pv;
	const double *_J;
	double kparm[NRRD_KERNEL_PARMS_NUM];

	ctx = gageContextNew();
	kparm[0] = 1.0;
	kparm[1] = 1.0;
	kparm[2] = 0.0;

	int E = 0;
	if (!E) E |= !(pv = gagePerVolumeNew(ctx, nrrd, gageKindVec));
	if (!E) E |= gagePerVolumeAttach(ctx, pv);

	// BC Cubic
	if (!E) E |= gageKernelSet(ctx, gageKernel00, nrrdKernelBCCubic, kparm);
	if (!E) E |= gageKernelSet(ctx, gageKernel11, nrrdKernelBCCubicD, kparm);
	if (!E) E |= gageKernelSet(ctx, gageKernel22, nrrdKernelBCCubicDD, kparm);

	if (!E) E |= gageQueryItemOn(ctx, pv, gageVecJacobian);
	if (!E) E |= gageUpdate(ctx);
	if (E) {
		std::cout << "tokamak_nimrod_parametric::compute_jacobian:" << std::endl
		          << biffGetDone(GAGE)
		          << std::endl;
		return;
	}
	_J = gageAnswerPointer(ctx, pv, gageVecJacobian);

	jdata = new nvis::fmat3[nx*ny*nz];
	for (unsigned int k = 0 ; k < nz ; ++k) {
		for (unsigned int j = 0 ; j < ny ; ++j) {
			for (unsigned int i = 0 ; i < nx ; ++i) {
				if (gageProbe(ctx, (double)i, (double)(j + offset), (double)(k + offset))) {
					std::cerr << "invalid position at " << nvis::vec3(i, j, k) << '\n';
					continue;
				}
				for (unsigned int r = 0 ; r < 3 ; ++r) {
					for (unsigned int c = 0 ; c < 3 ; ++c) {
						jdata[i+nx*(j+ny*k)][r][c] = _J[3*r+c];
					}
				}
			}
		}
	}
	nrrdNuke(nrrd);

	// Nrrd *nout = nrrdNew();
	// float *blah = (float*)calloc(9 * nx * ny * nz, sizeof(float));
	// for (unsigned int i = 0 ; i < nx*ny*nz ; ++i) {
	// 	for (unsigned int j = 0 ; j < 9 ; ++j) blah[9*i+j] = jdata[i][j/3][j%3];
	// }
	// if (nrrdWrap_va(nout, blah, nrrdTypeFloat, 4, 9, nx, ny, nz)) {
	// 	std::cerr << biffGetDone(NRRD) << std::endl;
	// 	exit(-1);
	// }
	// nrrdAxisInfoSet_va(nout, nrrdAxisInfoSpacing, AIR_NAN, 1.0, 1.0, 1.0);
	// nrrdAxisInfoSet_va(nout, nrrdAxisInfoCenter, nrrdCenterUnknown,
	//                    nrrdCenterNode, nrrdCenterNode, nrrdCenterNode);
	// if (nrrdSave("test.nrrd", nout, NULL)) {
	// 	std::cerr << "precompute_jacobian: " << biffGetDone(NRRD) << std::endl;
	// 	exit(-1);
	// }
	// std::cout << "saved test.nrrd\n";

#endif
}

void tokamak_nimrod_parametric::load(const std::string& h5file, const std::string& step,
                                     fvec3*& points, fvec3*& data, uvec3& size)
{
	hid_t file_id = H5Fopen(h5file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

#ifdef TWO_ARGS
	hid_t dataset_id = H5Dopen(file_id, "cartGrid");
#else
	hid_t dataset_id = H5Dopen(file_id, "cartGrid", NULL);
#endif

	if (dataset_id < 0)
		throw std::runtime_error("couldn't get \"cartGrid\" dataset");

	hid_t dataspace_id = H5Dget_space(dataset_id);

	if (dataspace_id < 0)
		throw std::runtime_error("couldn't get \"cartGrid\" dataspace");


	int ndims = H5Sget_simple_extent_ndims(dataspace_id);

	if (ndims != 4)
		throw std::runtime_error("dataspace \"cartGrid\" is not of dimension 4");

	hsize_t dims[4];
	H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
	H5Sclose(dataspace_id);

	std::reverse(dims, dims + 4);
	std::copy(dims + 1, dims + 4, size.begin());
	std::reverse(dims, dims + 4);

	std::cout << "grid dimensions: " << size[0] << 'x' << size[1] << 'x' << size[2] << '\n';
	hsize_t npoints = dims[0] * dims[1] * dims[2];

	std::cout << "npoints = " << npoints << '\n';

	points = new fvec3[npoints];
	H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, points);
	H5Dclose(dataset_id);

	// ---

	char stepname[256];
	snprintf(stepname, 256, "step_%07u/B", atoi(step.c_str()));

#ifdef TWO_ARGS
	dataset_id = H5Dopen(file_id, stepname);
#else
	dataset_id = H5Dopen(file_id, stepname, NULL);
#endif

	if (dataset_id < 0)
		throw std::runtime_error("couldn't get step dataset");

	data = new fvec3[npoints];
	H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(dataset_id);
	H5Fclose(file_id);

	// compute the poloidal slice

	float maxnorm = 0.0;

	for (unsigned int z = 0; z < size[2]; ++z) {
		for (unsigned int y = 0; y < size[1]; ++y) {
			for (unsigned int x = 0; x < size[0]; ++x) {
				mat3 J;

				if (x == 0)
					J[0] = points[_L(1, y, z)] - points[_L(0, y, z)];
				else if (x == size[0] - 1)
					J[0] = points[_L(size[0] - 1, y, z)] - points[_L(size[0] - 2, y, z)];
				else
					J[0] = 0.5 * (points[_L(x+1, y, z)] - points[_L(x-1, y, z)]);

				if (y == 0)
					J[1] = 0.5 * (points[_L(x, 1, z)] - points[_L(x, size[1] - 2, z)]);
				else if (y == size[1] - 1)
					J[1] = 0.5 * (points[_L(x, 1, z)] - points[_L(x, size[1] - 2, z)]);
				else
					J[1] = 0.5 * (points[_L(x, y+1, z)] - points[_L(x, y-1, z)]);

				if (z == 0)
					J[2] = 0.5 * (points[_L(x, y, 1)] - points[_L(x, y, size[2] - 2)]);
				else if (z == size[2] - 1)
					J[2] = 0.5 * (points[_L(x, y, 1)] - points[_L(x, y, size[2] - 2)]);
				else
					J[2] = 0.5 * (points[_L(x, y, z+1)] - points[_L(x, y, z-1)]);

				vec3 r = data[_L(x, y, z)];

				float det = det3(J[0], J[1], J[2]);

				if (fabs(det) < 1e-20)
					data[_L(x, y, z)] = 0;
				else {
					data[_L(x, y, z)][0] = det3(r, J[1], J[2]) / det;
					data[_L(x, y, z)][1] = det3(J[0], r, J[2]) / det;
					data[_L(x, y, z)][2] = det3(J[0], J[1], r) / det;
				}

				maxnorm = std::max(norm(data[_L(x, y, z)]), (double)maxnorm);
			}
		}
	}

	std::cout << maxnorm << '\n';

	//
	// {
	//     std::ofstream out( "points" );
	//
	//     for( unsigned int i=size[0]-1; i<size[0]*size[1]; i+=size[0] )
	//         out << points[i][0] << ' ' << points[i][1] << ' ' << points[i][2] << '\n';
	//
	//     out.close();
	// }
}

// --------------------------------------------------------------------------


double mod(double a, double b)
{
	return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}

vec2 tokamak_nimrod_parametric::project(const vec3& p) const
{
	return poloidal_project(p);
}

vec3 tokamak_nimrod_parametric::unproject(const vec2& p) const
{
	return poloidal_unproject(p);
}

vec2 tokamak_nimrod_parametric::poloidal_project(const nvis::vec3& v) const
{
	return (periodic_coords ?
	        nvis::vec2(mod(v[1], size[1] - 1),
	                   mod(v[0], size[0] - 1)) :
	        nvis::vec2(v[1], v[0]));
}

vec3 tokamak_nimrod_parametric::poloidal_unproject(const nvis::vec2& v) const
{
	return nvis::vec3(v[1], v[0], 0.0);
}

mat2 tokamak_nimrod_parametric::project(const mat3& m) const
{
	return poloidal_project(m);
}

mat2 tokamak_nimrod_parametric::poloidal_project(const nvis::mat3& m) const
{
	mat2 r;
	r[0][0] = m[1][1];
	r[0][1] = m[1][0];
	r[1][0] = m[0][1];
	r[1][1] = m[0][0];

	return r;
}

void tokamak_nimrod_parametric::get_slice(std::vector<nvis::vec2>& parm,
        std::vector<nvis::vec2>& phys,
        unsigned int& n1, unsigned int& n2) const
{
	n1 = size[0];
	n2 = size[1];

	phys.resize(n1*n2);
	parm.resize(n1*n2);

	for (unsigned int i = 0; i < n1*n2; ++i) {
		phys[i][0] = points[i][0];
		phys[i][1] = points[i][2];

		parm[i][0] = (double)(i / n1) / n2;
		parm[i][1] = (double)(i % n1) / n1;
	}
}

// --------------------------------------------------------------------------

bool tokamak_nimrod_parametric::intersect(const nvis::vec3& p0, const nvis::vec3& p1) const
{
	int z0 = (int)trunc(p0[2] / (size[2] - 1));
	int z1 = (int)trunc(p1[2] / (size[2] - 1));

	return z0 != z1;
}

// --------------------------------------------------------------------------

bool tokamak_nimrod_parametric::interpolate(const vec3& p, vec3& r, const fvec3* variable) const
{
	if (p[0] > size[0] - 1 || p[0] < 0)
		return false;

	int cx = (int)floor(p[0]);
	int cy = (int)floor(mod(p[1], size[1] - 1));
	int cz = (int)floor(mod(p[2], size[2] - 1));

	const float c[3] = { static_cast<float>(p[0] - floor(p[0])), 
                         static_cast<float>(p[1] - floor(p[1])), 
                         static_cast<float>(p[2] - floor(p[2])) };

	unsigned int base = cx + size[0] * cy + cz * size[0] * size[1];

	unsigned int indices[8] = {
		base, base + 1, base + size[0] + 1, base + size[0],
		base + size[0]*size[1], base + 1 + size[0]*size[1], base + size[0] + 1 + size[0]*size[1], base + size[0] + size[0]*size[1],
	};

	const float d[3] = { 1.0f - c[0], 1.0f - c[1], 1.0f - c[2] };

	if (variable == 0)
		variable = data;

	r = variable[indices[0]] * d[0] * d[1] * d[2] +
	    variable[indices[1]] * c[0] * d[1] * d[2] +
	    variable[indices[2]] * c[0] * c[1] * d[2] +
	    variable[indices[3]] * d[0] * c[1] * d[2] +
	    variable[indices[4]] * d[0] * d[1] * c[2] +
	    variable[indices[5]] * c[0] * d[1] * c[2] +
	    variable[indices[6]] * c[0] * c[1] * c[2] +
	    variable[indices[7]] * d[0] * c[1] * c[2];

#ifdef THETA_NORMALIZED
	r /= r[2];
#endif

	return true;
}

// --------------------------------------------------------------------------

bool tokamak_nimrod_parametric::jacobian(const vec3& p, mat3& r, const fmat3* variable) const
{
	if (p[0] > size[0] - 1 || p[0] < 0)
		return false;

	int cx = (int)floor(p[0]);
	int cy = (int)floor(mod(p[1], size[1] - 1));
	int cz = (int)floor(mod(p[2], size[2] - 1));

	const float c[3] = { static_cast<float>(p[0] - floor(p[0])), 
                         static_cast<float>(p[1] - floor(p[1])), 
                         static_cast<float>(p[2] - floor(p[2])) };

	unsigned int base = cx + size[0] * cy + cz * size[0] * size[1];

	unsigned int indices[8] = {
		base, base + 1, base + size[0] + 1, base + size[0],
		base + size[0]*size[1], base + 1 + size[0]*size[1], base + size[0] + 1 + size[0]*size[1], base + size[0] + size[0]*size[1],
	};

	const float d[3] = { 1.0f - c[0], 1.0f - c[1], 1.0f - c[2] };

#ifndef USE_TRILINEAR
	if (variable == 0)
		variable = jdata;

	r = (d[0] * d[1] * d[2]) * variable[indices[0]] +
	    (c[0] * d[1] * d[2]) * variable[indices[1]] +
	    (c[0] * c[1] * d[2]) * variable[indices[2]] +
	    (d[0] * c[1] * d[2]) * variable[indices[3]] +
	    (d[0] * d[1] * c[2]) * variable[indices[4]] +
	    (c[0] * d[1] * c[2]) * variable[indices[5]] +
	    (c[0] * c[1] * c[2]) * variable[indices[6]] +
	    (d[0] * c[1] * c[2]) * variable[indices[7]];
#else
	const fvec3* _variable = data;

	nvis::fvec3 ddx = (-1 * d[1] * d[2]) * _variable[indices[0]] +
	                  (1 * d[1] * d[2]) * _variable[indices[1]] +
	                  (1 * c[1] * d[2]) * _variable[indices[2]] +
	                  (-1 * c[1] * d[2]) * _variable[indices[3]] +
	                  (-1 * d[1] * c[2]) * _variable[indices[4]] +
	                  (1 * d[1] * c[2]) * _variable[indices[5]] +
	                  (1 * c[1] * c[2]) * _variable[indices[6]] +
	                  (-1 * c[1] * c[2]) * _variable[indices[7]];

	nvis::fvec3 ddy = (d[0] * -1 * d[2]) * _variable[indices[0]] +
	                  (c[0] * -1 * d[2]) * _variable[indices[1]] +
	                  (c[0] * 1 * d[2]) * _variable[indices[2]] +
	                  (d[0] * 1 * d[2]) * _variable[indices[3]] +
	                  (d[0] * -1 * c[2]) * _variable[indices[4]] +
	                  (c[0] * -1 * c[2]) * _variable[indices[5]] +
	                  (c[0] * 1 * c[2]) * _variable[indices[6]] +
	                  (d[0] * 1 * c[2]) * _variable[indices[7]];

	nvis::fvec3 ddz = (d[0] * d[1] * -1) * _variable[indices[0]] +
	                  (c[0] * d[1] * -1) * _variable[indices[1]] +
	                  (c[0] * c[1] * -1) * _variable[indices[2]] +
	                  (d[0] * c[1] * -1) * _variable[indices[3]] +
	                  (d[0] * d[1] * 1) * _variable[indices[4]] +
	                  (c[0] * d[1] * 1) * _variable[indices[5]] +
	                  (c[0] * c[1] * 1) * _variable[indices[6]] +
	                  (d[0] * c[1] * 1) * _variable[indices[7]];

	for (int i = 0 ; i < 3 ; ++i) {
		r[i][0] = ddx[i];
		r[i][1] = ddy[i];
		r[i][2] = ddz[i];
	}
#endif

#ifdef THETA_NORMALIZED
	nvis::vec3 f = (d[0] * d[1] * d[2]) * data[indices[0]] +
	               (c[0] * d[1] * d[2]) * data[indices[1]] +
	               (c[0] * c[1] * d[2]) * data[indices[2]] +
	               (d[0] * c[1] * d[2]) * data[indices[3]] +
	               (d[0] * d[1] * c[2]) * data[indices[4]] +
	               (c[0] * d[1] * c[2]) * data[indices[5]] +
	               (c[0] * c[1] * c[2]) * data[indices[6]] +
	               (d[0] * c[1] * c[2]) * data[indices[7]];

	nvis::mat3 J(r);
	J /= f[2];
	f /= f[2];
	r = J - nvis::outer(f, J[2]);
	// std::cout << "theta-normalized velocity(" << p << ") = " << f << '\n';
	// std::cout << "theta-normalized jacobian(" << p << ") = " << r << '\n';
#endif

	return true;
}

// --------------------------------------------------------------------------

vec3 tokamak_nimrod_parametric::operator()(const double&,
        const vec3& p) const
{
	vec3 r;

	if (!interpolate(p, r)) {		//
		// std::cerr << "field at " << p << " is undefined\n";
		throw undefined_point();
	}

	return r;
}

// --------------------------------------------------------------------------

mat3 tokamak_nimrod_parametric::derivative(const double&,
        const vec3& p) const
{
	mat3 r;

	if (!jacobian(p, r)) {		//
		// std::cerr << "jacobian field at " << p << " is undefined\n";
		throw undefined_point();
	}

	return r;
}

// --------------------------------------------------------------------------

Nrrd* tokamak_nimrod_parametric::to_nrrd() const
{
	Nrrd *nrrd = nrrdNew();
	double _spc[4] = { airNaN(), 1., 1., 1. };
	size_t _size[4] = { 3, size[0], size[1], size[2] };
	if (nrrdWrap_nva(nrrd, data, nrrdTypeFloat, 4, _size)) {
		std::cerr << biffGetDone(NRRD) << std::endl;
		exit(-1);
	}
	nrrdAxisInfoSet_nva(nrrd, nrrdAxisInfoSpacing, _spc);
	return nrrd;
}




























































































