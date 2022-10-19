#include <teem/gage.h>
#include <teem/nrrd.h>

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <fstream>

#include <math/fixed_matrix.hpp>
#include "tokamak_nrrd.hpp"

using namespace nvis;


#define USE_TRILINEAR
// #define THETA_NORMALIZED

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

tokamak_nrrd::tokamak_nrrd(const std::string& file) : periodic_coords(true)
{
	load(file, points, data, size);
	precompute_jacobian();
}

void tokamak_nrrd::precompute_jacobian()
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
		std::cout << "tokamak_nrrd::compute_jacobian:" << std::endl
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

void tokamak_nrrd::load(const std::string& file,
                        fvec3*& points, fvec3*& data, uvec3& size)
{
	Nrrd *nin = nrrdNew();
	if (nrrdLoad(nin, file.c_str(), NULL)) {
		std::cerr << "tokamak_nrrd::load: " << biffGetDone(NRRD) << std::endl;
		throw;
	}
	float *__data = (float*)nin->data;

	size = uvec3(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
	double hx = nin->axis[1].spacing;
	double hy = nin->axis[2].spacing;
	double hz = nin->axis[3].spacing;
	
	points = new fvec3[size[0]*size[1]*size[2]];
	data = new fvec3[size[0]*size[1]*size[2]];
	for (int k=0 ; k<size[2] ; ++k) {
		for (int j=0 ; j<size[1] ; ++j) {
			for (int i=0 ; i<size[0] ; ++i) {
				int id = i+size[0]*(j+size[1]*k);
				points[id] = fvec3(i*hx, j*hy, k*hz);
				data[id][0] = __data[3*id];
				data[id][1] = __data[3*id+1];
				data[id][2] = __data[3*id+2];
			}
		}
	}
}

// --------------------------------------------------------------------------

double nrrd_mod(double a, double b)
{
	return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}

vec2 tokamak_nrrd::project(const vec3& p) const
{
	return poloidal_project(p);
}

vec3 tokamak_nrrd::unproject(const vec2& p) const
{
	return poloidal_unproject(p);
}

vec2 tokamak_nrrd::poloidal_project(const nvis::vec3& v) const
{
	return (periodic_coords ?
	        nvis::vec2(nrrd_mod(v[1], size[1] - 1),
	                   nrrd_mod(v[0], size[0] - 1)) :
	        nvis::vec2(v[1], v[0]));
}

vec3 tokamak_nrrd::poloidal_unproject(const nvis::vec2& v) const
{
	return nvis::vec3(v[1], v[0], 0.0);
}

mat2 tokamak_nrrd::project(const mat3& m) const
{
	return poloidal_project(m);
}

mat2 tokamak_nrrd::poloidal_project(const nvis::mat3& m) const
{
	mat2 r;
	r[0][0] = m[1][1];
	r[0][1] = m[1][0];
	r[1][0] = m[0][1];
	r[1][1] = m[0][0];

	return r;
}

void tokamak_nrrd::get_slice(std::vector<nvis::vec2>& parm,
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

bool tokamak_nrrd::intersect(const nvis::vec3& p0, const nvis::vec3& p1) const
{
	int z0 = (int)trunc(p0[2] / (size[2] - 1));
	int z1 = (int)trunc(p1[2] / (size[2] - 1));

	return z0 != z1;
}

// --------------------------------------------------------------------------

bool tokamak_nrrd::interpolate(const vec3& p, vec3& r, const fvec3* variable) const
{
	if (std::isnan(p[0])) {
		throw std::runtime_error("invalid coordinate");
	}

	if (p[0] > size[0] - 1 || p[0] < 0)
		return false;

	int cx = (int)floor(p[0]);
	int cy = (int)floor(nrrd_mod(p[1], size[1] - 1));
	int cz = (int)floor(nrrd_mod(p[2], size[2] - 1));

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

bool tokamak_nrrd::jacobian(const vec3& p, mat3& r, const fmat3* variable) const
{
	if (p[0] > size[0] - 1 || p[0] < 0)
		return false;

	int cx = (int)floor(p[0]);
	int cy = (int)floor(nrrd_mod(p[1], size[1] - 1));
	int cz = (int)floor(nrrd_mod(p[2], size[2] - 1));

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

vec3 tokamak_nrrd::operator()(const double&,
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

mat3 tokamak_nrrd::derivative(const double&,
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





























































































