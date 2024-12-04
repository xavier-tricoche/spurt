#ifndef __MLS_HPP__
#define __MLS_HPP__

#include <vector>
#include <set>
#include <list>
#include <math.h>
#include <iostream>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>

// using LAPACK for SVD-based LS solution
#include <clapack.h>

namespace MLS
{
// polynomial precision
static unsigned int dof[2][4] = {
	{ 1, 3, 6, 10 },
	{ 1, 4, 10, 20 }
};

// smooth (C^2) in dim 2 and 3
inline double wendland(double r, double radius)
{
	double __r = r / radius;
	if (__r >= 1.) return 0;
	return pow(1. - __r, 4.)*(4.*__r + 1);
}

inline double franke_nielson(double r, double radius)
{
	double d = radius - r;
	if (d <= 0) return 0;
	d /= radius * r;
}

template< typename V >
static inline void linear_2d(std::vector<double>& b, const V& p)
{
	b[0] = 1.;
	b[1] = p[0];
	b[2] = p[1];
}

template< typename V >
static inline void linear_3d(std::vector<double>& b, const V& p)
{
	b[0] = 1.;
	b[1] = p[0];
	b[2] = p[1];
	b[3] = p[2];
}

template< typename V >
static inline void quadratic_2d(std::vector<double>& b, const V& p)
{
	linear_2d(b, p);
	b[3] = p[0] * p[0];
	b[4] = p[0] * p[1];
	b[5] = p[1] * p[1];
}

template< typename V >
static inline void quadratic_3d(std::vector<double>& b, const V& p)
{
	linear_3d(b, p);
	b[4] = p[0] * p[0];
	b[5] = p[0] * p[1];
	b[6] = p[0] * p[2];
	b[7] = p[1] * p[1];
	b[8] = p[1] * p[2];
	b[9] = p[2] * p[2];
}

template< typename V >
inline void cubic_2d(std::vector<double>& b, const V& p)
{
	quadratic_2d(b, p);
	b[6] = p[0] * p[0] * p[0];
	b[7] = p[0] * p[0] * p[1];
	b[8] = p[0] * p[1] * p[1];
	b[9] = p[1] * p[1] * p[1];
}

template< typename V >
inline void cubic_3d(std::vector<double>& b, const V& p)
{
	quadratic_3d(b, p);
	b[10] = p[0] * p[0] * p[0];
	b[11] = p[0] * p[0] * p[1];
	b[12] = p[0] * p[0] * p[2];
	b[13] = p[0] * p[1] * p[1];
	b[14] = p[0] * p[1] * p[2];
	b[15] = p[0] * p[2] * p[2];
	b[16] = p[1] * p[1] * p[1];
	b[17] = p[1] * p[1] * p[2];
	b[18] = p[1] * p[2] * p[2];
	b[19] = p[2] * p[2] * p[2];
}

template<typename V>
inline void set_basis(std::vector<double>& b, const V& p, int dim, int prec)
{
	if (prec == 0) {
		b[0] = 1.;
	}
	else if (dim == 2) {
		switch (prec) {
		case 1:
			linear_2d(b, p);
			break;
		case 2:
			quadratic_2d(b, p);
			break;
		case 3:
			cubic_2d(b, p);
			break;
		default:
			assert(false);
		}
	}
	else {
		switch (prec) {
		case 1:
			linear_3d(b, p);
			break;
		case 2:
			quadratic_3d(b, p);
			break;
		case 3:
			cubic_3d(b, p);
			break;
		default:
			assert(false);
		}
	}
}

};

namespace MLS
{
typedef __CLPK_integer 		integer_type;
typedef __CLPK_doublereal 	float_type;

struct CLAPACK_helper {
	float_type *A, *b, *work, *s, rcond;
	integer_type ncols, nrows, nrhs, lwork, info, rank;
	bool svd;
	char trans;
	unsigned int dim, prec;

	CLAPACK_helper(unsigned int max_nvals,
	               unsigned int _dim,
	               unsigned int _prec, unsigned int _nrhs,
	               bool use_svd = false) {
		dim = _dim;
		prec = _prec;
		nrows = max_nvals;
		ncols = MLS::dof[dim-2][prec];
		nrhs = _nrhs;
		lwork = -1;
		trans = 'N';
		svd = use_svd;
		A = (float_type*)malloc(nrows * ncols * sizeof(float_type));
		b = (float_type*)malloc(nrows * nrhs * sizeof(float_type));
		if (use_svd) {
			s = (float_type*)malloc(ncols * sizeof(float_type));
			rcond = 1.0e-4;
		}
		work = (float_type*)malloc(10 * sizeof(float_type));
		if (!use_svd)
			dgels_(&trans, &nrows, &ncols, &nrhs, A, &nrows, b, &nrows, work, &lwork, &info);
		else
			dgelss_(&nrows, &ncols, &nrhs, A, &nrows, b, &nrows, s, &rcond, &rank,
			        work, &lwork, &info);
		lwork = (integer_type)work[0];
		std::cout << "optimal work space is " << lwork << std::endl;
		delete[] work;
		work = (float_type*)malloc(lwork * sizeof(float_type));
	}

	~CLAPACK_helper() {
		delete[] A;
		delete[] b;
		delete[] work;
		delete[] s;
	}
};

}

namespace MLS
{

template< typename T >
inline int MLS(std::vector<double>& coef,
               const std::vector<nvis::vec3>& points,
               const std::vector<T>& values, const nvis::vec3& x0,
               CLAPACK_helper& helper, double radius)
{
	typedef T value_type;

	unsigned int dim = helper.dim;
	unsigned int prec = helper.prec;
	integer_type nbasisfn = helper.ncols;
	integer_type npts = points.size();

	std::vector< double > basis(nbasisfn);

	assert(npts <= helper.nrows); // we have enough pre-allocated space to compute
	assert(points.size() <= values.size()); // we have enough data
	assert(dim >= 2 && dim <= 3); // dimensions make sense
	assert(prec <= 3); // precision is supported

	nvis::vec3 x;
	for (unsigned int i = 0 ; i < points.size() ; ++i) {
		x = points[i] - x0;
		double r = nvis::norm(x);
		MLS::set_basis(basis, x, dim, prec);
		double theta = MLS::wendland(r, radius);

		for (unsigned int n = 0 ; n < nbasisfn ; ++n) {
			helper.A[i+n*npts] = theta * basis[n];
		}
		for (unsigned int n = 0 ; n < helper.nrhs ; ++n) {
			helper.b[i+n*npts] = theta * values[i][n];
		}
	}

	// solve with QR/LQ assuming full rank of A
	helper.rank = helper.ncols; // full rank assumed by default
	if (!helper.svd) {
		dgels_(&helper.trans, &npts, &helper.ncols, &helper.nrhs,
		       helper.A, &npts, helper.b, &npts, helper.work, &helper.lwork,
		       &helper.info);
	}
	else {
		dgelss_(&npts, &helper.ncols, &helper.nrhs, helper.A, &npts, helper.b, &npts,
		        helper.s, &helper.rcond, &helper.rank,
		        helper.work, &helper.lwork, &helper.info);
	}

	if (helper.info) {
		std::cout << "unable to solve LS problem at " << x0 << std::endl;
		std::cout << "we had " << points.size() << " neighbors" << std::endl;
	}
	if (helper.rank < helper.ncols) {
		std::cout << "SVD found LS matrix to be rank deficient: " << helper.rank << std::endl;
		std::cout << npts << " points were available" << std::endl;
	}

	coef.resize(helper.nrhs*dof[dim-2][prec]);
	std::fill(coef.begin(), coef.end(), 0);
	for (unsigned int j = 0 ; j < helper.nrhs ; ++j) {
		for (unsigned int i = 0 ; i < helper.rank ; ++i) {
			coef[j*nbasisfn+i] = helper.b[j*npts+i];
		}
	}
	// remaining coefficients are set to 0

	return helper.rank;
}



template< typename T, typename V >
inline int LeastSquares(std::vector< std::vector<double> >& coef,
                        const std::vector<V>& points,
                        const std::vector<T>& values, const V& x0,
                        CLAPACK_helper& helper)
{
	typedef T								value_type;

	unsigned int dim = helper.dim;
	unsigned int prec = helper.prec;
	integer_type nbasisfn = helper.ncols;
	integer_type npts = points.size();

	std::vector< double > basis(nbasisfn);

	assert(npts <= helper.nrows); // we have enough pre-allocated space to compute
	assert(points.size() <= values.size()); // we have enough data
	assert(dim >= 2 && dim <= 3); // dimensions make sense
	assert(prec <= 3); // precision is supported

	V x;
	for (unsigned int i = 0 ; i < points.size() ; ++i) {
		x = points[i] - x0;
		MLS::set_basis(basis, x, dim, prec);

		for (unsigned int n = 0 ; n < nbasisfn ; ++n) {
			helper.A[i+n*npts] = basis[n];
		}
		for (unsigned int n = 0 ; n < helper.nrhs ; ++n) {
			helper.b[i+n*npts] = values[i][n];
		}
	}

	// solve with QR/LQ assuming full rank of A
	helper.rank = helper.ncols; // full rank assumed by default
	if (!helper.svd) {
		dgels_(&helper.trans, &npts, &helper.ncols, &helper.nrhs,
		       helper.A, &npts, helper.b, &npts, helper.work, &helper.lwork,
		       &helper.info);
	}
	else {
		dgelss_(&npts, &helper.ncols, &helper.nrhs, helper.A, &npts, helper.b, &npts,
		        helper.s, &helper.rcond, &helper.rank,
		        helper.work, &helper.lwork, &helper.info);
	}

	if (helper.info) {
		std::cout << "unable to solve LS problem at " << x0 << std::endl;
		std::cout << "we had " << points.size() << " neighbors" << std::endl;
	}
	if (helper.rank < helper.ncols) {
		std::cout << "SVD found LS matrix to be rank deficient: " << helper.rank << std::endl;
		std::cout << npts << " points were available" << std::endl;
	}

	for (unsigned int j = 0 ; j < helper.nrhs ; ++j) {
		for (unsigned int i = 0 ; i < nbasisfn ; ++i) {
			coef[j][i] = helper.b[j*npts+i];
		}
	}

	return helper.rank;
}





}



#endif




























































































