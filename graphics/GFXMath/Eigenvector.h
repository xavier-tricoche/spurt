#ifndef EIGENVECTOR_H
#define EIGENVECTOR_H

#include "VectorN.h"

class Eigenvector {
public:
	Eigenvector() { isReal = false; }
	template <class T> Eigenvector(const T eigenvalue, const T *eigenvector, const unsigned int vectorDim) {
		set(eigenvalue, eigenvector, vectorDim);
	}

	template <class T> inline void set(const T eigenvalue, const T *eigenvector, const unsigned int vectorDim) {
		isReal = true;
		value = eigenvalue;

		dim = vectorDim;
		for (unsigned int i=0; i<vectorDim; i++) vec[i] = eigenvector[i];
	}

	template <class T> inline void setEigenvector(const T *eigenvector)
		{ vec[0] = eigenvector[0];  vec[1] = eigenvector[1];  vec[2] = eigenvector[2]; }
	template <class T> inline void setEigenvector(const T val0, const T val1, const T val2)
		{ vec[0] = val0;  vec[1] = val1;  vec[2] = val2; }

	bool isReal;
	double value;

	unsigned int dim;
	double vec[3];
};

#endif