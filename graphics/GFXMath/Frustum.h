#ifndef FRUSTUM_H
#define FRUSTUM_H

#include "VectorN.h"

#ifndef FRUSTUM_CONSTANTS
#define FRUSTUM_CONSTANTS
	enum FrustumSides { FRUSTUM_TOP = 0, FRUSTUM_BOTTOM, FRUSTUM_LEFT, FRUSTUM_RIGHT, FRUSTUM_NEAR, FRUSTUM_FAR };
#endif

// to do: take a look at the LightHouse3D tutorial on frustum culling

class Frustum {
public:
	Frustum() { }
	Frustum(const Vector3f &left, const Vector3f &right, const Vector3f &bottom, const Vector3f &top, const Vector3f &_near, const Vector3f &_far)
		{ set(left, right, bottom, top, _near, _far); }

	inline void set(const Vector3f &left, const Vector3f &right, const Vector3f &bottom, const Vector3f &top, const Vector3f &_near, const Vector3f &_far, bool doNormalization = true) {
		pts[FRUSTUM_TOP] = top;
		pts[FRUSTUM_BOTTOM] = bottom;
		pts[FRUSTUM_LEFT] = left;
		pts[FRUSTUM_RIGHT] = right;
		pts[FRUSTUM_NEAR] = _near;
		pts[FRUSTUM_FAR] = _far;
		
		if (doNormalization) normalize();
	}

	inline void normalize() {
		pts[FRUSTUM_TOP].normalize();
		pts[FRUSTUM_BOTTOM].normalize();
		pts[FRUSTUM_LEFT].normalize();
		pts[FRUSTUM_RIGHT].normalize();
		pts[FRUSTUM_NEAR].normalize();
		pts[FRUSTUM_FAR].normalize();
	}

	Vector3f pts[6];
};

#endif