#ifndef RAY_H
#define RAY_H

#include "VectorN.h"
#include "DirectionConstants.h"

#ifndef RAY_CLASSIFICATION
#define RAY_CLASSIFICATION
enum RayClassification
{	MMM = 0, MMP = 1, MMO = 2, MPM = 4, MPP = 5, MPO = 6, MOM = 8, MOP = 9, MOO = 10,
	PMM = 16, PMP = 17, PMO = 18, PPM = 20, PPP = 21, PPO = 22, POM = 24, POP = 25, POO = 26,
	OMM = 32, OMP = 33, OMO = 34, OPM = 36, OPP = 37, OPO = 38, OOM = 40, OOP = 41 };
#endif

template <class T>
class Ray {
public:
	Ray() {
		origin.set(0,0,0);
		inverseComputed = false;
		setDirection(Vector3<T>(1,1,1));
	}

	Ray(const Vector3<T> &rayOrigin, const Vector3<T> &rayDirection) {
		origin = rayOrigin;
		inverseComputed = false;
		setDirection(rayDirection);
	}
	Ray(const Ray<T> &r) { *this = r; }

	inline Ray<T>& operator= (const Ray<T> &r) {
		origin = r.origin;
		direction = r.direction;
		inverseComputed = r.inverseComputed;
		if (inverseComputed) {
			invDirection = r.invDirection;
			ibyj = r.ibyj;  jbyi = r.jbyi;  jbyk = r.jbyk;  kbyj = r.kbyj;  ibyk = r.ibyk;  kbyi = r.kbyi;
			c_xy = r.c_xy;  c_xz = r.c_xz;  c_yx = r.c_yx;  c_yz = r.c_yz;  c_zx = r.c_zx;  c_zy = r.c_zy;
			classification = r.classification;
		}
		
		return *this;
	}

	inline void setDirection(const Vector3<T> &dir) { direction = dir.unit(); }
	inline Vector3<T> getDirection() const { return direction; }

	inline void computeDirectionInverse() {
		inverseComputed = true;

		invDirection.x = (T)1.0/direction.x;  invDirection.y = (T)1.0/direction.y;  invDirection.z = (T)1.0/direction.z;

		//ray slope
		ibyj = direction.x * invDirection.y;  jbyi = direction.y * invDirection.x;  jbyk = direction.y * invDirection.z;
		kbyj = direction.z * invDirection.y;  ibyk = direction.x * invDirection.z;  kbyi = direction.z * invDirection.x;
		c_xy = origin.y - jbyi * origin.x;    c_xz = origin.z - kbyi * origin.x;    c_yx = origin.x - ibyj * origin.y;
		c_yz = origin.z - kbyj * origin.y;    c_zx = origin.x - ibyk * origin.z;    c_zy = origin.y - jbyk * origin.z;	

		//ray slope classification
		classification = 0;
		if (direction.z > 0) classification += 1;
		else if (direction.z == 0) classification += 2;
		if (direction.y > 0) classification += 4;
		else if (direction.y == 0) classification += 8;
		if (direction.x > 0) classification += 16;
		else if (direction.x == 0) classification += 32;
	}

	inline void reverseRayDirection(bool computeInverse) {
		setDirection(-direction);
		if (computeInverse) computeDirectionInverse();
		else inverseComputed = false;
	}

	Vector3<T> origin;

	template <class S> inline bool intersectBox(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, float &tNear, float &tFar,
												int &intersectedNearFace, int &intersectedFarFace) const {
		if (inverseComputed) return intersectBox_Eisemann(aabbMin, aabbMax, tNear, tFar, intersectedNearFace, intersectedFarFace);
		else return intersectBox_Williams(aabbMin, aabbMax, tNear, tFar, intersectedNearFace, intersectedFarFace);
	}

	template <class S> inline bool intersectBox(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, float &tNear, float &tFar) const {
		if (inverseComputed) return intersectBox_Eisemann(aabbMin, aabbMax, tNear, tFar);
		else return intersectBox_Williams(aabbMin, aabbMax, tNear, tFar);
	}

	template <class S> inline bool intersectBox(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax) const {
		float tNear, tFar;
		return intersectBox(aabbMin, aabbMax, tNear, tFar);
	}

	template <class S> inline bool intersectBox(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, Vector3<T> &intersectionPt, bool clampToOrigin = true) const {
		float tNear, tFar;
		if (intersectsBox(aabbMin, aabbMax, tNear, tFar)) {
			if (clampToOrigin && tNear < 0) tNear = 0;
			intersectionPt = origin + direction*tNear;
			return true;
		}
		else return false;
	}

	template <class S> inline bool getEndOfBoxIntersection(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, Vector3<T> &intersectionPt, const T perturbation = 0.0001) const {
		float tNear, tFar;
		intersectBox(aabbMin, aabbMax, tNear, tFar);

		if (tFar < 0) return false;
		else if (tFar < 0.001) tFar = 0.001; // if tFar is small, the point might not be changed due to floating point error

		intersectionPt = origin + direction*(tFar+perturbation);
		return true;
	}

	template <class S> inline bool intersectSphere(const Vector3<S> &sphereCenter, const S sphereRadiusSqr, float &tNear, float &tFar) const {
		Vector3f diff = origin - sphereCenter;
		float A = diff.dot(diff) - sphereRadiusSqr;
		float B = diff.dot(direction);
		if (A <= 0) { // inside sphere
			float C = B*B - A;
			tNear = 0;
			tFar = -B + sqrt(C);
			return true;
		}
		else if (B >= 0) return false; // ray origin is behind sphere
		
		float C = B*B - A;
		if (C < 0) return false; // ray does not intersect
		else {
			float root = sqrt(C);
			tNear = -B - root;
			tFar = -B + root;
			return true;
		}
	}

	//template <class S> bool inline intersectPlane(const Vector3<S> &planePt1, const Vector3<S> &planePt2, const Vector3<S> &planePt3, Vector3<T> &intersectionPt) {
	//	Vector3<S> planeNormal = ((planePt2-planePt1).cross(planePt3-planePt1)).unit();
	//	return intersectPlane(planePt1, planePt2, planePt3, planeNormal, true, intersectionPt);
	//}

	//template <class S> bool inline intersectPlane(const Vector3<S> &planePt1, const Vector3<S> &planePt2, const Vector3<S> &planePt3,
	//											  const Vector3<S> &planeNormal, bool normalized, Vector3<T> &intersectionPt) {
	//	float t;
	//	if (intersectPlane(planePt1, planePt2, planePt3, planeNormal, normalized, t)) {
	//		if (t < 0) return false;
	//		else if (t < 0.001) t = 0.001;

	//		intersectionPt = origin + direction*t;
	//		return true;
	//	}
	//	else return false;
	//}

	//template <class S> bool inline intersectPlane(const Vector3<S> &planePt1, const Vector3<S> &planePt2, const Vector3<S> &planePt3,
	//											  const Vector3<S> &planeNormal, bool normalized, float &t) {
	//	float d1 = direction.dot(planeNormal);
	//	float d2;
	//	if (normalized) d2 = origin.signedDistanceToPlane(planePt1, planeNormal);
	//	else d2 = origin.signedDistanceToPlane(planePt1, planeNormal.unit());

	//	const float eps = 0.0001;
	//	if (fabs(d1) > eps) { // ray not parallel to plane
	//		t = -d2/d1;
	//		return true;
	//	}
	//	else if (fabs(d2) <= eps) { // ray and plane are parallel
	//		t = 0;
	//		return true;
	//	}
	//	else return false;


	//	//Real fDdN = m_pkLine->Direction.Dot(m_pkPlane->Normal);
	//	//Real fSDistance = m_pkPlane->DistanceTo(m_pkLine->Origin);
	//	//if (Math<Real>::FAbs(fDdN) > Math<Real>::ZERO_TOLERANCE)
	//	//{
	//	//	// The line is not parallel to the plane, so they must intersect.
	//	//	m_fLineT = -fSDistance/fDdN;
	//	//	m_iIntersectionType = IT_POINT;
	//	//	return true;
	//	//}

	//	//// The Line and plane are parallel.  Determine if they are numerically
	//	//// close enough to be coincident.
	//	//if (Math<Real>::FAbs(fSDistance) <= Math<Real>::ZERO_TOLERANCE)
	//	//{
	//	//	// The line is coincident with the plane, so choose t = 0 for the
	//	//	// parameter.
	//	//	m_fLineT = (Real)0.0;
	//	//	m_iIntersectionType = IT_LINE;
	//	//	return true;
	//	//}

	//	//m_iIntersectionType = IT_EMPTY;
	//	//return false;

	//}


	// http://www.lighthouse3d.com/opengl/maths/index.php?raytriint
	template <class S> inline bool intersectTriangle(const Vector3<S> &triPt1, const Vector3<S> &triPt2, const Vector3<S> &triPt3, float &t) {
		Vector3<S> e1 = triPt2-triPt1;
		Vector3<S> e2 = triPt3-triPt1;

		Vector3<S> h = direction.cross(e2);
		float a = e1.dot(h);

		if (fabs(a) < 0.0001) return false;

		float f = 1/a;

		Vector3<S> s = origin-triPt1;
		float u = f * s.dot(h);

		if (u < 0.0 || u > 1.0) return false;

		Vector3<S> q = s.cross(e1);
		float v = f * direction.dot(q);

		if (v < 0.0 || u + v > 1.0) return false;

		t = f * e2.dot(q);

		if (t > 0.0001) return true;  // ray intersection
		else return false;  // this means that there is a line intersection but not a ray intersection
	}

private:
	Vector3<T> direction;

	bool inverseComputed;
	Vector3<T> invDirection;
	T ibyj, jbyi, jbyk, kbyj, ibyk, kbyi, c_xy, c_xz, c_yx, c_yz, c_zx, c_zy;
	int classification;


	// Ray-box intersection using IEEE numerical properties to ensure that the test is both robust and efficient, as described in:
	//      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley, "An Efficient and Robust Ray-Box Intersection Algorithm"
	//      this returns the correct intersection test and tNear, but fails with tFar when one of the ray directions == 0
	template <class S> inline bool intersectBox_Williams(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, float &tNear, float &tFar) const {
		float tMin, tMax;
		T invDirection = (T)1.0 / direction.x;
		if (invDirection >= 0) { tMin = (aabbMin.x - origin.x) * invDirection;  tMax = (aabbMax.x - origin.x) * invDirection; }
		else				   { tMin = (aabbMax.x - origin.x) * invDirection;  tMax = (aabbMin.x - origin.x) * invDirection; }

		float tMin2, tMax2;
		invDirection = (T)1.0 / direction.y;
		if (invDirection >= 0) { tMin2 = (aabbMin.y - origin.y) * invDirection;  tMax2 = (aabbMax.y - origin.y) * invDirection; }
		else				   { tMin2 = (aabbMax.y - origin.y) * invDirection;  tMax2 = (aabbMin.y - origin.y) * invDirection; }

		if (tMin > tMax2 || tMin2 > tMax) return false;
		if (tMin2 > tMin) tMin = tMin2;
		if (tMax2 < tMax) tMax = tMax2;

		invDirection = (T)1.0 / direction.z;
		if (invDirection >= 0) { tMin2 = (aabbMin.z - origin.z) * invDirection;  tMax2 = (aabbMax.z - origin.z) * invDirection; }
		else				   { tMin2 = (aabbMax.z - origin.z) * invDirection;  tMax2 = (aabbMin.z - origin.z) * invDirection; }

		if (tMin > tMax2 || tMin2 > tMax) return false;
		if (tMin2 > tMin) tMin = tMin2;
		if (tMax2 < tMax) tMax = tMax2;

		tNear = tMin;
		tFar = tMax;

		return tMax > tMin;
	}

	template <class S> inline bool intersectBox_Williams(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, float &tNear, float &tFar,
														 int &intersectedNearFace, int &intersectedFarFace) const {
		float tMin, tMax;
		T invDirection = (T)1.0 / direction.x;
		if (invDirection >= 0) { tMin = (aabbMin.x - origin.x) * invDirection;  tMax = (aabbMax.x - origin.x) * invDirection;
										intersectedNearFace = FACE_LEFT;  intersectedFarFace = FACE_RIGHT;					  }
		else				   { tMin = (aabbMax.x - origin.x) * invDirection;  tMax = (aabbMin.x - origin.x) * invDirection;
										intersectedNearFace = FACE_RIGHT;  intersectedFarFace = FACE_LEFT;					  }

		float tMin2, tMax2;
		invDirection = (T)1.0 / direction.y;
		if (invDirection >= 0) { tMin2 = (aabbMin.y - origin.y) * invDirection;  tMax2 = (aabbMax.y - origin.y) * invDirection; }
		else				   { tMin2 = (aabbMax.y - origin.y) * invDirection;  tMax2 = (aabbMin.y - origin.y) * invDirection; }

		if (tMin > tMax2 || tMin2 > tMax) return false;

		if (tMin2 > tMin) {
			tMin = tMin2;
			if (invDirection >= 0) intersectedNearFace = FACE_DOWN;
			else intersectedNearFace = FACE_UP;
		}
		else if (tMin2 == tMin) {
			if (invDirection >= 0) intersectedNearFace += FACE_DOWN;
			else intersectedNearFace += FACE_UP;
		}

		if (tMax2 < tMax) {
			tMax = tMax2;
			if (invDirection >= 0) intersectedFarFace = FACE_UP;
			else intersectedFarFace = FACE_DOWN;
		}
		else if (tMin2 == tMin) {
			if (invDirection >= 0) intersectedNearFace += FACE_UP;
			else intersectedNearFace += FACE_DOWN;
		}

		invDirection = (T)1.0 / direction.z;
		if (invDirection >= 0) { tMin2 = (aabbMin.z - origin.z) * invDirection;  tMax2 = (aabbMax.z - origin.z) * invDirection; }
		else				   { tMin2 = (aabbMax.z - origin.z) * invDirection;  tMax2 = (aabbMin.z - origin.z) * invDirection; }

		if (tMin > tMax2 || tMin2 > tMax) return false;

		if (tMin2 > tMin) {
			tMin = tMin2;
			if (invDirection >= 0) intersectedNearFace = FACE_BACK;
			else intersectedNearFace = FACE_FRONT;
		}
		else if (tMin2 == tMin) {
			if (invDirection >= 0) intersectedNearFace += FACE_BACK;
			else intersectedNearFace += FACE_FRONT;
		}

		if (tMax2 < tMax) {
			tMax = tMax2;
			if (invDirection >= 0) intersectedFarFace = FACE_FRONT;
			else intersectedFarFace = FACE_BACK;
		}
		else if (tMin2 == tMin) {
			if (invDirection >= 0) intersectedNearFace += FACE_FRONT;
			else intersectedNearFace += FACE_BACK;
		}

		tNear = tMin;
		tFar = tMax;

		return tMax > tMin;
	}

	// "Fast Ray/Axis-Aligned Bounding Box Overlap Tests using Ray Slopes" Martin Eisemann, Thorsten Grosch, Stefan Mueller, Marcus Magnor
	//    this monstrosity of code is modified from the code provided in the above submission
	//    supposedly is faster than the other code above... but at the very least it treats the case of ray direction = 0 correctly
	template <class S> bool inline intersectBox_Eisemann(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, float &tNear, float &tFar) const {
		switch (classification) {
			case MMM: {
				if ((origin.x < aabbMin.x) || (origin.y < aabbMin.y) || (origin.z < aabbMin.z) || (jbyi * aabbMin.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMin.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMin.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMax.x + c_zx > 0)) return false;
				
				// compute the intersection distance
				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}


			case MMP: {		
				if ((origin.x < aabbMin.x) || (origin.y < aabbMin.y) || (origin.z > aabbMax.z) || (jbyi * aabbMin.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMax.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMin.z + c_yz < 0)
					|| (kbyi * aabbMin.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case MPM: {		
				if ((origin.x < aabbMin.x) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (jbyi * aabbMin.x - aabbMin.y + c_xy < 0) 
					|| (ibyj * aabbMax.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMin.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMin.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case MPP: {
				if ((origin.x < aabbMin.x) || (origin.y > aabbMax.y) || (origin.z > aabbMax.z) || (jbyi * aabbMin.x - aabbMin.y + c_xy < 0) 
					|| (ibyj * aabbMax.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMax.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMin.z + c_yz < 0) 
					|| (kbyi * aabbMin.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case PMM: {
				if ((origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z < aabbMin.z) || (jbyi * aabbMax.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMin.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMax.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case PMP: {
				if ((origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z > aabbMax.z) || (jbyi * aabbMax.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMax.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMin.z + c_yz < 0)
					|| (kbyi * aabbMax.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMin.x + c_zx < 0)) return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case PPM: {
				if ((origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (jbyi * aabbMax.x - aabbMin.y + c_xy < 0)
					|| (ibyj * aabbMax.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMin.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMax.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case PPP: {
				if ((origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z > aabbMax.z) || (jbyi * aabbMax.x - aabbMin.y + c_xy < 0)
					|| (ibyj * aabbMax.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMax.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMin.z + c_yz < 0)
					|| (kbyi * aabbMax.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;
				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) tNear = t2;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) tFar = t2;
				return true;
			}

			case OMM: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z < aabbMin.z)
					|| (jbyk * aabbMin.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMax.z + c_yz > 0)) return false;

				tNear = (aabbMax.y - origin.y) * invDirection.y;
				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMin.y - origin.y) * invDirection.y;
				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case OMP: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z > aabbMax.z)
					|| (jbyk * aabbMax.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMin.z + c_yz < 0)) return false;

				tNear = (aabbMax.y - origin.y) * invDirection.y;
				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMin.y - origin.y) * invDirection.y;
				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case OPM: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z)
					|| (jbyk * aabbMin.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMax.z + c_yz > 0)) return false;

				tNear = (aabbMin.y - origin.y) * invDirection.y;		
				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMax.y - origin.y) * invDirection.y;		
				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case OPP: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z > aabbMax.z)
					|| (jbyk * aabbMax.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMin.z + c_yz < 0)) return false;
				
				tNear = (aabbMin.y - origin.y) * invDirection.y;		
				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMax.y - origin.y) * invDirection.y;		
				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case MOM: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x < aabbMin.x) || (origin.z < aabbMin.z) 
					|| (kbyi * aabbMin.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case MOP: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x < aabbMin.x) || (origin.z > aabbMax.z) 
					|| (kbyi * aabbMin.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMax.x + c_zx > 0)) return false;

				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case POM: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x > aabbMax.x) || (origin.z < aabbMin.z)
					|| (kbyi * aabbMax.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}	

			case POP: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x > aabbMax.x) || (origin.z > aabbMax.z)
					|| (kbyi * aabbMax.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMin.x + c_zx < 0)) return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) tFar = t1;
				return true;
			}	

			case MMO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x < aabbMin.x) || (origin.y < aabbMin.y)  
					|| (jbyi * aabbMin.x - aabbMax.y + c_xy > 0) || (ibyj * aabbMin.y - aabbMax.x + c_yx > 0)) return false;

				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				return true;
			}	

			case MPO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x < aabbMin.x) || (origin.y > aabbMax.y) 
					|| (jbyi * aabbMin.x - aabbMin.y + c_xy < 0) || (ibyj * aabbMax.y - aabbMax.x + c_yx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMin.x - origin.x) * invDirection.x;
				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case PMO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) 
					|| (jbyi * aabbMax.x - aabbMax.y + c_xy > 0) || (ibyj * aabbMin.y - aabbMin.x + c_yx < 0)) return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case PPO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x > aabbMax.x) || (origin.y > aabbMax.y) 
					|| (jbyi * aabbMax.x - aabbMin.y + c_xy < 0) || (ibyj * aabbMax.y - aabbMin.x + c_yx < 0)) return false;
			
				tNear = (aabbMin.x - origin.x) * invDirection.x;
				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) tNear = t1;

				tFar = (aabbMax.x - origin.x) * invDirection.x;
				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) tFar = t1;
				return true;
			}

			case MOO: {
				if((origin.x < aabbMin.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;

				tNear = (aabbMax.x - origin.x) * invDirection.x;
				tFar = (aabbMin.x - origin.x) * invDirection.x;
				return true;
			}

			case POO: {
				if((origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;
				tFar = (aabbMax.x - origin.x) * invDirection.x;
				return true;
			}

			case OMO: {
				if((origin.y < aabbMin.y) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;
				
				tNear = (aabbMax.y - origin.y) * invDirection.y;
				tFar = (aabbMin.y - origin.y) * invDirection.y;
				return true;
			}

			case OPO: {
				if((origin.y > aabbMax.y) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;

				tNear = (aabbMin.y - origin.y) * invDirection.y;
				tFar = (aabbMax.y - origin.y) * invDirection.y;
				return true;
			}

			case OOM: {
				if((origin.z < aabbMin.z) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y))
					return false;

				tNear = (aabbMax.z - origin.z) * invDirection.z;
				tFar = (aabbMin.z - origin.z) * invDirection.z;
				return true;
			}

			case OOP: {
				if((origin.z > aabbMax.z) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y))
					return false;

				tNear = (aabbMin.z - origin.z) * invDirection.z;
				tFar = (aabbMax.z - origin.z) * invDirection.z;
				return true;
			}
		}

		return false;
	}

	template <class S> bool inline intersectBox_Eisemann(const Vector3<S> &aabbMin, const Vector3<S> &aabbMax, float &tNear, float &tFar, int &intersectedNearFace, int &intersectedFarFace) const {
		switch (classification) {
			case MMM: {
				if ((origin.x < aabbMin.x) || (origin.y < aabbMin.y) || (origin.z < aabbMin.z) || (jbyi * aabbMin.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMin.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMin.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMax.x + c_zx > 0)) return false;
				
				// compute the intersection distances
				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_UP; }
				else if (t1 == tNear) intersectedNearFace += FACE_UP;

				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_FRONT; }
				else if (t2 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_DOWN; }
				else if (t1 == tFar) intersectedFarFace += FACE_DOWN;

				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_BACK; }
				else if (t2 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}


			case MMP: {		
				if ((origin.x < aabbMin.x) || (origin.y < aabbMin.y) || (origin.z > aabbMax.z) || (jbyi * aabbMin.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMax.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMin.z + c_yz < 0)
					|| (kbyi * aabbMin.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_UP; }
				else if (t1 == tNear) intersectedNearFace += FACE_UP;

				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_BACK; }
				else if (t2 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_DOWN; }
				else if (t1 == tFar) intersectedFarFace += FACE_DOWN;

				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_FRONT; }
				else if (t2 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case MPM: {		
				if ((origin.x < aabbMin.x) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (jbyi * aabbMin.x - aabbMin.y + c_xy < 0) 
					|| (ibyj * aabbMax.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMin.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMin.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_DOWN; }
				else if (t1 == tNear) intersectedNearFace += FACE_DOWN;

				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_FRONT; }
				else if (t2 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_UP; }
				else if (t1 == tFar) intersectedFarFace += FACE_UP;

				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_BACK; }
				else if (t2 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}

			case MPP: {
				if ((origin.x < aabbMin.x) || (origin.y > aabbMax.y) || (origin.z > aabbMax.z) || (jbyi * aabbMin.x - aabbMin.y + c_xy < 0) 
					|| (ibyj * aabbMax.y - aabbMax.x + c_yx > 0) || (jbyk * aabbMax.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMin.z + c_yz < 0) 
					|| (kbyi * aabbMin.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_DOWN; }
				else if (t1 == tNear) intersectedNearFace += FACE_DOWN;

				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_BACK; }
				else if (t2 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_UP; }
				else if (t1 == tFar) intersectedFarFace += FACE_UP;

				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_FRONT; }
				else if (t2 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case PMM: {
				if ((origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z < aabbMin.z) || (jbyi * aabbMax.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMin.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMax.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_UP; }
				else if (t1 == tNear) intersectedNearFace += FACE_UP;

				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_FRONT; }
				else if (t2 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_DOWN; }
				else if (t1 == tFar) intersectedFarFace += FACE_DOWN;

				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_BACK; }
				else if (t2 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}

			case PMP: {
				if ((origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z > aabbMax.z) || (jbyi * aabbMax.x - aabbMax.y + c_xy > 0)
					|| (ibyj * aabbMin.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMax.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMin.z + c_yz < 0)
					|| (kbyi * aabbMax.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMin.x + c_zx < 0)) return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_UP; }
				else if (t1 == tNear) intersectedNearFace += FACE_UP;

				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_BACK; }
				else if (t2 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_DOWN; }
				else if (t1 == tFar) intersectedFarFace += FACE_DOWN;

				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_FRONT; }
				else if (t2 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case PPM: {
				if ((origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (jbyi * aabbMax.x - aabbMin.y + c_xy < 0)
					|| (ibyj * aabbMax.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMin.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMax.z + c_yz > 0)
					|| (kbyi * aabbMax.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_DOWN; }
				else if (t1 == tNear) intersectedNearFace += FACE_DOWN;

				float t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_FRONT; }
				else if (t2 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_UP; }
				else if (t1 == tFar) intersectedFarFace += FACE_UP;

				t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_BACK; }
				else if (t2 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}

			case PPP: {
				if ((origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z > aabbMax.z) || (jbyi * aabbMax.x - aabbMin.y + c_xy < 0)
					|| (ibyj * aabbMax.y - aabbMin.x + c_yx < 0) || (jbyk * aabbMax.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMin.z + c_yz < 0)
					|| (kbyi * aabbMax.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_DOWN; }
				else if (t1 == tNear) intersectedNearFace += FACE_DOWN;

				float t2 = (aabbMin.z - origin.z) * invDirection.z;
				if (t2 > tNear) { tNear = t2;  intersectedNearFace = FACE_BACK; }
				else if (t2 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_UP; }
				else if (t1 == tFar) intersectedFarFace += FACE_UP;

				t2 = (aabbMax.z - origin.z) * invDirection.z;
				if (t2 < tFar) { tFar = t2;  intersectedFarFace = FACE_FRONT; }
				else if (t2 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case OMM: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z < aabbMin.z)
					|| (jbyk * aabbMin.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMax.z + c_yz > 0)) return false;

				tNear = (aabbMax.y - origin.y) * invDirection.y;  intersectedNearFace = FACE_UP;

				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_FRONT; }
				else if (t1 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMin.y - origin.y) * invDirection.y;  intersectedFarFace = FACE_DOWN;

				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_BACK; }
				else if (t1 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}

			case OMP: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.z > aabbMax.z)
					|| (jbyk * aabbMax.z - aabbMax.y + c_zy > 0) || (kbyj * aabbMin.y - aabbMin.z + c_yz < 0)) return false;

				tNear = (aabbMax.y - origin.y) * invDirection.y;  intersectedNearFace = FACE_UP;

				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_BACK; }
				else if (t1 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMin.y - origin.y) * invDirection.y;  intersectedFarFace = FACE_DOWN;

				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_FRONT; }
				else if (t1 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case OPM: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z)
					|| (jbyk * aabbMin.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMax.z + c_yz > 0)) return false;

				tNear = (aabbMin.y - origin.y) * invDirection.y;  intersectedNearFace = FACE_DOWN;

				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_FRONT; }
				else if (t1 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMax.y - origin.y) * invDirection.y;  intersectedFarFace = FACE_UP;

				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_BACK; }
				else if (t1 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}

			case OPP: {
				if((origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y > aabbMax.y) || (origin.z > aabbMax.z)
					|| (jbyk * aabbMax.z - aabbMin.y + c_zy < 0) || (kbyj * aabbMax.y - aabbMin.z + c_yz < 0)) return false;
				
				tNear = (aabbMin.y - origin.y) * invDirection.y;  intersectedNearFace = FACE_DOWN;

				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_BACK; }
				else if (t1 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMax.y - origin.y) * invDirection.y;  intersectedFarFace = FACE_UP;

				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_FRONT; }
				else if (t1 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case MOM: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x < aabbMin.x) || (origin.z < aabbMin.z) 
					|| (kbyi * aabbMin.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMax.x + c_zx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_FRONT; }
				else if (t1 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_BACK; }
				else if (t1 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}

			case MOP: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x < aabbMin.x) || (origin.z > aabbMax.z) 
					|| (kbyi * aabbMin.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMax.x + c_zx > 0)) return false;

				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_BACK; }
				else if (t1 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_FRONT; }
				else if (t1 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}

			case POM: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x > aabbMax.x) || (origin.z < aabbMin.z)
					|| (kbyi * aabbMax.x - aabbMax.z + c_xz > 0) || (ibyk * aabbMin.z - aabbMin.x + c_zx < 0)) return false;
				
				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_FRONT; }
				else if (t1 == tNear) intersectedNearFace += FACE_FRONT;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_BACK; }
				else if (t1 == tFar) intersectedFarFace += FACE_BACK;

				return true;
			}	

			case POP: {
				if((origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.x > aabbMax.x) || (origin.z > aabbMax.z)
					|| (kbyi * aabbMax.x - aabbMin.z + c_xz < 0) || (ibyk * aabbMax.z - aabbMin.x + c_zx < 0)) return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMin.z - origin.z) * invDirection.z;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_BACK; }
				else if (t1 == tNear) intersectedNearFace += FACE_BACK;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMax.z - origin.z) * invDirection.z;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_FRONT; }
				else if (t1 == tFar) intersectedFarFace += FACE_FRONT;

				return true;
			}	

			case MMO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x < aabbMin.x) || (origin.y < aabbMin.y)  
					|| (jbyi * aabbMin.x - aabbMax.y + c_xy > 0) || (ibyj * aabbMin.y - aabbMax.x + c_yx > 0)) return false;

				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_UP; }
				if (t1 == tNear) intersectedNearFace += FACE_UP;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_DOWN; }
				else if (t1 == tFar) intersectedFarFace += FACE_DOWN;

				return true;
			}	

			case MPO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x < aabbMin.x) || (origin.y > aabbMax.y) 
					|| (jbyi * aabbMin.x - aabbMin.y + c_xy < 0) || (ibyj * aabbMax.y - aabbMax.x + c_yx > 0)) return false;
				
				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;

				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_DOWN; }
				else if (t1 == tNear) intersectedNearFace += FACE_DOWN;

				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;

				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_UP; }
				else if (t1 == tFar) intersectedFarFace += FACE_UP;

				return true;
			}

			case PMO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) 
					|| (jbyi * aabbMax.x - aabbMax.y + c_xy > 0) || (ibyj * aabbMin.y - aabbMin.x + c_yx < 0)) return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_UP; }
				if (t1 == tNear) intersectedNearFace += FACE_UP;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_DOWN; }
				else if (t1 == tFar) intersectedFarFace += FACE_DOWN;

				return true;
			}

			case PPO: {
				if((origin.z < aabbMin.z) || (origin.z > aabbMax.z) || (origin.x > aabbMax.x) || (origin.y > aabbMax.y) 
					|| (jbyi * aabbMax.x - aabbMin.y + c_xy < 0) || (ibyj * aabbMax.y - aabbMin.x + c_yx < 0)) return false;
			
				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;

				float t1 = (aabbMin.y - origin.y) * invDirection.y;
				if (t1 > tNear) { tNear = t1;  intersectedNearFace = FACE_DOWN; }
				else if (t1 == tNear) intersectedNearFace += FACE_DOWN;

				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;

				t1 = (aabbMax.y - origin.y) * invDirection.y;
				if (t1 < tFar) { tFar = t1;  intersectedFarFace = FACE_UP; }
				else if (t1 == tFar) intersectedFarFace += FACE_UP;

				return true;
			}

			case MOO: {
				if((origin.x < aabbMin.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;

				tNear = (aabbMax.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_RIGHT;
				tFar = (aabbMin.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_LEFT;
				return true;
			}

			case POO: {
				if((origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;

				tNear = (aabbMin.x - origin.x) * invDirection.x;  intersectedNearFace = FACE_LEFT;
				tFar = (aabbMax.x - origin.x) * invDirection.x;  intersectedFarFace = FACE_RIGHT;
				return true;
			}

			case OMO: {
				if((origin.y < aabbMin.y) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;
				
				tNear = (aabbMax.y - origin.y) * invDirection.y;  intersectedNearFace = FACE_UP;
				tFar = (aabbMin.y - origin.y) * invDirection.y;  intersectedFarFace = FACE_DOWN;
				return true;
			}

			case OPO: {
				if((origin.y > aabbMax.y) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.z < aabbMin.z) || (origin.z > aabbMax.z))
					return false;

				tNear = (aabbMin.y - origin.y) * invDirection.y;  intersectedNearFace = FACE_DOWN;
				tFar = (aabbMax.y - origin.y) * invDirection.y;  intersectedFarFace = FACE_UP;
				return true;
			}

			case OOM: {
				if((origin.z < aabbMin.z) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y))
					return false;

				tNear = (aabbMax.z - origin.z) * invDirection.z;  intersectedNearFace = FACE_FRONT;
				tFar = (aabbMin.z - origin.z) * invDirection.z;  intersectedFarFace = FACE_BACK;
				return true;
			}

			case OOP: {
				if((origin.z > aabbMax.z) || (origin.x < aabbMin.x) || (origin.x > aabbMax.x) || (origin.y < aabbMin.y) || (origin.y > aabbMax.y))
					return false;

				tNear = (aabbMin.z - origin.z) * invDirection.z;  intersectedNearFace = FACE_BACK;
				tFar = (aabbMax.z - origin.z) * invDirection.z;  intersectedFarFace = FACE_FRONT;
				return true;
			}
		}

		return false;
	}
};

#endif
