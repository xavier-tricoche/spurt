#ifndef GEOMETRIC_UTILS_H
#define GEOMETRIC_UTILS_H

#include "VectorN.h"

#include <vector>
using namespace std;

namespace NearestNeighbor {
	static inline unsigned int get(const Vector3f &pt, const vector<Vector3f> &neighborPts) {
		unsigned int nearestPoint = 0;
		float minDistSq = (pt-neighborPts[0]).magnitudeSquared();

		for (unsigned int i=1; i<neighborPts.size(); i++) {
			float distSq = (pt-neighborPts[i]).magnitudeSquared();
			if (distSq < minDistSq) {
				minDistSq = distSq;
				nearestPoint = i;
			}
		}

		return nearestPoint;
	}

	static inline unsigned int get(const Vector3f &pt, const vector<unsigned int> &neighborPtIndices, const vector<Vector3f> &allPts) {
		unsigned int nearestPoint = neighborPtIndices[0];
		float minDistSq = (pt-allPts[neighborPtIndices[0]]).magnitudeSquared();

		for (unsigned int i=1; i<neighborPtIndices.size(); i++) {
			float distSq = (pt-allPts[neighborPtIndices[i]]).magnitudeSquared();
			if (distSq < minDistSq) {
				minDistSq = distSq;
				nearestPoint = neighborPtIndices[i];
			}
		}

		return nearestPoint;
	}
}

namespace Intersection {
	bool boxPlaneIntersection(const Vector3f &boxMin, const Vector3f &boxMax, const Vector3f &pointOnPlane, const Vector3f &planeNormal) {
		// not the fastest way to do this...
		float dist1 = boxMin.signedDistanceToPlane(pointOnPlane, planeNormal);
		float dist2 = boxMax.signedDistanceToPlane(pointOnPlane, planeNormal);
		if (dist1 < 0 && dist2 < 0) {
			for (unsigned int i=0; i<6; i++) {
				Vector3f pt;
				if (i == 0) pt.set(boxMin.x, boxMin.y, boxMax.z);
				else if (i == 1) pt.set(boxMin.x, boxMax.y, boxMin.z);
				else if (i == 2) pt.set(boxMax.x, boxMin.y, boxMin.z);
				else if (i == 3) pt.set(boxMax.x, boxMax.y, boxMin.z);
				else if (i == 4) pt.set(boxMax.x, boxMin.y, boxMax.z);
				else/*i == 5)*/ pt.set(boxMin.x, boxMax.y, boxMax.z);

				float dist = pt.signedDistanceToPlane(pointOnPlane, planeNormal);
				if (dist >= 0) return true;
			}
			return false;
		}
		else if (dist1 > 0 && dist2 > 0) {
			for (unsigned int i=0; i<6; i++) {
				Vector3f pt;
				if (i == 0) pt.set(boxMin.x, boxMin.y, boxMax.z);
				else if (i == 1) pt.set(boxMin.x, boxMax.y, boxMin.z);
				else if (i == 2) pt.set(boxMax.x, boxMin.y, boxMin.z);
				else if (i == 3) pt.set(boxMax.x, boxMax.y, boxMin.z);
				else if (i == 4) pt.set(boxMax.x, boxMin.y, boxMax.z);
				else/*i == 5)*/ pt.set(boxMin.x, boxMax.y, boxMax.z);

				float dist = pt.signedDistanceToPlane(pointOnPlane, planeNormal);
				if (dist <= 0) return true;
			}
			return false;
		}
		else return true;
	}

	vector<Vector3f> getTrianglePlaneIntersectionPoints(const Vector3f &tri1, const Vector3f &tri2, const Vector3f &tri3,
														const Vector3f &pointOnPlane, const Vector3f &planeNormal, const float epsilon = 0.0001) {
		vector<Vector3f> intersectionPoints;

		float signedDist[3];
		signedDist[0] = tri1.signedDistanceToPlane(pointOnPlane, planeNormal);
		signedDist[1] = tri2.signedDistanceToPlane(pointOnPlane, planeNormal);
		signedDist[2] = tri3.signedDistanceToPlane(pointOnPlane, planeNormal);

		if (fabs(signedDist[0]) < epsilon) signedDist[0] = 0;
		if (fabs(signedDist[1]) < epsilon) signedDist[1] = 0;
		if (fabs(signedDist[2]) < epsilon) signedDist[2] = 0;

		if (signedDist[0] > 0) {
			if (signedDist[1] > 0) {
				if (signedDist[2] > 0) { } // +++
				else if (signedDist[2] < 0) { // ++-
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[2]))*(tri3-tri1);
					intersectionPoints[1] = tri2+(signedDist[1]/(signedDist[1]-signedDist[2]))*(tri3-tri2);
				}
				else { intersectionPoints.resize(1);  intersectionPoints[0] = tri3; } // ++0
			}
			else if (signedDist[1] < 0) {
				if (signedDist[2] > 0) { // +-+
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[1]))*(tri2-tri1);
					intersectionPoints[1] = tri2+(signedDist[1]/(signedDist[1]-signedDist[2]))*(tri3-tri2);
				}
				else if (signedDist[2] < 0) { // +--
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[1]))*(tri2-tri1);
					intersectionPoints[1] = tri1+(signedDist[0]/(signedDist[0]-signedDist[2]))*(tri3-tri1);
				}
				else { // +-0
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[1]))*(tri2-tri1);
					intersectionPoints[1] = tri3;
				}
			}
			else {
				if (signedDist[2] > 0) { intersectionPoints.resize(1);  intersectionPoints[0] = tri2; } // +0+
				else if (signedDist[2] < 0) { // +0-
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[2]))*(tri3-tri1);
					intersectionPoints[1] = tri2;
				}
				else { intersectionPoints.resize(2);  intersectionPoints[0] = tri2;  intersectionPoints[1]= tri3; } // +00
			}
		}
		else if (signedDist[0] < 0) {
			if (signedDist[1] > 0) {
				if (signedDist[2] > 0) { // -++
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[1]))*(tri2-tri1);
					intersectionPoints[1] = tri1+(signedDist[0]/(signedDist[0]-signedDist[2]))*(tri3-tri1);
				}
				else if (signedDist[2] < 0) { // -+-
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[1]))*(tri2-tri1);
					intersectionPoints[1] = tri2+(signedDist[1]/(signedDist[1]-signedDist[2]))*(tri3-tri2);
				}
				else { // -+0
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[1]))*(tri2-tri1);
					intersectionPoints[1] = tri3;
				}
			}
			else if (signedDist[1] < 0) {
				if (signedDist[2] > 0) { // --+
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[2]))*(tri3-tri1);
					intersectionPoints[1] = tri2+(signedDist[1]/(signedDist[1]-signedDist[2]))*(tri3-tri2);
				}
				else if (signedDist[2] < 0) { } // ---
				else { intersectionPoints.resize(1);  intersectionPoints[0] = tri3; } // --0
			}
			else {
				if (signedDist[2] > 0) { // -0+
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri1+(signedDist[0]/(signedDist[0]-signedDist[2]))*(tri3-tri1);
					intersectionPoints[1] = tri2;
				}
				else if (signedDist[2] < 0) { intersectionPoints.resize(1);  intersectionPoints[0] = tri2; } // -0-
				else { intersectionPoints.resize(2);  intersectionPoints[0] = tri2;  intersectionPoints[1] = tri3; } // -00
			}
		}
		else
		{
			if (signedDist[1] > 0) {
				if (signedDist[2] > 0) { intersectionPoints.resize(1);  intersectionPoints[0] = tri1; } // 0++
				else if (signedDist[2] < 0) { // 0+-
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri2+(signedDist[1]/(signedDist[1]-signedDist[2]))*(tri3-tri2);
					intersectionPoints[1] = tri1;
				}
				else { intersectionPoints.resize(2);  intersectionPoints[0] = tri1;  intersectionPoints[1] = tri3; } // 0+0
			}
			else if (signedDist[1] < 0) {
				if (signedDist[2] > 0) { // 0-+
					intersectionPoints.resize(2);
					intersectionPoints[0] = tri2+(signedDist[1]/(signedDist[1]-signedDist[2]))*(tri3-tri2);
					intersectionPoints[1] = tri1;
				}
				else if (signedDist[2] < 0) { intersectionPoints.resize(1);  intersectionPoints[0] = tri1; } // 0--
				else { intersectionPoints.resize(2);  intersectionPoints[0] = tri1;  intersectionPoints[1] = tri3; } // 0-0
			}
			else {
				if (signedDist[2] > 0) { intersectionPoints.resize(2);  intersectionPoints[0] = tri1;  intersectionPoints[1] = tri2; } // 00+
				else if (signedDist[2] < 0) { intersectionPoints.resize(2);  intersectionPoints[0] = tri1;  intersectionPoints[1] = tri2; } // 00-
				else { // 000
					intersectionPoints.resize(3);  intersectionPoints[0] = tri1;
					intersectionPoints[1] = tri2;  intersectionPoints[2] = tri3;
				}
			}
		}

		return intersectionPoints;
	}

	vector<Vector3f> getTrianglePlaneIntersectionPoints(const Vector3f &tri1, const Vector3f &tri2, const Vector3f &tri3,
														const Vector3f &pointOnPlane1, const Vector3f &pointOnPlane2,
														const Vector3f &pointOnPlane3, const float epsilon = 0.0001) {
		return getTrianglePlaneIntersectionPoints(tri1, tri2, tri3, pointOnPlane1,
												 ((pointOnPlane2-pointOnPlane1).cross(pointOnPlane3-pointOnPlane1)).unit(), epsilon);
	}
}

#endif