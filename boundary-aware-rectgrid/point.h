#pragma once

#include <Eigen/Eigen>

typedef Eigen::Vector<double, 3> Point3;

Point3 interpolate_linear(const Point3& p0, const Point3& p1, double u) {
	return (1.-u)*p0 + u*p1;
}

Point3 crossProduct(const Point3& p0, const Point3& p1) {
	return p0.cross(p1);
}

double dotProduct(const Point3& p0, const Point3& p1) {
	return p0.dot(p1);
}