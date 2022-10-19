#pragma once

#include <vector>
using namespace std;
struct Point {
	double x; double y; double z;
	Point() { x = 0; y = 0; z = 0; } //default constructor
	Point(double x_, double y_, double z_) { //initializing constructor
		x = x_; y = y_; z = z_;
	}
	Point(double* p) {
		x = p[0]; y = p[1]; z = p[2];
	}
	Point operator-(const Point& o) {
		Point v;
		v.x = x - o.x;
		v.y = y - o.y;
		v.z = z - o.z;
		return v;
	}
	Point operator/(const double& n) {
		Point v;
		v.x = x / n;
		v.y = y / n;
		v.z = z / n;
		return v;
	}
	void set(double x_, double y_, double z_) {
		x = x_; y = y_; z = z_;
	}
};

static Point crossProduct(const Point& A, const Point& B)
{
	Point C;
	C.x = (A.y * B.z) - (A.z * B.y);
	C.y = (A.z * B.x) - (A.x * B.z);
	C.z = (A.x * B.y) - (A.y * B.x);
	return C;
}

static double dotproduct(const Point& pt1, const Point& pt2)
{
	/*double dir;
	dir = pt1.x * pt2.x + pt1.y * pt2.y + pt1.z * pt2.z;
	return dir; */
	return pt1.x * pt2.x + pt1.y * pt2.y + pt1.z * pt2.z;
}

//interpolation of the two points 
static Point interpolate_linear(const Point& p1, const Point& p2, double u)
{
	Point inter_pt;
	inter_pt.x = (1.0 - u) * p1.x + u * p2.x;
	inter_pt.y = (1.0 - u) * p1.y + u * p2.y;
	inter_pt.z = (1.0 - u) * p1.z + u * p2.z;
	return inter_pt;
}
