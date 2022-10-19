/*  Nate Andrysco
	November 2005 */

#ifndef VECTOR_N_H
#define VECTOR_N_H

#include <cmath>
#include <iostream>
#include <float.h>

#include "DirectionConstants.h"

using namespace std;

template <class T> class Vector2;
template <class T> class Vector3;
template <class T> class Vector4;

// to do: figure out how to do template specialization when T is int and you use a divide

template <class T> class Vector2
{
public:
	typedef T	value_type;
	
	union {
		struct {
			T x, y;
		};
		struct {
			T s, t;
		};
		T data[2];
	};

	Vector2() {
		data[0] = data[1] = 0;
	}
	Vector2(const T v0, const T v1) {
		data[0] = v0;
		data[1] = v1;
	}
	template <class S> Vector2<T>(const Vector2<S>& copy) {
		data[0] = (T)copy.data[0];
		data[1] = (T)copy.data[1];
	}
	template <class S> Vector2<T>(S val[2]) {
		data[0] = (T)val[0];
		data[1] = (T)val[1];
	}

	inline unsigned int getDim() const {
		return 2;
	}

	inline T operator()(int pos) const {
		if (pos == 0) return data[0];
		if (pos == 1) return data[1];
		fprintf(stderr, "Vector: Illegal data for () operator\n");
		return 0;
	}

	inline void set(T v0, T v1) {
		data[0] = v0;
		data[1] = v1;
	}

	inline T magnitude() const {
		return T(sqrt(data[0]*data[0] + data[1]*data[1]));
	}
	inline T magnitudeSquared() const {
		return data[0]*data[0] + data[1]*data[1];
	}
	inline Vector2<T> unit() const {
		Vector2<T> returnVal(data[0], data[1]);
		returnVal.normalize();
		return returnVal;
	}
	inline void normalize() {
		T magSqr = magnitudeSquared();
		if (magSqr != 0) {
			T magInv = (T)1.0 / sqrt(magSqr);
			*this *= magInv;
		}
	}

	inline T dot(const Vector2<T> &right) const {
		return (data[0]*right.data[0] + data[1]*right.data[1]);
	}
	inline T cross(const Vector2<T> &right) const {
		return data[0]*right.data[1] - data[1]*right.data[0];
	}

	inline Vector2<T>  operator-	() const {
		return Vector2<T>(-data[0], -data[1]);
	}

	inline friend Vector2<T> operator+ (T scalar, const Vector2<T> &v) {
		return Vector2<T>(scalar + v.data[0], scalar + v.data[1]);
	}
	inline friend Vector2<T> operator- (T scalar, const Vector2<T> &v) {
		return Vector2<T>(scalar - v.data[0], scalar - v.data[1]);
	}
	inline friend Vector2<T> operator*(T scalar, const Vector2<T> &v) {
		return Vector2<T>(scalar*v.data[0], scalar*v.data[1]);
	}
	inline friend Vector2<T> operator/ (T scalar, const Vector2<T> &v) {
		return Vector2<T>(scalar / v.data[0], scalar / v.data[1]);
	}

	inline Vector2<T>  operator+	(const T scalar) const {
		return Vector2<T>(data[0] + scalar, data[1] + scalar);
	}
	inline Vector2<T>  operator-	(const T scalar) const {
		return Vector2<T>(data[0] - scalar, data[1] - scalar);
	}
	inline Vector2<T>  operator*(const T scalar) const {
		return Vector2<T>(data[0]*scalar, data[1]*scalar);
	}
	inline Vector2<T>  operator/	(const T scalar) const {
		T invScalar = (T)1.0 / scalar;
		return Vector2<T>(data[0]*invScalar, data[1]*invScalar);
	}

	inline Vector2<T>& operator=	(const T scalar) {
		data[0] = data[1] = scalar;
		return *this;
	}
	inline Vector2<T>& operator+=	(const T scalar) {
		data[0] += scalar;
		data[1] += scalar;
		return *this;
	}
	inline Vector2<T>& operator-=	(const T scalar) {
		data[0] -= scalar;
		data[1] -= scalar;
		return *this;
	}
	inline Vector2<T>& operator*=	(const T scalar) {
		data[0] *= scalar;
		data[1] *= scalar;
		return *this;
	}
	inline Vector2<T>& operator/=	(const T scalar) {
		T invScalar = (T)1.0 / scalar;
		data[0] *= invScalar;
		data[1] *= invScalar;
		return *this;
	}

	inline Vector2<T>  operator+	(const Vector2<T>& vector) const {
		return Vector2<T>(data[0] + vector.data[0], data[1] + vector.data[1]);
	}
	inline Vector2<T>  operator-	(const Vector2<T>& vector) const {
		return Vector2<T>(data[0] - vector.data[0], data[1] - vector.data[1]);
	}
	inline Vector2<T>  operator*(const Vector2<T>& vector) const {
		return Vector2<T>(data[0]*vector.data[0], data[1]*vector.data[1]);
	}
	inline Vector2<T>  operator/	(const Vector2<T>& vector) const {
		return Vector2<T>(data[0] / vector.data[0], data[1] / vector.data[1]);
	}

	inline Vector2<T>& operator=	(const Vector2<T>& vector) {
		data[0] = vector.data[0];
		data[1] = vector.data[1];
		return *this;
	}
	inline Vector2<T>& operator+=	(const Vector2<T>& vector) {
		data[0] += vector.data[0];
		data[1] += vector.data[1];
		return *this;
	}
	inline Vector2<T>& operator-=	(const Vector2<T>& vector) {
		data[0] -= vector.data[0];
		data[1] -= vector.data[1];
		return *this;
	}
	inline Vector2<T>& operator*=	(const Vector2<T>& vector) {
		data[0] *= vector.data[0];
		data[1] *= vector.data[1];
		return *this;
	}
	inline Vector2<T>& operator/=	(const Vector2<T>& vector) {
		data[0] /= vector.data[0];
		data[1] /= vector.data[1];
		return *this;
	}

	inline bool operator==(const Vector2<T> &vector) const {
		return (abs(data[0] - vector.data[0]) < FLT_EPSILON && abs(data[1] - vector.data[1]) < FLT_EPSILON);
	}
	inline bool operator!=(const Vector2<T> &vector) const {
		return (abs(data[0] - vector.data[0]) > FLT_EPSILON || abs(data[1] - vector.data[1]) > FLT_EPSILON);
	}
	inline bool operator< (const Vector2<T> &vector) const {
		if (data[0] < vector.data[0]) return true;
		if (data[0] == vector.data[0]) if (data[1] < vector.data[1]) return true;
		return false;
	}
	inline bool operator<=(const Vector2<T> &vector) const {
		if (data[0] < vector.data[0]) return true;
		if (data[0] == vector.data[0]) if (data[1] <= vector.data[1]) return true;
		return false;
	}
	inline bool operator> (const Vector2<T> &vector) const {
		if (data[0] > vector.data[0]) return true;
		if (data[0] == vector.data[0]) if (data[1] > vector.data[1]) return true;
		return false;
	}
	inline bool operator>=(const Vector2<T> &vector) const {
		if (data[0] > vector.data[0]) return true;
		if (data[0] == vector.data[0]) if (data[1] >= vector.data[1]) return true;
		return false;
	}

	inline Vector2<T> maxVector(const Vector2<T> &vector) {
		return Vector2<T>(max(data[0], vector.data[0]), max(data[1], vector.data[1]));
	}
	inline Vector2<T> minVector(const Vector2<T> &vector) {
		return Vector2<T>(min(data[0], vector.data[0]), min(data[1], vector.data[1]));
	}
	inline T maxValue() {
		return max(data[0], data[1]);
	}
	inline T minValue() {
		return min(data[0], data[1]);
	}

	inline void print() {
		cout << data[0] << " " << data[1] << endl;
	}
	inline bool isZero() {
		return ((abs(data[0]) < FLT_EPSILON) && (abs(data[1]) < FLT_EPSILON));
	}
	inline void absoluteValue() {
		data[0] = abs(data[0]);
		data[1] = abs(data[1]);
	}
	inline void roundUp() {
		data[0] = ceil(data[0]);
		data[1] = ceil(data[1]);
	}
	inline void roundDown() {
		data[0] = floor(data[0]);
		data[1] = floor(data[1]);
	}
	inline void clamp(T minVal, T maxVal) {
		if (data[0] <= minVal) data[0] = minVal;
		else if (data[0] >= maxVal) data[0] = maxVal;
		if (data[1] <= minVal) data[1] = minVal;
		else if (data[1] >= maxVal) data[1] = maxVal;
	}
	inline void clamp(Vector2<T> minVal, Vector2<T> maxVal) {
		if (data[0] <= minVal.data[0]) data[0] = minVal.data[0];
		else if (data[0] >= maxVal.data[0]) data[0] = maxVal.data[0];
		if (data[1] <= minVal.data[1]) data[1] = minVal.data[1];
		else if (data[1] >= maxVal.data[1]) data[1] = maxVal.data[1];
	}

	inline float angleBetween(const Vector2<T> &vector, bool normalized) {
		float dotProduct;

		if (normalized) dotProduct = this->dot(vector);
		else {
			Vector2<T> b = vector;
			b.normalize();
			dotProduct = ((*this).unit()).dot(b);
		}

		if (dotProduct < -1) dotProduct = -1;
		else if (dotProduct > 1) dotProduct = 1;

		return acosf(dotProduct);
	}

	inline Vector2<T> perpendicular() {
		return Vector2<T>(-data[1], data[0]);
	}

	inline void generateRandomColor() {
		generateRandom(0, 1);
	}
	inline void generateRandom(const T &minVal, const T &maxVal) {
		generateRandom(Vector2<T>(minVal, minVal), Vector2<T>(maxVal, maxVal));
	}
	inline void generateRandom(const Vector2<T> &minVal, const Vector2<T> &maxVal) {
		Vector2<T> dim = maxVal - minVal;
		data[0] = (((T)rand()) / (T)RAND_MAX * dim.data[0] + minVal.data[0]);
		data[1] = (((T)rand()) / (T)RAND_MAX * dim.data[1] + minVal.data[1]);
	}

	inline double signedDistanceToLine(const Vector2<T> &linePt1, const Vector2<T> &linePt2) const {
		return (linePt1 -linePt2).perpendicular().dot(linePt1 - *this);
	}

	inline bool isInsideBox(const Vector2<T> &boxMin, const Vector2<T> &boxMax) const {
		return data[0] >= boxMin.data[0] && data[1] >= boxMin.data[1] && data[0] <= boxMax.data[0] && data[1] <= boxMax.data[1];
	}
	inline bool isInsideCircle(const Vector2<T> &circleCenter, const T &circleRadius) const {
		return (*this - circleCenter).magnitudeSquared() <= circleRadius*circleRadius;
	}
};


template <class T> class Vector3
{
public:
	typedef T	value_type;
	
	union {
		struct {
			T x, y, z;
		};
		struct {
			T r, g, b;
		};
		T data[3];
	};

	Vector3() {
		data[0] = data[1] = data[2] = 0;
	}
	Vector3<T>(const T &v0, const T &v1, const T &v2) {
		data[0] = (T)v0;
		data[1] = (T)v1;
		data[2] = (T)v2;
	}
	Vector3<T>(const Vector3<T>& copy) {
		data[0] = copy.data[0];
		data[1] = copy.data[1];
		data[2] = copy.data[2];
	}
	template <class S> Vector3<T>(const Vector3<S>& copy) {
		data[0] = (T)copy.data[0];
		data[1] = (T)copy.data[1];
		data[2] = (T)copy.data[2];
	}
	template <class S> Vector3<T>(S val[3]) {
		data[0] = (T)val[0];
		data[1] = (T)val[1];
		data[2] = (T)val[2];
	}

	inline unsigned int getDim() const {
		return 3;
	}

	inline T operator()(int pos) const {
		if (pos == 0) return data[0];
		if (pos == 1) return data[1];
		if (pos == 2) return data[2];
		fprintf(stderr, "Vector: Illegal data for () operator\n");
		return 0;
	}

	inline void set(const T &v0, const T &v1, const T &v2) {
		data[0] = v0;
		data[1] = v1;
		data[2] = v2;
	}
	inline void set(const T v[3]) {
		set(v[0], v[1], v[2]);
	}

	template <class S> inline void copyTo(S *v) const {
		v[0] = data[0];
		v[1] = data[1];
		v[2] = data[2];
	}

	inline double magnitude() const {
		return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
	}
	inline T magnitudeSquared() const {
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
	}
	inline Vector3<T> unit() const {
		Vector3<T> returnVal(data[0], data[1], data[2]);
		returnVal.normalize();
		return returnVal;
	}
	inline void normalize() {
		T magSqr = magnitudeSquared();
		if (magSqr != 0) {
			T magInv = (T)1.0 / sqrt(magSqr);
			*this *= magInv;
		}
	}

	inline T dot(const Vector3<T> &right) const {
		return (data[0]*right.data[0] + data[1]*right.data[1] + data[2]*right.data[2]);
	}
	inline Vector3<T> cross(const Vector3<T> &right) const {
		return Vector3<T>(data[1]*right.data[2] - data[2]*right.data[1], data[2]*right.data[0] - data[0]*right.data[2], data[0]*right.data[1] - data[1]*right.data[0]);
	}

	inline Vector3<T>  operator-	() const {
		return Vector3<T>(-data[0], -data[1], -data[2]);
	}

	inline friend Vector3<T> operator+ (const T &scalar, const Vector3<T> &v) {
		return Vector3<T>(scalar + v.data[0], scalar + v.data[1], scalar + v.data[2]);
	}
	inline friend Vector3<T> operator- (const T &scalar, const Vector3<T> &v) {
		return Vector3<T>(scalar - v.data[0], scalar - v.data[1], scalar - v.data[2]);
	}
	inline friend Vector3<T> operator*(const T &scalar, const Vector3<T> &v) {
		return Vector3<T>(scalar*v.data[0], scalar*v.data[1], scalar*v.data[2]);
	}
	//template <class S, class V> inline friend Vector3<T> operator* (const S &scalar, const Vector3<V> &v) { return Vector3<T>((T)scalar*v.data[0], (T)scalar*v.data[1], (T)scalar*v.data[2]); }
	inline friend Vector3<T> operator/ (const T &scalar, const Vector3<T> &v) {
		return Vector3<T>(scalar / v.data[0], scalar / v.data[1], scalar / v.data[2]);
	}

	inline Vector3<T>  operator+ (const T &scalar) const {
		return Vector3<T>(data[0] + scalar, data[1] + scalar, data[2] + scalar);
	}
	inline Vector3<T>  operator- (const T &scalar) const {
		return Vector3<T>(data[0] - scalar, data[1] - scalar, data[2] - scalar);
	}
	inline Vector3<T>  operator*(const T &scalar) const {
		return Vector3<T>(data[0]*scalar, data[1]*scalar, data[2]*scalar);
	}
	//template <class S> inline Vector3<T>  operator* (const S &scalar) const { return Vector3<T>(data[0]*(T)scalar, data[1]*(T)scalar, data[2]*(T)scalar); }
	inline Vector3<T>  operator/ (const T &scalar) const {
		T invScalar = (T)1.0 / scalar;
		return Vector3<T>(data[0]*invScalar, data[1]*invScalar, data[2]*invScalar);
	}

	inline Vector3<T>& operator= (const T &scalar) {
		data[0] = data[1] = data[2] = scalar;
		return *this;
	}
	inline Vector3<T>& operator+=(const T &scalar) {
		data[0] += scalar;
		data[1] += scalar;
		data[2] += scalar;
		return *this;
	}
	inline Vector3<T>& operator-=(const T &scalar) {
		data[0] -= scalar;
		data[1] -= scalar;
		data[2] -= scalar;
		return *this;
	}
	inline Vector3<T>& operator*=(const T &scalar) {
		data[0] *= scalar;
		data[1] *= scalar;
		data[2] *= scalar;
		return *this;
	}
	inline Vector3<T>& operator/=(const T &scalar) {
		T invScalar = (T)1.0 / scalar;
		data[0] *= invScalar;
		data[1] *= invScalar;
		data[2] *= invScalar;
		return *this;
	}

	inline Vector3<T>  operator+ (const Vector3<T> &vector) const {
		return Vector3<T>(data[0] + vector.data[0], data[1] + vector.data[1], data[2] + vector.data[2]);
	}
	inline Vector3<T>  operator- (const Vector3<T> &vector) const {
		return Vector3<T>(data[0] - vector.data[0], data[1] - vector.data[1], data[2] - vector.data[2]);
	}
	inline Vector3<T>  operator*(const Vector3<T> &vector) const {
		return Vector3<T>(data[0]*vector.data[0], data[1]*vector.data[1], data[2]*vector.data[2]);
	}
	inline Vector3<T>  operator/ (const Vector3<T> &vector) const {
		return Vector3<T>(data[0] / vector.data[0], data[1] / vector.data[1], data[2] / vector.data[2]);
	}

	inline Vector3<T>& operator= (const Vector3<T> &vector) {
		data[0] =  vector.data[0];
		data[1] =  vector.data[1];
		data[2] =  vector.data[2];
		return *this;
	}
	inline Vector3<T>& operator+=(const Vector3<T> &vector) {
		data[0] += vector.data[0];
		data[1] += vector.data[1];
		data[2] += vector.data[2];
		return *this;
	}
	inline Vector3<T>& operator-=(const Vector3<T> &vector) {
		data[0] -= vector.data[0];
		data[1] -= vector.data[1];
		data[2] -= vector.data[2];
		return *this;
	}
	inline Vector3<T>& operator*=(const Vector3<T> &vector) {
		data[0] *= vector.data[0];
		data[1] *= vector.data[1];
		data[2] *= vector.data[2];
		return *this;
	}
	inline Vector3<T>& operator/=(const Vector3<T> &vector) {
		data[0] /= vector.data[0];
		data[1] /= vector.data[1];
		data[2] /= vector.data[2];
		return *this;
	}

	//inline bool operator==(const Vector3<T> &vector) const { return (abs(data[0]-vector.data[0]) < FLT_EPSILON && abs(data[1]-vector.data[1]) < FLT_EPSILON && abs(data[2]-vector.data[2]) < FLT_EPSILON); }
	inline bool operator==(const Vector3<T> &vector) const {
		return data[0] == vector.data[0] && data[1] == vector.data[1] && data[2] == vector.data[2];
	}
	//inline bool operator!=(const Vector3<T> &vector) const { return (abs(data[0]-vector.data[0]) > FLT_EPSILON || abs(data[1]-vector.data[1]) > FLT_EPSILON || abs(data[2]-vector.data[2]) > FLT_EPSILON); }
	inline bool operator!=(const Vector3<T> &vector) const {
		return data[0] != vector.data[0] || data[1] != vector.data[1] || data[2] != vector.data[2];
	}
	inline bool operator< (const Vector3<T> &vector) const {
		if (data[0] < vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] < vector.data[1]) return true;
			if (data[1] == vector.data[1]) if (data[2] < vector.data[2]) return true;
		}
		return false;
	}
	inline bool operator<=(const Vector3<T> &vector) const {
		if (data[0] < vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] < vector.data[1]) return true;
			if (data[1] == vector.data[1]) if (data[2] <= vector.data[2]) return true;
		}
		return false;
	}
	inline bool operator> (const Vector3<T> &vector) const {
		if (data[0] > vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] > vector.data[1]) return true;
			if (data[1] == vector.data[1]) if (data[2] > vector.data[2]) return true;
		}
		return false;
	}
	inline bool operator>=(const Vector3<T> &vector) const {
		if (data[0] > vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] > vector.data[1]) return true;
			if (data[1] == vector.data[1]) if (data[2] >= vector.data[2]) return true;
		}
		return false;
	}

	inline bool operator==(const T &scalar) const {
		return (abs(data[0] - scalar) < FLT_EPSILON && abs(data[1] - scalar) < FLT_EPSILON && abs(data[2] - scalar) < FLT_EPSILON);
	}
	inline bool operator!=(const T &scalar) const {
		return (abs(data[0] - scalar) > FLT_EPSILON || abs(data[1] - scalar) > FLT_EPSILON || abs(data[2] - scalar) > FLT_EPSILON);
	}
	inline bool operator< (const T &scalar) const {
		return data[0] <  scalar && data[1] <  scalar && data[2] <  scalar;
	}
	inline bool operator<=(const T &scalar) const {
		return data[0] <= scalar && data[1] <= scalar && data[2] <= scalar;
	}
	inline bool operator> (const T &scalar) const {
		return data[0] >  scalar && data[1] >  scalar && data[2] >  scalar;
	}
	inline bool operator>=(const T &scalar) const {
		return data[0] >= scalar && data[1] >= scalar && data[2] >= scalar;
	}

	inline Vector3<T> maxVector(const Vector3<T> &vector) const {
		return Vector3<T>(max(data[0], vector.data[0]), max(data[1], vector.data[1]), max(data[2], vector.data[2]));
	}
	inline Vector3<T> maxVector(T v0, T v1, T v2) const {
		return Vector3<T>(max(data[0], v0), max(data[1], v1), max(data[2], v2));
	}
	inline Vector3<T> minVector(const Vector3<T> &vector) const {
		return Vector3<T>(min(data[0], vector.data[0]), min(data[1], vector.data[1]), min(data[2], vector.data[2]));
	}
	inline Vector3<T> minVector(T v0, T v1, T v2) const {
		return Vector3<T>(min(data[0], v0), min(data[1], v1), min(data[2], v2));
	}
	inline T maxValue() const {
		return max(max(data[0], data[1]), data[2]);
	}
	inline T minValue() const {
		return min(min(data[0], data[1]), data[2]);
	}
	inline void updateMinMax(Vector3<T> &minVec, Vector3<T> &maxVec) const {
		if (data[0] < minVec.data[0]) minVec.data[0] = data[0];
		else if (data[0] > maxVec.data[0]) maxVec.data[0] = data[0];

		if (data[1] < minVec.data[1]) minVec.data[1] = data[1];
		else if (data[1] > maxVec.data[1]) maxVec.data[1] = data[1];

		if (data[2] < minVec.data[2]) minVec.data[2] = data[2];
		else if (data[2] > maxVec.data[2]) maxVec.data[2] = data[2];
	}

	inline void print() const {
		cout << data[0] << " " << data[1] << " " << data[2] << endl;
	}
	inline bool isZero() const {
		return ((abs(data[0]) < FLT_EPSILON) && (abs(data[1]) < FLT_EPSILON) && (abs(data[2]) < FLT_EPSILON));
	}
	inline void absoluteValue() {
		data[0] = abs(data[0]);
		data[1] = abs(data[1]);
		data[2] = abs(data[2]);
	}
	inline void roundUp() {
		data[0] = ceil(data[0]);
		data[1] = ceil(data[1]);
		data[2] = ceil(data[2]);
	}
	inline void roundDown() {
		data[0] = floor(data[0]);
		data[1] = floor(data[1]);
		data[2] = floor(data[2]);
	}
	inline void clamp(const T &minVal, const T &maxVal) {
		if (data[0] <= minVal) data[0] = minVal;
		else if (data[0] >= maxVal) data[0] = maxVal;
		if (data[1] <= minVal) data[1] = minVal;
		else if (data[1] >= maxVal) data[1] = maxVal;
		if (data[2] <= minVal) data[2] = minVal;
		else if (data[2] >= maxVal) data[2] = maxVal;
	}
	inline void clamp(const Vector3<T> &minVal, const Vector3<T> &maxVal) {
		if (data[0] <= minVal.data[0]) data[0] = minVal.data[0];
		else if (data[0] >= maxVal.data[0]) data[0] = maxVal.data[0];
		if (data[1] <= minVal.data[1]) data[1] = minVal.data[1];
		else if (data[1] >= maxVal.data[1]) data[1] = maxVal.data[1];
		if (data[2] <= minVal.data[2]) data[2] = minVal.data[2];
		else if (data[2] >= maxVal.data[2]) data[2] = maxVal.data[2];
	}

	inline Vector3<T> rotateX(const double &radians) {
		if (radians == 0) return *this;
		double sinAngle = sin(radians);
		double cosAngle = cos(radians);
		return Vector3<T>(data[0], (T)(data[1]*cosAngle - data[2]*sinAngle), (T)(data[1]*sinAngle + data[2]*cosAngle));
	}
	inline Vector3<T> rotateY(const double &radians) {
		if (radians == 0) return *this;
		double sinAngle = sin(radians);
		double cosAngle = cos(radians);
		return Vector3<T>((T)(data[0]*cosAngle + data[2]*sinAngle), data[1], -(T)(data[0]*sinAngle + data[2]*cosAngle));
	}
	inline Vector3<T> rotateZ(const double &radians) {
		if (radians == 0) return *this;
		double sinAngle = sin(radians);
		double cosAngle = cos(radians);
		return Vector3<T>((T)(data[0]*cosAngle - data[1]*sinAngle), (T)(data[0]*sinAngle + data[1]*cosAngle), data[2]);
	}
	template <class S> inline Vector3<T> rotate(const double &radians, const Vector3<S> &normalizedAxis) {
		if (radians == 0) return *this;
		T ux = (T)normalizedAxis.data[0] * data[0];
		T vy = (T)normalizedAxis.data[1] * data[1];
		T wz = (T)normalizedAxis.data[2] * data[2];
		double sinAngle = sin(radians);
		double cosAngle = cos(radians);
		T uu = (T)(normalizedAxis.data[0] * normalizedAxis.data[0]);
		T vv = (T)(normalizedAxis.data[1] * normalizedAxis.data[1]);
		T ww = (T)(normalizedAxis.data[2] * normalizedAxis.data[2]);
		T dotProduct = ux + vy + wz;
		return Vector3<T>((T)(normalizedAxis.data[0]*dotProduct + cosAngle*(data[0]*(vv + ww) - normalizedAxis.data[0]*(vy + wz)) +
		                      sinAngle*(normalizedAxis.data[1]*data[2] - normalizedAxis.data[2]*data[1])),
		                  (T)(normalizedAxis.data[1]*dotProduct + cosAngle*(data[1]*(uu + ww) - normalizedAxis.data[1]*(ux + wz)) +
		                      sinAngle*(normalizedAxis.data[2]*data[0] - normalizedAxis.data[0]*data[2])),
		                  (T)(normalizedAxis.data[2]*dotProduct + cosAngle*(data[2]*(uu + vv) - normalizedAxis.data[2]*(ux + vy)) +
		                      sinAngle*(normalizedAxis.data[0]*data[1] - normalizedAxis.data[1]*data[0])));
	}

	template <class S> inline double angleBetween(const Vector3<S> &vector, bool normalized) const {
		double dotProduct;
		if (normalized) dotProduct = this->dot(vector);
		else dotProduct = (this->unit()).dot(vector.unit());

		if (dotProduct < -1) dotProduct = -1;
		else if (dotProduct > 1) dotProduct = 1;

		return acos(dotProduct);
	}

	/*inline Vector3<T> getPerpendicular(bool normalized) const {
		Vector3<T> absolute(abs(data[0]), abs(data[1]), abs(data[2]));

		int min, other1, other2;
		if (absolute.data[0] <= absolute.data[1]) {
			if (absolute.data[0] <= absolute.data[2]) { min = 0; other1 = 1; other2 = 2; }
			else { min = 2; other1 = 0; other2 = 1; }
		}
		else {
			if (absolute.data[1] <= absolute.data[2]) { min = 1; other1 = 0; other2 = 2; }
			else { min = 2; other1 = 0; other2 = 1; }
		}

		Vector3<T> result;
		result.data[min]    = 0.0;
		result.data[other1] = data[other2];
		result.data[other2] = -data[other1];

		if (!normalized && data[min] != 0) result.normalize();
		return result;
	}*/

	// produces a perpendicular vector to the one given, output not necessarily unit length or the same length as input
	inline Vector3<T> getPerpendicular() const {
		int idx = 0;
		if (data[0]*data[0] < data[1]*data[1]) idx = 1;
		if (data[idx]*data[idx] < data[2]*data[2]) idx = 2;
		switch (idx) {
		case 0:
			return Vector3<T>(data[1] - data[2],	-data[0],			data[0]);
		case 1:
			return Vector3<T>(-data[1],			data[0] - data[2],	data[1]);
		case 2:
			return Vector3<T>(-data[2],			data[2],			data[0] - data[1]);
		}
		return *this;
	}

	inline Vector3<T> getSorted(const bool sortMinToMax = true) const {
		if (sortMinToMax) {
			if (data[0] < data[1]) {
				if (data[1] < data[2]) return *this;
				else if (data[0] < data[2]) return Vector3<T>(data[0], data[2], data[1]);
				else return Vector3<T>(data[2], data[0], data[1]);
			}
			else if (data[0] < data[2]) return Vector3<T>(data[1], data[0], data[2]);
			else if (data[1] < data[2]) return Vector3<T>(data[1], data[2], data[0]);
			else return Vector3<T>(data[2], data[1], data[0]);
		}
		else {
			if (data[0] > data[1]) {
				if (data[1] > data[2]) return *this;
				else if (data[0] > data[2]) return Vector3<T>(data[0], data[2], data[1]);
				else return Vector3<T>(data[2], data[0], data[1]);
			}
			else if (data[0] > data[2]) return Vector3<T>(data[1], data[0], data[2]);
			else if (data[1] > data[2]) return Vector3<T>(data[1], data[2], data[0]);
			else return Vector3<T>(data[2], data[1], data[0]);
		}
	}

	inline unsigned int getIndexOfMaxValue() const {
		return (data[0] > data[1] ? (data[1] > data[2] ? 0 : (data[0] > data[2] ? 0 : 2)) : (data[2] > data[1] ? 2 : 1));
	}
	inline unsigned int getIndexOfMinValue() const {
		return (data[0] < data[1] ? (data[1] < data[2] ? 0 : (data[0] < data[2] ? 0 : 2)) : (data[2] < data[1] ? 2 : 1));
	}

	inline Vector3<T> projectToRay(const Vector3<T> &rayPt, const Vector3<T> &rayDir) const {
		return projectToLine(rayPt -rayDir*(T)100000, rayPt + rayDir*(T)100000);
	}
	inline Vector3<T> projectToLine(const Vector3<T> &linePt1, const Vector3<T> &linePt2) const {
		Vector3<T> e1 = *this - linePt1;
		Vector3<T> e2 = linePt2 - linePt1;
		T t = e1.dot(e2);
		T mag = e2.magnitudeSquared();
		if (mag != 0) t /= mag;
		return linePt1 + t*e2;
	}

	inline void generateRandomColor() {
		generateRandom(0, 1);
	}
	inline void generateRandom(const T &minVal, const T &maxVal) {
		generateRandom(Vector3<T>(minVal, minVal, minVal), Vector3<T>(maxVal, maxVal, maxVal));
	}
	inline void generateRandom(const Vector3<T> &minVal, const Vector3<T> &maxVal) {
		Vector3<T> dim = maxVal - minVal;
		data[0] = (((T)rand()) / (T)RAND_MAX * dim.data[0] + minVal.data[0]);
		data[1] = (((T)rand()) / (T)RAND_MAX * dim.data[1] + minVal.data[1]);
		data[2] = (((T)rand()) / (T)RAND_MAX * dim.data[2] + minVal.data[2]);
	}

	inline Vector3<unsigned char> convertFloatColorToUnsignedChar() const {
		Vector3<T> scaled((T)(data[0]*255.0), (T)(data[1]*255.0), (T)(data[2]*255.0));
		scaled.clamp(0, 255);
		return Vector3<unsigned char>((unsigned char)scaled.x, (unsigned char)scaled.y, (unsigned char)scaled.z);
	}

	inline bool isInsideBox(const Vector3<T> &boxMin, const Vector3<T> &boxMax) const {
		return data[0] >= boxMin.data[0] && data[1] >= boxMin.data[1] && data[2] >= boxMin.data[2] &&
		       data[0] <= boxMax.data[0] && data[1] <= boxMax.data[1] && data[2] <= boxMax.data[2];
	}
	inline bool isInsideSphere(const Vector3<T> &sphereCenter, const T &sphereRadius) const {
		return (*this - sphereCenter).magnitudeSquared() <= sphereRadius*sphereRadius;
	}

	// the two points define the center of end caps
	inline bool isInsideCylinder(const Vector3<T> &cylinderPoint1, const Vector3<T> &cylinderPoint2, const T &cylinderRadius) const {
		Vector3<T> dc = cylinderPoint2 - cylinderPoint1;	// translate points to origin
		Vector3<T> dx = *this - cylinderPoint1;

		T dotProduct = dc.dot(dx);
		if (dotProduct < 0) return false;  // if dotProduct < 0, point is behind pt1 cap

		T dcMagSq = dc.magnitudeSquared();
		if (dotProduct > dcMagSq) return false;  // if dotProduct > cylinder's axis length then point is outside pt2 cap

		// find distance sq. from point to line = (1 - cos^2(dx to d)) * |dx|^2 = (1 - (dx * d)^2 / (|dx|^2 * |d|^2) ) * |dx|^2 = dx*dx - dot*dot / cylinderLengthSq
		if (dx.magnitudeSquared() - dotProduct*dotProduct / dcMagSq > cylinderRadius*cylinderRadius) return false;
		return true;  // if needed in the future, the distance sq. to cylinder boundary is given by the left hand side of the above if-statement
	}

	// works with irregular pyramids
	inline bool isInsideTetrahedron(const Vector3<T> &pyramidPoint1, const Vector3<T> &pyramidPoint2, const Vector3<T> &pyramidPoint3, const Vector3<T> &pyramidPoint4) const {
		const T eps = (T)1.0e-9;
		Vector3<T> a1 = pyramidPoint1 - pyramidPoint4, a2 = pyramidPoint2 - pyramidPoint4, a3 = pyramidPoint3 - pyramidPoint4;
		Vector3<T> B = *this - pyramidPoint4;  // "a"s and B represent matrix to be inverted
		Vector3<T> A = a2.cross(a3);  // 1st column of inverse matrix

		T denominator = A.dot(a1), bary = A.dot(B); // = bary.x * denominator
		T alpha = -eps * denominator, dMinusB = denominator - bary;

		if (denominator > 0 ? (bary > alpha && dMinusB > alpha) : (bary < alpha && dMinusB < alpha)) {
			bary = -(a1.cross(a3)).dot(B); // = bary.y * denominator
			dMinusB -= bary;

			if (denominator > 0 ? (bary > alpha && dMinusB > alpha) : (bary < alpha && dMinusB < alpha)) {
				bary = (a1.cross(a2)).dot(B); // = bary.z * denominator
				dMinusB -= bary;

				if (denominator > 0 ? (bary > alpha && dMinusB > alpha) : (bary < alpha && dMinusB < alpha)) return true;
			}
		}
		return false;
	}

	inline double signedDistanceToSphere(const Vector3<T> &sphereCenter, const T sphereRadius) const {
		return (*this - sphereCenter).magnitude() - sphereRadius;
	}

	inline double signedDistanceToBox(const Vector3<T> &boxCenter, const Vector3<T> &boxHalfSize) const {
		Vector3<double> diff = *this - boxCenter;
		diff.absoluteValue();
		diff -= boxHalfSize;
		if (diff.x > 0) {
			if (diff.y > 0) {
				if (diff.z > 0) return diff.magnitude();
				else return Vector2<double>(diff.x, diff.y).magnitude();
			}
			else if (diff.z > 0) return Vector2<double>(diff.x, diff.z).magnitude();
			else return diff.x;
		}
		else if (diff.y > 0) {
			if (diff.z > 0) return Vector2<double>(diff.y, diff.z).magnitude();
			else return diff.y;
		}
		else if (diff.z > 0) return diff.z;
		else return diff.maxValue(); // all vals are neg. so use max to get min dist.
	}

	// points should be in counter-clockwise order to obtain the correct normal
	inline double signedDistanceToPlane(const Vector3<T> &pointOnPlane1, const Vector3<T> &pointOnPlane2, const Vector3<T> &pointOnPlane3) const {
		return signedDistanceToPlane(pointOnPlane1, ((pointOnPlane2 - pointOnPlane1).cross(pointOnPlane3 - pointOnPlane1)).unit());
	}
	inline double signedDistanceToPlane(const Vector3<T> &pointOnPlane, const Vector3<T> &planeUnitNormal) const {
		return (*this - pointOnPlane).dot(planeUnitNormal);
	}
	inline double angleBetweenPlane(const Vector3<T> &pointOnPlane1, const Vector3<T> &pointOnPlane2, const Vector3<T> &pointOnPlane3, bool normalized) const {
		if (normalized) return this->angleBetween(((pointOnPlane2 - pointOnPlane1).cross(pointOnPlane3 - pointOnPlane1)).unit(), true);
		else return this->angleBetween((pointOnPlane2 - pointOnPlane1).cross(pointOnPlane3 - pointOnPlane1), false);
	}


	inline int getRelativeDirectionTo(const Vector3<T> &pt) const {
		if (x <= pt.x) {
			if (y <= pt.y) {
				if (z <= pt.z) return CORNER_LBB;
				else return CORNER_LBF;
			}
			else           {
				if (z <= pt.z) return CORNER_LTB;
				else return CORNER_LTF;
			}
		}
		else {
			if (y <= pt.y) {
				if (z <= pt.z) return CORNER_RBB;
				else return CORNER_RBF;
			}
			else           {
				if (z <= pt.z) return CORNER_RTB;
				else return CORNER_RTF;
			}
		}
	}
};

template <class T> class Vector4
{
public:
	typedef T	value_type;
	
	union {
		struct {
			T x, y, z, w;
		};
		struct {
			T r, g, b, a;
		};
		T data[4];
	};

	Vector4() {
		data[0] = data[1] = data[2] = data[3] = 0;
	}
	Vector4(const T v0, const T v1, const T v2, const T v3) {
		data[0] = v0;
		data[1] = v1;
		data[2] = v2;
		data[3] = v3;
	}
	Vector4(const Vector4<T>& copy) {
		data[0] = copy.data[0];
		data[1] = copy.data[1];
		data[2] = copy.data[2];
		data[3] = copy.data[3];
	}
	Vector4(T val[4]) {
		data[0] = val[0];
		data[1] = val[1];
		data[2] = val[2];
		data[3] = val[3];
	}

	inline unsigned int getDim() const {
		return 4;
	}

	inline T operator()(int pos) const {
		if (pos == 0) return data[0];
		if (pos == 1) return data[1];
		if (pos == 2) return data[2];
		if (pos == 3) return data[3];
		fprintf(stderr, "Vector: Illegal data for () operator\n");
		return 0;
	}

	inline void set(T v0, T v1, T v2, T v3) {
		data[0] = v0;
		data[1] = v1;
		data[2] = v2;
		data[3] = v3;
	}

	inline T magnitude() {
		return T(sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3]));
	}
	inline T magnitudeSquared() {
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3];
	}
	inline Vector4<T> unit() const {
		Vector4<T> returnVal(data[0], data[1], data[2], data[3]);
		returnVal.normalize();
		return returnVal;
	}
	inline void normalize() {
		T magSqr = magnitudeSquared();
		if (magSqr != 0) {
			T magInv = (T)1.0 / sqrt(magSqr);
			*this *= magInv;
		}
	}

	inline T dot(const Vector4<T> &right) {
		return (data[0]*right.data[0] + data[1]*right.data[1] + data[2]*right.data[2] + data[3]*right.data[3]);
	}

	inline Vector4<T>  operator-	() const {
		return Vector4<T>(-data[0], -data[1], -data[2], -data[3]);
	}

	inline friend Vector4<T> operator+ (T scalar, const Vector4<T> &v) {
		return Vector4<T>(scalar + v.data[0], scalar + v.data[1], scalar + v.data[2], scalar + v.data[3]);
	}
	inline friend Vector4<T> operator- (T scalar, const Vector4<T> &v) {
		return Vector4<T>(scalar - v.data[0], scalar - v.data[1], scalar - v.data[2], scalar - v.data[3]);
	}
	inline friend Vector4<T> operator*(T scalar, const Vector4<T> &v) {
		return Vector4<T>(scalar*v.data[0], scalar*v.data[1], scalar*v.data[2], scalar*v.data[3]);
	}
	inline friend Vector4<T> operator/ (T scalar, const Vector4<T> &v) {
		return Vector4<T>(scalar / v.data[0], scalar / v.data[1], scalar / v.data[2], scalar / v.data[3]);
	}

	inline Vector4<T>  operator+ (const T scalar) const {
		return Vector4<T>(data[0] + scalar, data[1] + scalar, data[2] + scalar, data[3] + scalar);
	}
	inline Vector4<T>  operator- (const T scalar) const {
		return Vector4<T>(data[0] - scalar, data[1] - scalar, data[2] - scalar, data[3] - scalar);
	}
	inline Vector4<T>  operator*(const T scalar) const {
		return Vector4<T>(data[0]*scalar, data[1]*scalar, data[2]*scalar, data[3]*scalar);
	}
	inline Vector4<T>  operator/ (const T scalar) const {
		T invScalar = (T)1.0 / scalar;
		return Vector4<T>(data[0]*invScalar, data[1]*invScalar, data[2]*invScalar, data[3]*invScalar);
	}

	inline Vector4<T>& operator= (const T scalar) {
		data[0] = data[1] = data[2] = data[3] = scalar;
		return *this;
	}
	inline Vector4<T>& operator+=(const T scalar) {
		data[0] += scalar;
		data[1] += scalar;
		data[2] += scalar;
		data[3] += scalar;
		return *this;
	}
	inline Vector4<T>& operator-=(const T scalar) {
		data[0] -= scalar;
		data[1] -= scalar;
		data[2] -= scalar;
		data[3] -= scalar;
		return *this;
	}
	inline Vector4<T>& operator*=(const T scalar) {
		data[0] *= scalar;
		data[1] *= scalar;
		data[2] *= scalar;
		data[3] *= scalar;
		return *this;
	}
	inline Vector4<T>& operator/=(const T scalar) {
		T invScalar = (T)1.0 / scalar;
		data[0] *= invScalar;
		data[1] *= invScalar;
		data[2] *= invScalar;
		data[3] *= invScalar;
		return *this;
	}

	inline Vector4<T>  operator+ (const Vector4<T>& vector) const {
		return Vector4<T>(data[0] + vector.data[0], data[1] + vector.data[1], data[2] + vector.data[2], data[3] + vector.data[3]);
	}
	inline Vector4<T>  operator- (const Vector4<T>& vector) const {
		return Vector4<T>(data[0] - vector.data[0], data[1] - vector.data[1], data[2] - vector.data[2], data[3] - vector.data[3]);
	}
	inline Vector4<T>  operator*(const Vector4<T>& vector) const {
		return Vector4<T>(data[0]*vector.data[0], data[1]*vector.data[1], data[2]*vector.data[2], data[3]*vector.data[3]);
	}
	inline Vector4<T>  operator/ (const Vector4<T>& vector) const {
		return Vector4<T>(data[0] / vector.data[0], data[1] / vector.data[1], data[2] / vector.data[2], data[3] / vector.data[3]);
	}

	inline Vector4<T>& operator= (const Vector4<T>& vector) {
		data[0] = vector.data[0];
		data[1] = vector.data[1];
		data[2] = vector.data[2];
		data[3] = vector.data[3];
		return *this;
	}
	inline Vector4<T>& operator+=(const Vector4<T>& vector) {
		data[0] += vector.data[0];
		data[1] += vector.data[1];
		data[2] += vector.data[2];
		data[3] += vector.data[3];
		return *this;
	}
	inline Vector4<T>& operator-=(const Vector4<T>& vector) {
		data[0] -= vector.data[0];
		data[1] -= vector.data[1];
		data[2] -= vector.data[2];
		data[3] -= vector.data[3];
		return *this;
	}
	inline Vector4<T>& operator*=(const Vector4<T>& vector) {
		data[0] *= vector.data[0];
		data[1] *= vector.data[1];
		data[2] *= vector.data[2];
		data[3] *= vector.data[3];
		return *this;
	}
	inline Vector4<T>& operator/=(const Vector4<T>& vector) {
		data[0] /= vector.data[0];
		data[1] /= vector.data[1];
		data[2] /= vector.data[2];
		data[3] /= vector.data[3];
		return *this;
	}

	//inline bool operator==(const Vector4<T> &vector) const { return (abs(data[0]-vector.data[0]) < FLT_EPSILON && abs(data[1]-vector.data[1]) < FLT_EPSILON && abs(data[2]-vector.data[2]) < FLT_EPSILON && abs(data[3]-vector.data[3]) < FLT_EPSILON); }
	inline bool operator==(const Vector4<T> &vector) const {
		return data[0] == vector.data[0] && data[1] == vector.data[1] && data[2] == vector.data[2] && data[3] == vector.data[3];
	}
	//inline bool operator!=(const Vector4<T> &vector) const { return (abs(data[0]-vector.data[0]) > FLT_EPSILON || abs(data[1]-vector.data[1]) > FLT_EPSILON || abs(data[2]-vector.data[2]) > FLT_EPSILON || abs(data[3]-vector.data[3]) > FLT_EPSILON); }
	inline bool operator!=(const Vector4<T> &vector) const {
		return data[0] != vector.data[0] || data[1] != vector.data[1] || data[2] != vector.data[2] || data[3] != vector.data[3];
	}
	inline bool operator< (const Vector4<T> &vector) const {
		if (data[0] < vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] < vector.data[1]) return true;
			if (data[1] == vector.data[1]) {
				if (data[2] < vector.data[2]) return true;
				if (data[2] == vector.data[2]) if (data[3] < vector.data[3]) return true;
			}
		}
		return false;
	}
	inline bool operator<=(const Vector4<T> &vector) const {
		if (data[0] < vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] < vector.data[1]) return true;
			if (data[1] == vector.data[1]) {
				if (data[2] < vector.data[2]) return true;
				if (data[2] == vector.data[2]) if (data[3] <= vector.data[3]) return true;
			}
		}
		return false;
	}
	inline bool operator> (const Vector4<T> &vector) const {
		if (data[0] > vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] > vector.data[1]) return true;
			if (data[1] == vector.data[1]) {
				if (data[2] < vector.data[2]) return true;
				if (data[2] == vector.data[2]) if (data[3] > vector.data[3]) return true;
			}
		}
		return false;
	}
	inline bool operator>=(const Vector4<T> &vector) const {
		if (data[0] > vector.data[0]) return true;
		if (data[0] == vector.data[0]) {
			if (data[1] > vector.data[1]) return true;
			if (data[1] == vector.data[1]) {
				if (data[2] < vector.data[2]) return true;
				if (data[2] == vector.data[2]) if (data[3] >= vector.data[3]) return true;
			}
		}
		return false;
	}

	inline Vector4<T> maxVector(const Vector4<T> &vector) {
		return Vector4<T>(max(data[0], vector.data[0]), max(data[1], vector.data[1]), max(data[2], vector.data[2]), max(data[3], vector.data[3]));
	}
	inline Vector4<T> minVector(const Vector4<T> &vector) {
		return Vector4<T>(min(data[0], vector.data[0]), min(data[1], vector.data[1]), min(data[2], vector.data[2]), min(data[3], vector.data[3]));
	}
	inline T maxValue() {
		return max(max(data[0], data[1]), max(data[2], data[3]));
	}
	inline T minValue() {
		return min(min(data[0], data[1]), max(data[2], data[3]));
	}

	inline void print() const {
		cout << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << endl;
	}
	inline bool isZero() {
		return ((abs(data[0]) < FLT_EPSILON) && (abs(data[1]) < FLT_EPSILON) && (abs(data[2]) < FLT_EPSILON) && (abs(data[3]) < FLT_EPSILON));
	}
	inline void absoluteValue() {
		data[0] = abs(data[0]);
		data[1] = abs(data[1]);
		data[2] = abs(data[2]);
		data[3] = abs(data[3]);
	}
	inline void roundUp() {
		data[0] = ceil(data[0]);
		data[1] = ceil(data[1]);
		data[2] = ceil(data[2]);
		data[3] = ceil(data[3]);
	}
	inline void roundDown() {
		data[0] = floor(data[0]);
		data[1] = floor(data[1]);
		data[2] = floor(data[2]);
		data[3] = floor(data[3]);
	}
	inline void clamp(T minVal, T maxVal) {
		if (data[0] <= minVal) data[0] = minVal;
		else if (data[0] >= maxVal) data[0] = maxVal;
		if (data[1] <= minVal) data[1] = minVal;
		else if (data[1] >= maxVal) data[1] = maxVal;
		if (data[2] <= minVal) data[2] = minVal;
		else if (data[2] >= maxVal) data[2] = maxVal;
		if (data[3] <= minVal) data[3] = minVal;
		else if (data[3] >= maxVal) data[3] = maxVal;
	}
	inline void clamp(Vector3<T> minVal, Vector3<T> maxVal) {
		if (data[0] <= minVal.data[0]) data[0] = minVal.data[0];
		else if (data[0] >= maxVal.data[0]) data[0] = maxVal.data[0];
		if (data[1] <= minVal.data[1]) data[1] = minVal.data[1];
		else if (data[1] >= maxVal.data[1]) data[1] = maxVal.data[1];
		if (data[2] <= minVal.data[2]) data[2] = minVal.data[2];
		else if (data[2] >= maxVal.data[2]) data[2] = maxVal.data[2];
		if (data[2] <= minVal.data[3]) data[2] = minVal.data[3];
		else if (data[3] >= maxVal.data[3]) data[3] = maxVal.data[3];
	}

	inline float angleBetween(const Vector4<T> &vector, bool normalized) {
		float dotProduct;

		if (normalized) dotProduct = this->dot(vector);
		else {
			Vector4<T> b = vector;
			b.normalize();
			dotProduct = ((*this).unit()).dot(b);
		}

		if (dotProduct < -1) dotProduct = -1;
		else if (dotProduct > 1) dotProduct = 1;

		return acosf(dotProduct);
	}

	inline void compositeColorUsingOverOperator(const Vector4<T> &nextColor) {
		T scale = nextColor.a * (1.0f - a);
		r += nextColor.r * scale;
		g += nextColor.g * scale;
		b += nextColor.b * scale;
		a += scale;
	}
};


typedef Vector2<float> Vector2f;
typedef Vector3<float> Vector3f;
typedef Vector4<float> Vector4f;

typedef Vector2<double> Vector2d;
typedef Vector3<double> Vector3d;
typedef Vector4<double> Vector4d;

typedef Vector2<int> Vector2i;
typedef Vector3<int> Vector3i;
typedef Vector4<int> Vector4i;

typedef Vector2<unsigned int> Vector2ui;
typedef Vector3<unsigned int> Vector3ui;
typedef Vector4<unsigned int> Vector4ui;

typedef Vector2<unsigned short int> Vector2usi;
typedef Vector3<unsigned short int> Vector3usi;
typedef Vector4<unsigned short int> Vector4usi;

typedef Vector2<char> Vector2c;
typedef Vector3<char> Vector3c;
typedef Vector4<char> Vector4c;

typedef Vector2<unsigned char> Vector2uc;
typedef Vector3<unsigned char> Vector3uc;
typedef Vector4<unsigned char> Vector4uc;

typedef Vector2<size_t> Vector2st;
typedef Vector3<size_t> Vector3st;
typedef Vector4<size_t> Vector4st;

#endif
