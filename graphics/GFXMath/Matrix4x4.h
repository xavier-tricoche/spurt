#ifndef MATRIX_4_X_4_H
#define MATRIX_4_X_4_H

#include <math.h>
#include "VectorN.h"

// column ordered, i.e.
// | 0  4  8  12 |
// | 1  5  9  13 |
// | 2  6  10 14 |
// | 3  7  11 15 |

#ifndef M_PI
#define M_PI 3.14159265358979323846 
#endif

#ifndef DEGREE_TO_RADIAN
#define DEGREE_TO_RADIAN M_PI/(double)180.0;
#endif

template <class T>
class Matrix4x4 {
public:
	Matrix4x4() { sizeOfArray = 16*sizeof(T);  setAsIdentity(); }
	Matrix4x4(T val) { sizeOfArray = 16*sizeof(T);  *this = val; }
	Matrix4x4(T mat[16]) { sizeOfArray = 16*sizeof(T);  memcpy(data, mat, sizeOfArray); }
	Matrix4x4(const Matrix4x4 &m) { sizeOfArray = 16*sizeof(T);  *this = m; }

	inline T* pointer() { return data; }
	inline T get(int x, int y) const { return data[y*4+x]; }
	inline void set(int x, int y, T val) { data[y*4+x] = val; }

	inline void transpose() {
		T temp;
		temp = data[1];		data[1]  = data[4];		data[4]  = temp;
		temp = data[2];		data[2]  = data[8];		data[8]  = temp;
		temp = data[3];		data[3]  = data[12];	data[12] = temp;
		temp = data[6];		data[6]  = data[9];		data[9]  = temp;
		temp = data[7];		data[7]  = data[13];	data[13] = temp;
		temp = data[11];	data[11] = data[14];	data[14] = temp;
	}
	inline Matrix4x4<T> getTranspose() { Matrix4x4<T> result(*this);  result.transpose();  return result; }

	inline Matrix4x4<T> operator- () const { return *this * -1; }

	inline Matrix4x4<T> operator+ (const Matrix4x4<T>& mat) const {
		Matrix4x4 m(*this);
		m.data[0]  += mat.data[0];	m.data[1]  += mat.data[1];	m.data[2]  += mat.data[2];	m.data[3]  += mat.data[3];
		m.data[4]  += mat.data[4];	m.data[5]  += mat.data[5];	m.data[6]  += mat.data[6];	m.data[7]  += mat.data[7];
		m.data[8]  += mat.data[8];	m.data[9]  += mat.data[9];	m.data[10] += mat.data[10];	m.data[11] += mat.data[11];
		m.data[12] += mat.data[12];	m.data[13] += mat.data[13];	m.data[14] += mat.data[14];	m.data[15] += mat.data[15];
		return m;
	}
	inline Matrix4x4<T> operator- (const Matrix4x4<T>& mat) const {
		Matrix4x4 m(*this);
		m.data[0]  -= mat.data[0];	m.data[1]  -= mat.data[1];	m.data[2]  -= mat.data[2];	m.data[3]  -= mat.data[3];
		m.data[4]  -= mat.data[4];	m.data[5]  -= mat.data[5];	m.data[6]  -= mat.data[6];	m.data[7]  -= mat.data[7];
		m.data[8]  -= mat.data[8];	m.data[9]  -= mat.data[9];	m.data[10] -= mat.data[10];	m.data[11] -= mat.data[11];
		m.data[12] -= mat.data[12];	m.data[13] -= mat.data[13];	m.data[14] -= mat.data[14];	m.data[15] -= mat.data[15];
		return m;
	}
	inline Matrix4x4<T> operator* (const Matrix4x4<T>& mat) const {
		Matrix4x4<T> result;
		result.data[0]  = data[0]*mat.data[0] +data[4]*mat.data[1] +data[8] *mat.data[2] +data[12]*mat.data[3];
		result.data[1]  = data[1]*mat.data[0] +data[5]*mat.data[1] +data[9] *mat.data[2] +data[13]*mat.data[3];
		result.data[2]  = data[2]*mat.data[0] +data[6]*mat.data[1] +data[10]*mat.data[2] +data[14]*mat.data[3];
		result.data[3]  = data[3]*mat.data[0] +data[7]*mat.data[1] +data[11]*mat.data[2] +data[15]*mat.data[3];
		result.data[4]  = data[0]*mat.data[4] +data[4]*mat.data[5] +data[8] *mat.data[6] +data[12]*mat.data[7];
		result.data[5]  = data[1]*mat.data[4] +data[5]*mat.data[5] +data[9] *mat.data[6] +data[13]*mat.data[7];
		result.data[6]  = data[2]*mat.data[4] +data[6]*mat.data[5] +data[10]*mat.data[6] +data[14]*mat.data[7];
		result.data[7]  = data[3]*mat.data[4] +data[7]*mat.data[5] +data[11]*mat.data[6] +data[15]*mat.data[7];
		result.data[8]  = data[0]*mat.data[8] +data[4]*mat.data[9] +data[8] *mat.data[10]+data[12]*mat.data[11];
		result.data[9]  = data[1]*mat.data[8] +data[5]*mat.data[9] +data[9] *mat.data[10]+data[13]*mat.data[11];
		result.data[10] = data[2]*mat.data[8] +data[6]*mat.data[9] +data[10]*mat.data[10]+data[14]*mat.data[11];
		result.data[11] = data[3]*mat.data[8] +data[7]*mat.data[9] +data[11]*mat.data[10]+data[15]*mat.data[11];
		result.data[12] = data[0]*mat.data[12]+data[4]*mat.data[13]+data[8] *mat.data[14]+data[12]*mat.data[15];
		result.data[13] = data[1]*mat.data[12]+data[5]*mat.data[13]+data[9] *mat.data[14]+data[13]*mat.data[15];
		result.data[14] = data[2]*mat.data[12]+data[6]*mat.data[13]+data[10]*mat.data[14]+data[14]*mat.data[15];
		result.data[15] = data[3]*mat.data[12]+data[7]*mat.data[13]+data[11]*mat.data[14]+data[15]*mat.data[15];
		return result;
	}
	inline Matrix4x4<T>& operator=  (const Matrix4x4<T>& mat) { memcpy(data, mat.data, sizeOfArray);  return *this; }
	inline Matrix4x4<T>& operator+= (const Matrix4x4<T>& mat) { *this = *this + mat;  return *this; }
	inline Matrix4x4<T>& operator-= (const Matrix4x4<T>& mat) { *this = *this - mat;  return *this; }
	inline Matrix4x4<T>& operator*= (const Matrix4x4<T>& mat) { *this = *this * mat;  return *this; }
	
	inline Vector3<T> operator* (const Vector3<T> &v) const {
		Vector3<T> result;
		for (int i=0; i<3; i++) result.data[i] = v.x*get(i,0) + v.y*get(i,1) + v.z*get(i,2) + get(i,3);
		return result;
	}

	inline Vector4<T> operator* (const Vector4<T> &v) const {
		Vector4<T> result;
		for (int i=0; i<4; i++) result.data[i] = v.x*get(i,0) + v.y*get(i,1) + v.z*get(i,2) + v.w*get(i,3);
		return result;
	}

	inline Matrix4x4<T>  operator+ (const T scalar) const { return *this + Matrix4x4(scalar); }
	inline Matrix4x4<T>  operator- (const T scalar) const { return *this - Matrix4x4(scalar); }
	inline Matrix4x4<T>  operator* (const T scalar) const {
		Matrix4x4<T> m(*this);
		m.data[0]  *= scalar;	m.data[1]  *= scalar;	m.data[2]  *= scalar;	m.data[3]  *= scalar;
		m.data[4]  *= scalar;	m.data[5]  *= scalar;	m.data[6]  *= scalar;	m.data[7]  *= scalar;
		m.data[8]  *= scalar;	m.data[9]  *= scalar;	m.data[10] *= scalar;	m.data[11] *= scalar;
		m.data[12] *= scalar;	m.data[13] *= scalar;	m.data[14] *= scalar;	m.data[15] *= scalar;
		return m;
	}
	inline Matrix4x4<T>  operator/ (const T scalar) const { return *this * (T)(1.0/scalar); }

	inline Matrix4x4<T>& operator= (const T &scalar) {
		data[0] = data[1] = data[2]  = data[3]  = data[4]  = data[5]  = data[6]  = data[7]  =
		data[8] = data[9] = data[10] = data[11] = data[12] = data[13] = data[14] = data[15] = scalar;
		return *this;
	}
	inline Matrix4x4<T>& operator+= (const T &scalar) { return *this += Matrix4x4(scalar); }
	inline Matrix4x4<T>& operator-= (const T &scalar) { return *this -= Matrix4x4(scalar); }
	inline Matrix4x4<T>& operator*= (const T &scalar) { *this = *this * scalar;  return *this; }
	inline Matrix4x4<T>& operator/= (const T &scalar) { *this = *this / scalar;  return *this; }

	inline bool operator==(const Matrix4x4<T> &mat) const {
		if (data[0]  == mat.data[0]  && data[1]  == mat.data[1]  && data[2]  == mat.data[2]  && data[3]  == mat.data[3]  &&
			data[4]  == mat.data[4]  && data[5]  == mat.data[5]  && data[6]  == mat.data[6]  && data[7]  == mat.data[7]  &&
			data[8]  == mat.data[8]  && data[9]  == mat.data[9]  && data[10] == mat.data[10] && data[11] == mat.data[11] &&
			data[12] == mat.data[12] && data[13] == mat.data[13] && data[14] == mat.data[14] && data[15] == mat.data[15])
				return true;

		else return false;
	}
	inline bool operator!=(const Matrix4x4<T> &mat) const { return !(*this == mat); }
	
	inline void setAsZero() { memset(data, 0, sizeOfArray); }
	inline void setAsIdentity() {
		memset(data, 0, sizeOfArray);
		data[0] = data[5] = data[10] = data[15] = 1;
	}
	inline bool isIdentity() const {
		if (data[1] || data[2] || data[3] || data[4] || data[6] || data[7] ||
			data[8] || data[9] || data[11] || data[12] || data[13] || data[14]) return false;
		if (data[0] != 1 || data[5] != 1 || data[10] != 1 || data[15] != 1) return false;
		return true;
	}
	inline void setAsTranslate(const Vector3<T> &t) { setAsTranslate(t.x, t.y, t.z); }
	inline void setAsTranslate(T tx, T ty, T tz) {
		setAsIdentity();
		data[12] = tx;  data[13] = ty;  data[14] = tz;
	}
	inline Vector3<T> getTranslate() const { return Vector3<T>(data[12], data[13], data[14]); }

	inline void setAsScale(const Vector3<T> &s) { setAsScale(s.x, s.y, s.z); }
	inline void setAsScale(const T &sx, const T &sy, const T &sz) {
		setAsZero();
		data[0] = sx;  data[5] = sy;  data[10] = sz;  data[15] = 1;
	}
	inline Vector3<T> getScale() const { return Vector3<T>(data[0], data[5], data[10]); }

	inline void setAsRotateX(const double degrees) {
		double rad = degrees*DEGREE_TO_RADIAN;
		setAsIdentity();
		data[5] = (T)cos(rad);	data[9] = (T)sin(rad);
		data[6] = -data[9];	data[10] = data[5];
	}
	inline void setAsRotateY(const double degrees) {
		double rad = degrees*DEGREE_TO_RADIAN;
		setAsIdentity();
		data[0] = (T)cos(rad);	data[2] = (T)sin(rad);
		data[8] = -data[2];	data[10] = data[0];
	}
	inline void setAsRotateZ(const double degrees) {
		double rad = degrees*DEGREE_TO_RADIAN;
		setAsIdentity();
		data[0] = (T)cos(rad);	data[4] = (T)sin(rad);
		data[1] = -data[4];	data[5] = data[0];
	}
	inline void setAsRotate(const Vector3<T> &axis, const double degrees) { setAsRotate(axis.x, axis.y, axis.z, degrees); }
	inline void setAsRotate(const double axisX, const double axisY, const double axisZ, const double degrees) {
		if (degrees == 0) { setAsIdentity();  return; }

		// compute quaternion from given axis and angle
		double magInv = 1.0 / sqrt(axisX*axisX + axisY*axisY + axisZ*axisZ);

		double quaternion[4];
		quaternion[0] = axisX*magInv;
		quaternion[1] = axisY*magInv;
		quaternion[2] = axisZ*magInv;

		double angle = degrees*DEGREE_TO_RADIAN;

		// conversion of a rotation angle and axis to a quaternion
		T sinScale = (T)sin(angle*0.5);
		quaternion[0] *= sinScale;  quaternion[1] *= sinScale;  quaternion[2] *= sinScale;

		quaternion[3] = (T)cos(angle*0.5);

		// compute the actual rotation matrix
		data[0] = (T)(1.0 - 2.0 * (quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2]));
		data[4] = (T)(2.0 * (quaternion[0] * quaternion[1] - quaternion[2] * quaternion[3]));
		data[8] = (T)(2.0 * (quaternion[2] * quaternion[0] + quaternion[1] * quaternion[3]));

		data[1] = (T)(2.0 * (quaternion[0] * quaternion[1] + quaternion[2] * quaternion[3]));
		data[5] = (T)(1.0 - 2.0 * (quaternion[2] * quaternion[2] + quaternion[0] * quaternion[0]));
		data[9] = (T)(2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3]));

		data[2] =  (T)(2.0 * (quaternion[2] * quaternion[0] - quaternion[1] * quaternion[3]));
		data[6] =  (T)(2.0 * (quaternion[1] * quaternion[2] + quaternion[0] * quaternion[3]));
		data[10] = (T)(1.0 - 2.0 * (quaternion[1] * quaternion[1] + quaternion[0] * quaternion[0]));

		data[3] = data[7] = data[11] = data[12] = data[13] = data[14] = 0;
		data[15] = 1;
	}

	inline void print() const {
		cout << "| " << data[0] << " " << data[4] << " " << data[8]  << " " << data[12] << " |\n";
		cout << "| " << data[1] << " " << data[5] << " " << data[9]  << " " << data[13] << " |\n";
		cout << "| " << data[2] << " " << data[6] << " " << data[10] << " " << data[14] << " |\n";
		cout << "| " << data[3] << " " << data[7] << " " << data[11] << " " << data[15] << " |\n";
	}

	inline void invert() { *this = getInverse(); }
	inline Matrix4x4<T> getInverse() const {
		Matrix4x4<T> result = getInverseTranspose();
		result.transpose();
		return result;
	}

	inline void invertTranspose() { *this = getInverseTranspose(); }
	inline Matrix4x4<T> getInverseTranspose() const {
		Matrix4x4<T> result;

		// calculate pairs for first 8 elements (cofactors)
		float temp[12];	
		temp[0]  = data[10]*data[15];
		temp[1]  = data[11]*data[14];
		temp[2]  = data[9] *data[15];
		temp[3]  = data[11]*data[13];
		temp[4]  = data[9] *data[14];
		temp[5]  = data[10]*data[13];
		temp[6]  = data[8] *data[15];
		temp[7]  = data[11]*data[12];
		temp[8]  = data[8] *data[14];
		temp[9]  = data[10]*data[12];
		temp[10] = data[8]*data[13];
		temp[11] = data[9]*data[12];

		//calculate first 8 cofactors
		result.data[0] = temp[0]*data[5]+temp[3]*data[6]+temp[4] *data[7]-temp[1]*data[5]-temp[2]*data[6]-temp[5] *data[7];
		result.data[1] = temp[1]*data[4]+temp[6]*data[6]+temp[9] *data[7]-temp[0]*data[4]-temp[7]*data[6]-temp[8] *data[7];
		result.data[2] = temp[2]*data[4]+temp[7]*data[5]+temp[10]*data[7]-temp[3]*data[4]-temp[6]*data[5]-temp[11]*data[7];
		result.data[3] = temp[5]*data[4]+temp[8]*data[5]+temp[11]*data[6]-temp[4]*data[4]-temp[9]*data[5]-temp[10]*data[6];
		result.data[4] = temp[1]*data[1]+temp[2]*data[2]+temp[5] *data[3]-temp[0]*data[1]-temp[3]*data[2]-temp[4] *data[3];
		result.data[5] = temp[0]*data[0]+temp[7]*data[2]+temp[8] *data[3]-temp[1]*data[0]-temp[6]*data[2]-temp[9] *data[3];
		result.data[6] = temp[3]*data[0]+temp[6]*data[1]+temp[11]*data[3]-temp[2]*data[0]-temp[7]*data[1]-temp[10]*data[3];
		result.data[7] = temp[4]*data[0]+temp[9]*data[1]+temp[10]*data[2]-temp[5]*data[0]-temp[8]*data[1]-temp[11]*data[2];

		//calculate pairs for second 8 cofactors
		temp[0]  = data[2]*data[7];
		temp[1]  = data[3]*data[6];
		temp[2]  = data[1]*data[7];
		temp[3]  = data[3]*data[5];
		temp[4]  = data[1]*data[6];
		temp[5]  = data[2]*data[5];
		temp[6]  = data[0]*data[7];
		temp[7]  = data[3]*data[4];
		temp[8]  = data[0]*data[6];
		temp[9]  = data[2]*data[4];
		temp[10] = data[0]*data[5];
		temp[11] = data[1]*data[4];

		// calculate second 8 cofactors
		result.data[8]  = temp[0] *data[13]+temp[3]*data[14] +temp[4] *data[15]-temp[1] *data[13]-temp[2] *data[14]-temp[5] *data[15];
		result.data[9]  = temp[1] *data[12]+temp[6]*data[14] +temp[9] *data[15]-temp[0] *data[12]-temp[7] *data[14]-temp[8] *data[15];
		result.data[10] = temp[2] *data[12]+temp[7]*data[13] +temp[10]*data[15]-temp[3] *data[12]-temp[6] *data[13]-temp[11]*data[15];
		result.data[11] = temp[5] *data[12]+temp[8]*data[13] +temp[11]*data[14]-temp[4] *data[12]-temp[9] *data[13]-temp[10]*data[14];
		result.data[12] = temp[2] *data[10]+temp[5]*data[11] +temp[1] *data[9] -temp[4] *data[11]-temp[0] *data[9] -temp[3] *data[10];
		result.data[13] = temp[8] *data[11]+temp[0]*data[8]  +temp[7] *data[10]-temp[6] *data[10]-temp[9] *data[11]-temp[1] *data[8];
		result.data[14] = temp[6] *data[9] +temp[11]*data[11]+temp[3] *data[8] -temp[10]*data[11]-temp[2] *data[8] -temp[7] *data[9];
		result.data[15] = temp[10]*data[10]+temp[4]*data[8]  +temp[9] *data[9] -temp[8] *data[9] -temp[11]*data[10]-temp[5] *data[8];

		float determinant = data[0]*result.data[0]+data[1]*result.data[1]+data[2]*result.data[2]+data[3]*result.data[3];

		if (determinant == 0) { return Matrix4x4<T>(); }
		else return result/determinant;
	}

	inline void setAsPerspective(double left, double right, double bottom, double top, double zNear, double zFar) {
		double dx = right-left;
		double dy = top-bottom;
		double dz = zFar-zNear;

		if (dx == 0 || dy == 0 || dz == 0) return;

		setAsZero();

		data[0] = (T)((2.0*zNear)/dx);
		data[5] = (T)((2.0*zNear)/dy);
		data[8] = (T)((right+left)/dx);
		data[9] = (T)((top+bottom)/dy);
		data[10] = -(T)((zFar+zNear)/dz);
		data[11] = -1;
		data[14] = -(T)((2.0*zFar*zNear)/dz);
	}
	inline void setAsPerspective(double fovy_degrees, double aspect, double zNear, double zFar) {
		double rad = fovy_degrees * (double)M_PI/180.0;
		double top = zNear*tanf(rad * 0.5);
		double bottom = -top;
		setAsPerspective(aspect*bottom, aspect*top, bottom, top, zNear, zFar);
	}
	inline void setAsOrthographic(double left, double right, double bottom, double top, double zNear, double zFar) {
		setAsIdentity();
		double dx = right-left;
		double dy = top-bottom;
		double dz = zFar-zNear;
		data[0]  = (T)(2.0/dx);
		data[5]  = (T)(2.0/dy);
		data[10] = -(T)(2.0/dz);
		data[12] = -(T)((right+left)/dx);
		data[13] = -(T)((top+bottom)/dy);
		data[14] = -(T)((zFar+zNear)/dz);
	}

	template <class S> inline void setAsCamera(S eyex, S eyey, S eyez, S centerx, S centery, S centerz, S upx, S upy, S upz)
		{ setAsCamera(Vector3<S>(eyex, eyey, eyez), Vector3<S>(centerx, centery, centerz), Vector3<S>(upx, upy, upz)); }
	template <class S> inline void setAsCamera(Vector3<S> eye, Vector3<S> lookAtPos, Vector3<S> upVector) {
		Vector3<S> dir = lookAtPos-eye;
		Vector3<S> left = dir.cross(upVector);
		Vector3<S> up = left.cross(dir);

		dir.normalize();
		left.normalize();
		up.normalize();

		data[0] = (T)left.x;	data[4] = (T)left.y;	data[8]  = (T)left.z;	data[12] = -(T)left.dot(eye);
		data[1] = (T)up.x;		data[5] = (T)up.y;		data[9]  = (T)up.z;		data[13] = -(T)up.dot(eye);
		data[2] = (T)-dir.x;	data[6] = (T)-dir.y;	data[10] = (T)-dir.z;	data[14] =  (T)dir.dot(eye);
		data[3] = 0;			data[7] = 0;			data[11] = 0;			data[15] = 1;
	}

	template <class S> inline void getCameraProperties(Vector3<S> &eye, Vector3<S> &lookAtPos, Vector3<S> &upVector) const {
		upVector.set((S)get(1,0), (S)get(1,1), (S)get(1,2));
		upVector.normalize();
		Vector3<S> lookDir((S)-get(2,0), (S)-get(2,1), (S)-get(2,2));
		lookDir.normalize();

		Matrix4x4<T> invertedMat = getInverse();

		eye.set((S)invertedMat.get(0,3), (S)invertedMat.get(1,3), (S)invertedMat.get(2,3));
		lookAtPos = eye+lookDir;
	}

#ifdef __GL_H__
	inline void multiplyOpenGLMatrix() const { glMultMatrixf(data); }
	inline void setCurrentOpenGLMatrix() const { glLoadMatrixf(data); }
	inline void setAsCurrentOpenGLModelviewMatrix() const { glGetFloatv(GL_MODELVIEW_MATRIX, data); }
	inline void setAsCurrentOpenGLProjectionMatrix() const { glGetFloatv(GL_PROJECTION_MATRIX, data); }
	inline void setAsCurrentOpenGLTextureMatrix() const { glGetFloatv(GL_TEXTURE_MATRIX, data); }
#endif

private:
	T data[16];
	unsigned int sizeOfArray;
};

typedef Matrix4x4<float> Matrix4x4f;

#endif