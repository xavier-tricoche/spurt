#ifndef MATRIX_2_X_2_H
#define MATRIX_2_X_2_H

#include "VectorN.h"
#include "MathUtils.h"
#include "Eigenvector.h"

// column ordered, i.e.
// | 0  2 |
// | 1  3 |

template <class T>
class Matrix2x2 {
public:
	Matrix2x2() { sizeOfArray = 4*sizeof(T);  setAsIdentity(); }
	Matrix2x2(T val) { sizeOfArray = 4*sizeof(T);  *this = val; }
	Matrix2x2(T mat[16]) { sizeOfArray = 4*sizeof(T);  memcpy(data, mat, sizeOfArray); }
	Matrix2x2(const Matrix2x2 &m) { sizeOfArray = 4*sizeof(T);  *this = m; }

	T* pointer() { return data; }
	inline T get(int x, int y) const { return data[y*2+x]; }
	inline void set(int x, int y, T val) { data[y*2+x] = val; }

	inline void transpose() { T temp = data[1];  data[1] = data[2]; }
	inline Matrix2x2<T> getTranspose() const { Matrix2x2<T> result(*this);  result.transpose();  return result; }

	inline Matrix2x2<T> operator- () const { return *this * -1; }

	inline Matrix2x2<T> operator+ (const Matrix2x2<T>& mat) const {
		Matrix2x2 m(*this);
		m.data[0] += mat.data[0];  m.data[1] += mat.data[1];  m.data[2] += mat.data[2];  m.data[3] += mat.data[3];
		return m;
	}
	inline Matrix2x2<T> operator- (const Matrix2x2<T>& mat) const {
		Matrix2x2 m(*this);
		m.data[0] -= mat.data[0];  m.data[1] -= mat.data[1];  m.data[2] -= mat.data[2];  m.data[3] -= mat.data[3];
		return m;
	}
	inline Matrix2x2<T> operator* (const Matrix2x2<T>& mat) const {
		Matrix2x2<T> result;
		result.data[0]  = data[0]*mat.data[0] + data[2]*mat.data[1];
		result.data[1]  = data[1]*mat.data[0] + data[3]*mat.data[1];
		result.data[2]  = data[0]*mat.data[2] + data[2]*mat.data[3];
		result.data[3]  = data[1]*mat.data[2] + data[3]*mat.data[3];
		return result;
	}
	inline Matrix2x2<T>& operator=  (const Matrix2x2<T>& mat) { memcpy(data, mat.data, sizeOfArray);  return *this; }
	inline Matrix2x2<T>& operator+= (const Matrix2x2<T>& mat) { *this = *this + mat;  return *this; }
	inline Matrix2x2<T>& operator-= (const Matrix2x2<T>& mat) { *this = *this - mat;  return *this; }
	inline Matrix2x2<T>& operator*= (const Matrix2x2<T>& mat) { *this = *this * mat;  return *this; }
	
	inline Vector2<T> operator* (const Vector2<T>& v) const {
		return Vector2<T>(v.x*get(0,0) + v.y*get(0,1), v.x*get(1,0) + v.y*get(1,1));
	}

	inline Matrix2x2<T>  operator+ (const T scalar) const { return *this + Matrix2x2(scalar); }
	inline Matrix2x2<T>  operator- (const T scalar) const { return *this - Matrix2x2(scalar); }
	inline Matrix2x2<T>  operator* (const T scalar) const {
		Matrix2x2<T> m(*this);
		m.data[0] *= scalar;  m.data[1] *= scalar;  m.data[2] *= scalar;  m.data[3] *= scalar;
		return m;
	}
	inline Matrix2x2<T>  operator/ (const T scalar) const { return *this * (T)(1.0/scalar); }

	inline Matrix2x2<T>& operator= (const T &scalar) {
		data[0] = data[1] = data[2] = data[3] = scalar;
		return *this;
	}
	inline Matrix2x2<T>& operator+= (const T &scalar) { return *this += Matrix2x2(scalar); }
	inline Matrix2x2<T>& operator-= (const T &scalar) { return *this -= Matrix2x2(scalar); }
	inline Matrix2x2<T>& operator*= (const T &scalar) { *this = *this * scalar;  return *this; }
	inline Matrix2x2<T>& operator/= (const T &scalar) { *this = *this / scalar;  return *this; }

	inline bool operator==(const Matrix2x2<T> &mat) const {
		return data[0] == mat.data[0] && data[1] == mat.data[1] && data[2] == mat.data[2] && data[3] == mat.data[3];
	}
	inline bool operator!=(const Matrix2x2<T> &mat) const { return !(*this == mat); }

	inline void setAsZero() { memset(data, 0, sizeOfArray); }
	inline void setAsIdentity() {
		memset(data, 0, sizeOfArray);
		data[0] = data[3] = 1;
	}
	bool isIdentity() {
		if (data[1] || data[2]) return false;
		if (data[0] != 1 || data[3] != 1) return false;
		return true;
	}

	inline void setDiagonal(const Vector2<T> &diag) { data[0] = diag.x;  data[3] = diag.y; }
	inline Vector2<T> getDiagonal() const { return Vector2<T>(data[0], data[3]); }

	void print() {
		cout << "| " << data[0] << " " << data[2] << " " << " |\n";
		cout << "| " << data[1] << " " << data[3] << " " << " |\n";
	}

	inline void getColumnVectors(Vector2<T> *v1, Vector2<T> *v2) const
		{ v1->set(data[0], data[1]);  v2->set(data[2], data[3]); }

	inline void setColumnsFromVectors(const Vector2<T> &v1, const Vector2<T> &v2) {
		data[0] = v1.x;  data[1] = v1.y;  data[2] = v2.x;  data[3] = v2.y;
	}
	inline void setColumnFromVector(const Vector2<T> &v, unsigned int column) {
		switch (column) {
			case 0: { data[0] = v.x;  data[1] = v.y;  break; }
			case 1: { data[2] = v.x;  data[3] = v.y;  break; }
		}
	}

	inline T getFrobeniusNorm() const {
		return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3]);
	}

	inline Matrix2x2<T> getFrobeniusNormalization() const {
		double scale = getFrobeniusNorm();
		if (scale != 0) return *this * 1.0/scale;
		else return *this;
	}

	inline T getTrace() const { return data[0]*data[3]; }
	inline T getDeterminant() const { return data[0]*data[3] - data[1]*data[2]; }

	inline void getEigenvalues(Vector2d &eigenvalue1, Vector2d &eigenvalue2) const {
		double trace = getTrace();
		double determinant = getDeterminant();
		// l^2 - tr.l + det = 0
		// delta = tr^2 - 4*det
		// lambda_max = 0.5*(tr + sqrt(delta)) = 0.5*(tr + sqrt(tr^2 - 4*det))

		double delta = trace*trace - 4.0f*determinant;
		eigenvalue1.x = eigenvalue2.x = 0.5*trace;
		if (delta>=0) {
			eigenvalue1.y = eigenvalue2.y = 0;
			double root = 0.5*sqrt(delta);
			eigenvalue1.x += root;
			eigenvalue2.x -= root;
		}
		else {
			double root = 0.5*sqrt(-delta);
			eigenvalue1.y = root;
			eigenvalue2.y = -root;
		}
	}

	inline void getEigenvectors(Eigenvector &eigenvector1, Eigenvector &eigenvector2) const {
		Vector2d eigenvalue1, eigenvalue2;
		getEigenvalues(eigenvalue1, eigenvalue2);

		if (eigenvalue1.y != 0) {
			eigenvector1.isReal = eigenvector2.isReal = false;
			return;
		}

		if (data[1] != 0) {
			eigenvector1.set(eigenvalue1.x, Vector2d(eigenvalue1.x-data[3], data[1]).data, 2);
			eigenvector2.set(eigenvalue2.x, Vector2d(eigenvalue2.x-data[3], data[1]).data, 2);
		}
		else if (data[2] != 0) {
			eigenvector1.set(eigenvalue1.x, Vector2d(data[2], eigenvalue1.x-data[0]).data, 2);
			eigenvector2.set(eigenvalue2.x, Vector2d(data[2], eigenvalue2.x-data[0]).data, 2);
		}
		else {
			eigenvector1.set(eigenvalue1.x, Vector2d(1,0).data, 2);
			eigenvector2.set(eigenvalue2.x, Vector2d(0,1).data, 2);
		}
	}

private:
	T data[4];
	int sizeOfArray;
};

typedef Matrix2x2<float> Matrix2x2f;
typedef Matrix2x2<double> Matrix2x2d;

#endif