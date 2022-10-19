#ifndef MATRIX_3_X_3_H
#define MATRIX_3_X_3_H

#include "VectorN.h"
#include "MathUtils.h"
#include "Eigenvector.h"

// column ordered, i.e.
// | 0  3  6 |
// | 1  4  7 |
// | 2  5  8 |

template <class T>
class Matrix3x3 {
public:
	Matrix3x3() { sizeOfArray = 9*sizeof(T);  setAsIdentity(); }
	Matrix3x3(T val) { sizeOfArray = 9*sizeof(T);  *this = val; }
	Matrix3x3(T mat[16]) { sizeOfArray = 9*sizeof(T);  memcpy(data, mat, sizeOfArray); }
	Matrix3x3(const Matrix3x3 &m) { sizeOfArray = 9*sizeof(T);  *this = m; }

	T* pointer() { return data; }
	inline T get(int x, int y) const { return data[y*3+x]; }
	inline void set(int x, int y, T val) { data[y*3+x] = val; }

	inline void transpose() {
		T temp;
		temp = data[1];  data[1] = data[3];  data[3] = temp;
		temp = data[2];  data[2] = data[6];  data[6] = temp;
		temp = data[5];  data[5] = data[7];  data[7] = temp;
	}
	inline Matrix3x3<T> getTranspose() const { Matrix3x3<T> result(*this);  result.transpose();  return result; }

	inline Matrix3x3<T> operator- () const { return *this * -1; }

	inline Matrix3x3<T> operator+ (const Matrix3x3<T>& mat) const {
		Matrix3x3 m(*this);
		m.data[0] += mat.data[0];	m.data[1] += mat.data[1];	m.data[2] += mat.data[2];
		m.data[3] += mat.data[3];	m.data[4] += mat.data[4];	m.data[5] += mat.data[5];
		m.data[6] += mat.data[6];	m.data[7] += mat.data[7];	m.data[8] += mat.data[8];
		return m;
	}
	inline Matrix3x3<T> operator- (const Matrix3x3<T>& mat) const {
		Matrix3x3 m(*this);
		m.data[0] -= mat.data[0];	m.data[1] -= mat.data[1];	m.data[2] -= mat.data[2];
		m.data[3] -= mat.data[3];	m.data[4] -= mat.data[4];	m.data[5] -= mat.data[5];
		m.data[6] -= mat.data[6];	m.data[7] -= mat.data[7];	m.data[8] -= mat.data[8];
		return m;
	}
	inline Matrix3x3<T> operator* (const Matrix3x3<T>& mat) const {
		Matrix3x3<T> result;
		result.data[0]  = data[0]*mat.data[0] + data[3]*mat.data[1] + data[6]*mat.data[2];
		result.data[1]  = data[1]*mat.data[0] + data[4]*mat.data[1] + data[7]*mat.data[2];
		result.data[2]  = data[2]*mat.data[0] + data[5]*mat.data[1] + data[8]*mat.data[2];
		result.data[3]  = data[0]*mat.data[3] + data[3]*mat.data[4] + data[6]*mat.data[5];
		result.data[4]  = data[1]*mat.data[3] + data[4]*mat.data[4] + data[7]*mat.data[5];
		result.data[5]  = data[2]*mat.data[3] + data[5]*mat.data[4] + data[8]*mat.data[5];
		result.data[6]  = data[0]*mat.data[6] + data[3]*mat.data[7] + data[6]*mat.data[8];
		result.data[7]  = data[1]*mat.data[6] + data[4]*mat.data[7] + data[7]*mat.data[8];
		result.data[8]  = data[2]*mat.data[6] + data[5]*mat.data[7] + data[8]*mat.data[8];
		return result;
	}
	inline Matrix3x3<T>& operator=  (const Matrix3x3<T>& mat) { memcpy(data, mat.data, sizeOfArray);  return *this; }
	inline Matrix3x3<T>& operator+= (const Matrix3x3<T>& mat) { *this = *this + mat;  return *this; }
	inline Matrix3x3<T>& operator-= (const Matrix3x3<T>& mat) { *this = *this - mat;  return *this; }
	inline Matrix3x3<T>& operator*= (const Matrix3x3<T>& mat) { *this = *this * mat;  return *this; }
	
	inline Vector3<T> operator* (const Vector3<T>& v) const {
		Vector3<T> result;
		for (int i=0; i<3; i++) { result.data[i] = v.x*get(i,0) + v.y*get(i,1) + v.z*get(i,2); }
		return result;
	}

	inline Matrix3x3<T>  operator+ (const T scalar) const { return *this + Matrix3x3(scalar); }
	inline Matrix3x3<T>  operator- (const T scalar) const { return *this - Matrix3x3(scalar); }
	inline Matrix3x3<T>  operator* (const T scalar) const {
		Matrix3x3<T> m(*this);
		m.data[0]  *= scalar;	m.data[1]  *= scalar;	m.data[2]  *= scalar;
		m.data[3]  *= scalar;	m.data[4]  *= scalar;	m.data[5]  *= scalar;
		m.data[6]  *= scalar;	m.data[7]  *= scalar;	m.data[8]  *= scalar;
		return m;
	}
	inline Matrix3x3<T>  operator/ (const T scalar) const { return *this * (T)(1.0/scalar); }

	inline Matrix3x3<T>& operator= (const T &scalar) {
		data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = data[6] = data[7] = data[8] = scalar;
		return *this;
	}
	inline Matrix3x3<T>& operator+= (const T &scalar) { return *this += Matrix3x3(scalar); }
	inline Matrix3x3<T>& operator-= (const T &scalar) { return *this -= Matrix3x3(scalar); }
	inline Matrix3x3<T>& operator*= (const T &scalar) { *this = *this * scalar;  return *this; }
	inline Matrix3x3<T>& operator/= (const T &scalar) { *this = *this / scalar;  return *this; }

	inline bool operator==(const Matrix3x3<T> &mat) const {
		if (data[0] == mat.data[0] && data[1] == mat.data[1] && data[2] == mat.data[2] &&
			data[3] == mat.data[3] && data[4] == mat.data[4] && data[5] == mat.data[5] &&
			data[6] == mat.data[6] && data[7] == mat.data[7] && data[8] == mat.data[8])
				return true;

		else return false;
	}
	inline bool operator!=(const Matrix3x3<T> &mat) const { return !(*this == mat); }

	inline void setAsZero() { memset(data, 0, sizeOfArray); }
	inline bool isZero() { return !(data[0] || data[1] || data[2] || data[3] || data[4] || data[5] || data[6] || data[7] || data[8]); }

	inline void setAsIdentity() {
		memset(data, 0, sizeOfArray);
		data[0] = data[4] = data[8] = 1;
	}
	inline bool isIdentity() {
		if (data[1] || data[2] || data[3] || data[5] || data[6] || data[7]) return false;
		if (data[0] != 1 || data[4] != 1 || data[8] != 1) return false;
		return true;
	}

	inline void setDiagonal(const Vector3<T> &diag) { data[0] = diag.x;  data[4] = diag.y;  data[8] = diag.z; }
	inline Vector3<T> getDiagonal() const { return Vector3<T>(data[0], data[4], data[8]); }

	void print() {
		cout << "| " << data[0] << " " << data[3] << " " << data[6]  << " " << " |\n";
		cout << "| " << data[1] << " " << data[4] << " " << data[7]  << " " << " |\n";
		cout << "| " << data[2] << " " << data[5] << " " << data[8]  << " " << " |\n";
	}

	inline void getColumnVectors(Vector3<T> *v1, Vector3<T> *v2, Vector3<T> *v3) const
		{ v1->set(data[0], data[1], data[2]);  v2->set(data[3], data[4], data[5]);  v3->set(data[6], data[7], data[8]); }
	template <class S> inline void getColumnVectors(S *v1, S *v2, S *v3) const {
		v1[0] = data[0];  v1[1] = data[1];  v1[2] = data[2];
		v2[0] = data[3];  v2[1] = data[4];  v2[2] = data[5];
		v3[0] = data[6];  v3[1] = data[7];  v3[2] = data[8];
	}

	inline void setColumnsFromVectors(const Vector3<T> &v1, const Vector3<T> &v2, const Vector3<T> &v3) {
		data[0] = v1.x;  data[1] = v1.y;  data[2] = v1.z;
		data[3] = v2.x;  data[4] = v2.y;  data[5] = v2.z;
		data[6] = v3.x;  data[7] = v3.y;  data[8] = v3.z;
	}
	inline void setColumnFromVector(const Vector3<T> &v, unsigned int column) {
		switch (column) {
			case 0: { data[0] = v.x;  data[1] = v.y;  data[2] = v.z;  break; }
			case 1: { data[3] = v.x;  data[4] = v.y;  data[5] = v.z;  break; }
			case 2: { data[6] = v.x;  data[7] = v.y;  data[8] = v.z;  break; }
		}
	}

	// nullspace stuff and eigen-value/vector stuff taken from http://teem.sourceforge.net/
	inline void allignColumnVectors() {
		Vector3<T> v1, v2, v3;
		getColumnVectors(&v1, &v2, &v3);
		
		unsigned int maxIndex = Vector3<T>(v1.magnitudeSquared(), v2.magnitudeSquared(), v3.magnitudeSquared()).getIndexOfMaxValue();

		Vector3<T> *t1, *t2, *t3;
		if (maxIndex == 0)		{ t1 = &v1;  t2 = &v2;  t3 = &v3; }
		else if (maxIndex == 1) { t1 = &v2;  t2 = &v3;  t3 = &v1; }
		else					{ t1 = &v3;  t2 = &v1;  t3 = &v2; }

		if (t1->dot(*t2) < 0) *t2 *= -1;
		if (t1->dot(*t3) < 0) *t3 *= -1;

		setColumnsFromVectors(v1, v2, v3);
	}

	// makes sure column vectors are mutually orthogonal, first vector is unchanged
	inline void enforceOrthogonality() {
		Vector3<T> v0, v1, v2;
		getColumnVectors(&v0, &v1, &v2);

		double d00 = v0.dot(v0), d10 = v1.dot(v0), d11 = v1.dot(v1);
		Vector3<T> temp = v1 + -d10/d00 * v0;
		v1 = temp * sqrt(d11 / temp.dot(temp));
		
		double d20 = v2.dot(v0), d21 = v2.dot(v1), d22 = v2.dot(v2);
		temp = v2 + -d20/d00 * v0 + -d21/d00 * v1;
		v2 = temp * sqrt(d22 / temp.dot(temp));

		setColumnFromVector(v1, 1);
		setColumnFromVector(v2, 2);
	}

	// makes sure 3rd column vector has a positive dot product with cross product of first 2 column vectors
	inline void makeRightHanded() {
		Vector3<T> v0, v1, v2;
		getColumnVectors(&v0, &v1, &v2);
		if (v2.dot(v0.cross(v1)) < 0) setColumnFromVector(v2*-1, 2);
	}

	// get normalized vector that spans the nullspace (matrix assumed to have nullspace of dimension 1)
	inline Vector3<T> getNullSpace1D() const {
		Matrix3x3<T> n = getTranspose();
		
		Vector3<T> v0, v1, v2;
		n.getColumnVectors(&v0, &v1, &v2);

		Vector3<T> t1 = v0.cross(v1), t2 = v0.cross(v2), t3 = v1.cross(v2);
		n.setColumnsFromVectors(t1, t2, t3);
		n.allignColumnVectors();

		n.getColumnVectors(&v0, &v1, &v2);
		return (v0+v1+v2).unit();  // add them up - largest mag. (the more accurate) should dominate
	}

	// get 2 normalized vectors that spans the nullspace (matrix assumed to have nullspace of dimension 2)
	inline void getNullSpace2D(Vector3<T> *nullSpace1, Vector3<T> *nullSpace2) const {
		Matrix3x3<T> n = getTranspose();
		n.allignColumnVectors();

		Vector3<T> v0, v1, v2;
		n.getColumnVectors(&v0, &v1, &v2);

		Vector3<T> t = v0+v1+v2;
		t.normalize();
		
		// any two vectors which are perpendicular to the (supposedly 1D) span of the column vectors span the nullspace
		*nullSpace1 = t.getPerpendicular(); 
		nullSpace1->normalize();
		*nullSpace2 = t.cross(*nullSpace1);
	}
	template <class S> inline void getNullSpace2D(S *nullSpace1, S *nullSpace2) const {
		Vector3<T> n1, n2;			getNullSpace2D(&n1, &n2);
		n1.copyTo(nullSpace1);		n2.copyTo(nullSpace2);
	}

	inline T getFrobeniusNorm() const {
		return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3] + data[4]*data[4] +
					data[5]*data[5] + data[6]*data[6] + data[7]*data[7] + data[8]*data[8]);
	}

	inline Matrix3x3<T> getFrobeniusNormalization() const {
		double scale = getFrobeniusNorm();
		if (scale != 0) return *this * 1.0/scale;
		else return *this;
	}

	// finds eigenvalues of given matrix, returns cubic root information from getCubicRoots
	inline int getEigenvalues(double &eigenvalue0, double &eigenvalue1, double &eigenvalue2, const bool refineWithNewtonRaphson = true, const double eps = 1.0E-11) const {
		// frobenius normalization avoids the creation of imaginary eigenvalues when the coefficients of the matrix are large (> 50000)
		double scale = getFrobeniusNorm();
		Matrix3x3<T> m(scale != 0 ? *this * 1.0/scale : *this);
		
		double A = -m.data[0] - m.data[4] - m.data[8];
		double B = m.data[0]*m.data[4] - m.data[3]*m.data[1] + m.data[0]*m.data[8] - m.data[6]*m.data[2] + m.data[4]*m.data[8] - m.data[7]*m.data[5];
		double C = (m.data[6]*m.data[4] - m.data[3]*m.data[7])*m.data[2] + (m.data[0]*m.data[7] - m.data[6]*m.data[1])*m.data[5] + (m.data[3]*m.data[1] - m.data[0]*m.data[4])*m.data[8];

		int roots = MathUtils::getCubicRoots(eigenvalue0, eigenvalue1, eigenvalue2, A, B, C, refineWithNewtonRaphson, eps);
		
		double invScale = 1.0/scale;
		eigenvalue0 *= invScale;
		eigenvalue1 *= invScale;
		eigenvalue2 *= invScale;
		return roots;
	}

	// finds eigenvalues and eigenvectors of given matrix m, eigenvalues (and associated eigenvectors) are sorted in descending order
	inline void getEigenvectors(Eigenvector &eigenvector0, Eigenvector &eigenvector1, Eigenvector &eigenvector2,
								const bool refineWithNewtonRaphson = true, const double eps = 1.0E-11) const {
		Vector3<T> eigenvalues;
		int roots = getEigenvalues(eigenvector0.value, eigenvector1.value, eigenvector2.value, refineWithNewtonRaphson, eps);

		Vector3<T> eval(eigenvector0.value, eigenvector1.value, eigenvector2.value);

		Matrix3x3 n(*this);

		switch (roots) {
			case CUBIC_ROOT_THREE: {
				Vector3<T> diag = getDiagonal();
				n.setDiagonal(diag-eval.x);  n.getNullSpace1D().copyTo(eigenvector0.vec);
				n.setDiagonal(diag-eval.y);  n.getNullSpace1D().copyTo(eigenvector1.vec);
				n.setDiagonal(diag-eval.z);  n.getNullSpace1D().copyTo(eigenvector2.vec);

				n.setColumnsFromVectors(eigenvector0.vec, eigenvector1.vec, eigenvector2.vec);
				n.enforceOrthogonality();
				n.makeRightHanded();
				n.getColumnVectors(eigenvector0.vec, eigenvector1.vec, eigenvector2.vec);

				eigenvector0.isReal = true;  eigenvector1.isReal = true;  eigenvector2.isReal = true;
				break;
		   }

			case CUBIC_ROOT_SINGLE_DOUBLE: {
				Vector3<T> diag = getDiagonal();

				if (eval.x != eval.y) { // eval.x big, rest small - cigar shape
					n.setDiagonal(diag-eval.x);  n.getNullSpace1D().copyTo(eigenvector0.vec);
					n.setDiagonal(diag-eval.y);  n.getNullSpace2D(eigenvector1.vec, eigenvector2.vec);
				}
				else { // eval.z small, rest big - pancake shape
					n.setDiagonal(diag-eval.x);  n.getNullSpace2D(eigenvector0.vec, eigenvector1.vec);
					n.setDiagonal(diag-eval.z);  n.getNullSpace1D().copyTo(eigenvector2.vec);
				}

				n.setColumnsFromVectors(eigenvector0.vec, eigenvector1.vec, eigenvector2.vec);
				n.enforceOrthogonality();
				n.makeRightHanded();
				n.getColumnVectors(eigenvector0.vec, eigenvector1.vec, eigenvector2.vec);

				eigenvector0.isReal = true;  eigenvector1.isReal = true;  eigenvector2.isReal = true;
				break;
			}

			case CUBIC_ROOT_TRIPLE:
				eigenvector0.setEigenvector(1,0,0);  eigenvector1.setEigenvector(0,1,0);  eigenvector2.setEigenvector(0,0,1);
				eigenvector0.isReal = true;  eigenvector1.isReal = true;  eigenvector2.isReal = true;
				break;

			case CUBIC_ROOT_SINGLE:
				Vector3<T> diag = getDiagonal();
				n.setDiagonal(diag-eval.x);  n.getNullSpace1D().copyTo(eigenvector0.vec);
				eigenvector0.isReal = true;  eigenvector1.isReal = false;  eigenvector2.isReal = false;
				break;
		}
	}

private:
	T data[9];
	int sizeOfArray;
};

typedef Matrix3x3<float> Matrix3x3f;
typedef Matrix3x3<double> Matrix3x3d;

#endif