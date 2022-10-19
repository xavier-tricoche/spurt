#ifndef QUATERNION_H
#define QUATERNION_H

#include <math.h>

#ifndef M_PI
	#define M_PI 3.14159265358979
#endif

class Quaternion {
public:
	float quaternion[4];

	Quaternion() {
		reset();
	}

	Quaternion(Quaternion *q) {
		quaternion[0] = q->quaternion[0];
		quaternion[1] = q->quaternion[1];
		quaternion[2] = q->quaternion[2];
		quaternion[3] = q->quaternion[3];
	}

	void reset() {
		quaternion[0] = 0;
		quaternion[1] = 0;
		quaternion[2] = 0;
		quaternion[3] = 1;
	}

	void compute(float angle_degrees, float axis[3]) {
		compute(angle_degrees, axis[0], axis[1], axis[2]);
	}

	void compute(float angle_degrees, float axis_x, float axis_y, float axis_z) {
		// normalize axis
		double div = sqrt(axis_x*axis_x + axis_y*axis_y + axis_z*axis_z);
		if (div == 0) div = 10e-10;
		float magInv = (float)(1.0 / div);

		quaternion[0] = axis_x*magInv;
		quaternion[1] = axis_y*magInv;
		quaternion[2] = axis_z*magInv;

		double angle_radians = angle_degrees*M_PI/180.0;

		// conversion of a rotation angle and axis to a quaternion
		float sinScale = (float)sin(angle_radians*0.5);
		quaternion[0] *= sinScale;  quaternion[1] *= sinScale;  quaternion[2] *= sinScale;

		quaternion[3] = (float)cos(angle_radians*0.5);
	}

	inline Quaternion operator+ (const Quaternion& q) const {
		Quaternion result;
		result.quaternion[0] = quaternion[0]*q.quaternion[3] + q.quaternion[0]*quaternion[3] + q.quaternion[1]*quaternion[2] - q.quaternion[2]*quaternion[1];
		result.quaternion[1] = quaternion[1]*q.quaternion[3] + q.quaternion[1]*quaternion[3] + q.quaternion[2]*quaternion[0] - q.quaternion[0]*quaternion[2];
		result.quaternion[2] = quaternion[2]*q.quaternion[3] + q.quaternion[2]*quaternion[3] + q.quaternion[0]*quaternion[1] - q.quaternion[1]*quaternion[0];
		result.quaternion[3] = quaternion[3]*q.quaternion[3] - (quaternion[0]*q.quaternion[0] + quaternion[1]*q.quaternion[1] + quaternion[2]*q.quaternion[2]);

		float magInv = (float)(1.0 / sqrt(result.quaternion[0]*result.quaternion[0] + result.quaternion[1]*result.quaternion[1] + 
									 result.quaternion[2]*result.quaternion[2] + result.quaternion[3]*result.quaternion[3]));

		result.quaternion[0] *= magInv;  result.quaternion[1] *= magInv;  result.quaternion[2] *= magInv;  result.quaternion[3] *= magInv;

		return result;
	}

	inline Quaternion& operator+=(const Quaternion& q) {
		*this = *this + q;
		return *this;
	}

	inline Quaternion& operator= (const Quaternion& q) {
		quaternion[0] = q.quaternion[0];  quaternion[1] = q.quaternion[1];  quaternion[2] = q.quaternion[2];  quaternion[3] = q.quaternion[3];
		return *this;
	}

	float* getRotationMatrix() {
		float *rotationMatrix = new float[16];
		rotationMatrix[0] = (float)1.0 - (float)2.0 * (quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2]);
		rotationMatrix[1] = (float)2.0 * (quaternion[0] * quaternion[1] - quaternion[2] * quaternion[3]);
		rotationMatrix[2] = (float)2.0 * (quaternion[2] * quaternion[0] + quaternion[1] * quaternion[3]);
		rotationMatrix[3] = 0.0;

		rotationMatrix[4] = (float)2.0 * (quaternion[0] * quaternion[1] + quaternion[2] * quaternion[3]);
		rotationMatrix[5] = (float)1.0 - (float)2.0 * (quaternion[2] * quaternion[2] + quaternion[0] * quaternion[0]);
		rotationMatrix[6] = (float)2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3]);
		rotationMatrix[7] = 0.0;

		rotationMatrix[8] = (float)2.0 * (quaternion[2] * quaternion[0] - quaternion[1] * quaternion[3]);
		rotationMatrix[9] = (float)2.0 * (quaternion[1] * quaternion[2] + quaternion[0] * quaternion[3]);
		rotationMatrix[10] = (float)1.0 - (float)2.0 * (quaternion[1] * quaternion[1] + quaternion[0] * quaternion[0]);
		rotationMatrix[11] = 0.0;

		rotationMatrix[12] = 0.0;
		rotationMatrix[13] = 0.0;
		rotationMatrix[14] = 0.0;
		rotationMatrix[15] = 1.0;

		return rotationMatrix;
	}

	void getRotationAxisAndAngle(float &angle_degrees, float axis[3]) {
		getRotationAxisAndAngle(angle_degrees, axis[0], axis[1], axis[2]);
	}

	void getRotationAxisAndAngle(float &angle_degrees, float &axis_x, float &axis_y, float &axis_z) {
		float mag = sqrt(quaternion[0]*quaternion[0] + quaternion[1]*quaternion[1] + quaternion[2]*quaternion[2]);
	
		// if the magnitude is zero, just make the angle zero
		if (mag == 0) angle_degrees = 0;

		else {
			mag = (float)1.0 / mag;
			angle_degrees = 2 * acos(quaternion[3]);
			angle_degrees = (float)(180.0 * angle_degrees / M_PI);
		}

		axis_x = quaternion[0] * mag;
		axis_y = quaternion[1] * mag;
		axis_z = quaternion[2] * mag;
	}
};

#endif