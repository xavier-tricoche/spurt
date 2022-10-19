/* Created by Nate Andrysco
   October 2005            */

#ifndef VIRTUAL_TRACKBALL_H
#define VIRTUAL_TRACKBALL_H

#include "graphics/GFXMath/Quaternion.h"

#include <math.h>

#if 0
#include <GL/gl.h>
#include "GL/glut.h"
#else
#include <glut.h>
#endif

class VirtualTrackball {
public:
	VirtualTrackball() { reset(); }

	inline VirtualTrackball& operator= (const VirtualTrackball &trackball) {
		reset();

		compositeAngle = trackball.compositeAngle;
		compositeAxis[0] = trackball.compositeAxis[0];
		compositeAxis[1] = trackball.compositeAxis[1];
		compositeAxis[2] = trackball.compositeAxis[2];
		compositeQuaternion = trackball.compositeQuaternion;
		rotationSensitivity[0] = trackball.rotationSensitivity[0];
		rotationSensitivity[1] = trackball.rotationSensitivity[1];

		scale = trackball.scale;

		translation[0] = trackball.translation[0];
		translation[1] = trackball.translation[1];
		translation[2] = trackball.translation[2];
		translationSensitivity[0] = trackball.translationSensitivity[0];
		translationSensitivity[1] = trackball.translationSensitivity[1];
		
		return *this;
	}

	inline void reset() { resetRotation();  resetScale();  resetTranslate(); }

	// Rotation Procedures
	inline void resetRotation() {
		rotationMouseDown = false;
		rotationStartPosition[0] = rotationStartPosition[1] = rotationStartPosition[2] = 0;
		movementQuaternion.reset();
		compositeQuaternion.reset();
		compositeAngle = compositeAxis[0] = compositeAxis[1] = compositeAxis[2] = 0;
		rotationSensitivity[0] = rotationSensitivity[1] = 1;
	}

	inline void trackballRecordMouseDown(const float x, const float y, const float screenWidth, const float screenHeight) {
		rotationMouseDown = true;

		movementQuaternion.reset();

		// normalize mouse position to [-1, 1]
		rotationStartPosition[0] = (2.0*x - screenWidth) / screenWidth * rotationSensitivity[0];
		rotationStartPosition[1] = (2.0*y - screenHeight) / screenHeight * rotationSensitivity[1];

		// compute arcball distances
		float d = sqrt(rotationStartPosition[0]*rotationStartPosition[0] + rotationStartPosition[1]*rotationStartPosition[1]);
		rotationStartPosition[2] = cos((M_PI/2.0) * ((d < 1.0) ? d : 1.0));
		float a = 1.0 / sqrt(rotationStartPosition[0]*rotationStartPosition[0] + rotationStartPosition[1]*rotationStartPosition[1]
							 + rotationStartPosition[2]*rotationStartPosition[2]);
		rotationStartPosition[0] *= a;
		rotationStartPosition[1] *= a;
		rotationStartPosition[2] *= a;
	}

	inline void trackballRecordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		if (!rotationMouseDown) return;

		float rotationMovementPosition[3];
		rotationMovementPosition[0] = (2.0*x - screenWidth) / screenWidth * rotationSensitivity[0];
		rotationMovementPosition[1] = (2.0*y - screenHeight) / screenHeight * rotationSensitivity[1];

		float d = sqrt(rotationMovementPosition[0]*rotationMovementPosition[0] + rotationMovementPosition[1]*rotationMovementPosition[1]);
		rotationMovementPosition[2] = cos((M_PI/2.0) * ((d < 1.0) ? d : 1.0));
		float a = 1.0 / sqrt(rotationMovementPosition[0]*rotationMovementPosition[0] + rotationMovementPosition[1]*rotationMovementPosition[1]
							 + rotationMovementPosition[2]*rotationMovementPosition[2]);
		rotationMovementPosition[0] *= a;
		rotationMovementPosition[1] *= a;
		rotationMovementPosition[2] *= a;

		float dx = rotationMovementPosition[0]-rotationStartPosition[0];
		float dy = rotationMovementPosition[1]-rotationStartPosition[1];
		float dz = rotationMovementPosition[2]-rotationStartPosition[2];

		// if the mouse has moved
		if (dx || dy || dz) {
			// compute rotation angle and axis
			float movementAngle = 90.0 * sqrt(dx*dx + dy*dy + dz*dz);
			if (movementAngle > 180) movementAngle -= 360.0;
			else if (movementAngle < -180) movementAngle += 360.0;

			float movementAxis[3];
			movementAxis[0] = rotationStartPosition[1]*rotationMovementPosition[2] - rotationStartPosition[2]*rotationMovementPosition[1];
			movementAxis[1] = rotationStartPosition[2]*rotationMovementPosition[0] - rotationStartPosition[0]*rotationMovementPosition[2];
			movementAxis[2] = rotationStartPosition[0]*rotationMovementPosition[1] - rotationStartPosition[1]*rotationMovementPosition[0];

			// compute quarternion from rotation angle and axis
			movementQuaternion.compute(movementAngle, movementAxis);

			// add movement quaternion with the quaternion that is the composite of previous rotations
			Quaternion tempQuat = compositeQuaternion+movementQuaternion;

			// compute the new rotation angle and axis
			tempQuat.getRotationAxisAndAngle(compositeAngle, compositeAxis);
		}
	}

	inline void trackballRecordMouseRelease(const float x, const float y, const float screenWidth, const float screenHeight) {
		rotationMouseDown = false;

		// mouse was released, so make the movement quaternion part of overall rotation quaternion
		compositeQuaternion += movementQuaternion;
		compositeQuaternion.getRotationAxisAndAngle(compositeAngle, compositeAxis);

		movementQuaternion.reset();
	}

	inline void doGlRotate() const { glRotatef(compositeAngle, compositeAxis[0], compositeAxis[1], compositeAxis[2]); }
	inline void reverseRotateX() { rotationSensitivity[0] *= -1; }
	inline void reverseRotateY() { rotationSensitivity[1] *= -1; }
	inline void setRotateSensitivityX(const float sensitivity) { rotationSensitivity[0] = (rotationSensitivity[0] < 0) ? -sensitivity : sensitivity; }
	inline void setRotateSensitivityY(const float sensitivity) { rotationSensitivity[1] = (rotationSensitivity[1] < 0) ? -sensitivity : sensitivity; }

	inline void getRotationMatrix(float rotationMatrix[16]) const { rotationMatrix = (compositeQuaternion + movementQuaternion).getRotationMatrix(); }
	inline void getRotationAxisAndAngle(float &angle, float axis[3]) const { (compositeQuaternion + movementQuaternion).getRotationAxisAndAngle(angle, axis); }

	// Scale Procedures
	inline void resetScale() {
		scaleMouseDown = false;
		scale = 1.0;
		scaleStartPosition[0] = scaleStartPosition[1] = 0;
	}

	inline void scaleRecordMouseDown(const float x, const float y, const float screenWidth, const float screenHeight) {
		scaleMouseDown = true;
		scaleStartPosition[0] = x;
		scaleStartPosition[1] = y;
	}

	inline void scaleRecordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		if (!scaleMouseDown) return;

		// compute new scale based on movement in y direction
		float oldScale = scale;
		scale /= (1.0+(y-scaleStartPosition[1])/(screenHeight*0.5));
		
		if (scale < 0) scale = oldScale;
		scaleStartPosition[0] = x;
		scaleStartPosition[1] = y;
	}

	inline void scaleRecordMouseRelease(const float x, const float y, const float screenWidth, const float screenHeight) { scaleMouseDown = false; }
	
	inline void doGlScale() const {
		glTranslatef(-translation[0], -translation[1], -translation[2]);
		glScalef(scale, scale, scale);
		glTranslatef(translation[0], translation[1], translation[2]);
	}

	inline void getScale(float &scaleVal) const { scaleVal = scale; }
	inline void setScale(const float scaleVal) { scale = scaleVal; }

	// Translation Procedures
	inline void resetTranslate() {
		translationStartPosition[0] = translationStartPosition[1] = 0;
		translation[0] = translation[1] = translation[2] = 0;
		translationSensitivity[0] = translationSensitivity[1] = 0.1;
		translationMouseDown = false;
	}

	inline void translateSetSensitivity(const float widthSensitivity, const float heightSensitivity) {
		translationSensitivity[0] = widthSensitivity * 0.1;
		translationSensitivity[1] = heightSensitivity * 0.1;
	}

	inline void translateRecordMouseDown(const float x, const float y, const float screenWidth, const float screenHeight) {
		translationMouseDown = true;
		translationStartPosition[0] = x;
		translationStartPosition[1] = y;
	}

	inline void translateRecordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		if (!translationMouseDown) return;

		// add to translation in direction that the mouse moved
		float dx = x - translationStartPosition[0];
		float dy = y - translationStartPosition[1];

		translation[0] += (dx/(screenWidth*scale*translationSensitivity[0]));
		translation[1] += (dy/(screenHeight*scale*translationSensitivity[1]));

		translationStartPosition[0] = x;
		translationStartPosition[1] = y;
	}

	inline void translateRecordMouseRelease(const float x, const float y, const float screenWidth, const float screenHeight) { translationMouseDown = false; }
	inline void doGlTranslate() const { glTranslatef(translation[0], translation[1], translation[2]); }

	inline void getTranslate(float &x, float &y, float &z) const { x = translation[0];  y = translation[1];  z = translation[2]; }
	inline void setTranslate(const float x, const float y, const float z) { translation[0] = x;  translation[1] = y;  translation[2] = z; }

	inline void doGlTransformations() const {
		doGlTranslate();
		doGlRotate();
		doGlScale();
	}

private:
	// Rotation Vars
	float compositeAngle;  // in degrees
	float compositeAxis[3];
	Quaternion compositeQuaternion;
	float rotationStartPosition[3];
	Quaternion movementQuaternion;
	float rotationSensitivity[2];
	bool rotationMouseDown;

	// Scale Vars
	float scaleStartPosition[2];
	float scale;
	bool scaleMouseDown;

	// Translation Vars
	float translationStartPosition[2];
	float translation[3];
	float translationSensitivity[2];
	bool translationMouseDown;
};

#endif