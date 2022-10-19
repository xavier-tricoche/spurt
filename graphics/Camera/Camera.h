#ifndef CAMERA_H_9376283
#define CAMERA_H_9376283

#include "graphics/GFXMath/Matrix4x4.h"
#include "graphics/GFXMath/VectorN.h"
#include "graphics/GFXMath/MathUtils.h"

using namespace MathUtils;

#if 0
#include <GL/gl.h>
#include "GL/glut.h"
#else
#include <glut.h>
#endif

class Camera
{
public:
	Camera() {
		set(Vector3f(0, 0, 10), Vector3f(0, 0, 0), Vector3f(0, 1, 0));
	}
	Camera(const Vector3f &location, const Vector3f &lookAtPos, const Vector3f &upVector) {
		set(location, lookAtPos, upVector);
	}
	Camera(const Vector3f &location, const float hRotation, const float vRotation) {
		set(location, hRotation, vRotation);
	}
	Camera(const Camera &camera) {
		*this = camera;
	}

	inline void set(const Vector3f &location, const float yaw, const float pitch) {
		reset();

		eye = originalEye = location;
		up.set(0, 1, 0);
		absoluteUp = up;
		left.set(-1, 0, 0);
		lookDir.set(0, 0, 1);

		rotateCameraAboutCenter(yaw, pitch);
	}

	inline void set(const Vector3f &location, const Vector3f &lookAtPos, const Vector3f &upVector) {
		reset();

		eye = originalEye = location;
		lookDir = (lookAtPos - eye).unit();
		absoluteUp = upVector.unit();

		left = lookDir.cross(absoluteUp);
		up = left.cross(lookDir);

		updateCamera();
	}

	inline void reset(bool alsoResetPerspective = true) {
		eye = originalEye;

		resetPan();
		resetDolly();
		resetRotate();
		if (alsoResetPerspective) resetPerspective();
	}
	inline void resetPan() {
		doPan = false;
		panSensitivity[0] = panSensitivity[1] = 10;
	}
	inline void resetDolly() {
		doDolly = false;
		dollySensitivity = 20;
	}
	inline void resetRotate() {
		doRotate = false;
		rotateSensitivity[0] = rotateSensitivity[1] = 5;
	}

	inline Camera& operator= (const Camera &camera) {
		reset();
		eye = originalEye = camera.getPosition();
		lookDir = camera.getLookAtDirection();
		left = camera.getLeftVector();
		up = camera.getUpVector();

		if (camera.isOrthographic()) setOrthographic(camera.getLeftClipPlane(), camera.getRightClipPlane(), camera.getBottomClipPlane(),
			        camera.getTopClipPlane(), camera.getNearClipPlane(), camera.getFarClipPlane());
		else setPerspective(camera.getLeftClipPlane(), camera.getRightClipPlane(), camera.getBottomClipPlane(),
			                    camera.getTopClipPlane(), camera.getNearClipPlane(), camera.getFarClipPlane());

		return *this;
	}

	inline bool operator==(Camera &camera) const {
		Vector3f cEye = camera.getPosition();
		Vector3f cLookDir = camera.getLookAtDirection();
		Vector3f cUpVec = camera.getUpVector();
		Vector3f cLeftVec = camera.getLeftVector();

		if (eye != cEye) return false;
		if (lookDir != cLookDir) return false;
		if (up != cUpVec) return false;
		if (left != cLeftVec) return false;

		return true;
	}
	inline bool operator!=(Camera &camera) const {
		return !(*this == camera);
	}

	inline void setGluLookAt() const {
		Vector3f lookAt = getLookAtPos();
		gluLookAt(eye.x, eye.y, eye.z, lookAt.x, lookAt.y, lookAt.z, up.x, up.y, up.z);
	}

    // inline void glMatrixMultiplyByCameraMatrix() {
    //  gluLookAt(cameraMat.pointer());
    // }
	inline Matrix4x4f getCameraMatrix() const {
		return cameraMat;
	}

	inline void recordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		rotateRecordMouseMotion(x, y, screenWidth, screenHeight);
		dollyRecordMouseMotion(x, y, screenWidth, screenHeight);
		panRecordMouseMotion(x, y, screenWidth, screenHeight);
	}

	inline void rotateRecordMouseDown(const float x, const float y, const float screenWidth, const float screenHeight) {
		doRotate = true;
		rotateStartPosition[0] = x;
		rotateStartPosition[1] = y;
	}
	inline void rotateRecordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		if (!doRotate) return;

		float dx = (rotateStartPosition[0] - x) / screenWidth * rotateSensitivity[0];
		float dy = (rotateStartPosition[1] - y) / screenHeight * rotateSensitivity[1];

		rotateCameraAboutCenter(dx, dy);

		rotateStartPosition[0] = x;
		rotateStartPosition[1] = y;
	}
	inline void rotateRecordMouseRelease(const float x, const float y, const float screenWidth, const float screenHeight) {
		doRotate = false;
	}

	inline void dollyRecordMouseDown(const float x, const float y, const float screenWidth, const float screenHeight) {
		doDolly = true;
		dollyStartPosition[0] = x;
		dollyStartPosition[1] = y;
	}
	inline void dollyRecordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		if (!doDolly) return;

		float dist = (dollyStartPosition[1] - y) * dollySensitivity / screenHeight;

		dollyCamera(dist);

		dollyStartPosition[0] = x;
		dollyStartPosition[1] = y;
	}
	inline void dollyRecordMouseRelease(const float x, const float y, const float screenWidth, const float screenHeight) {
		doDolly = false;
	}

	inline void panRecordMouseDown(const float x, const float y, const float screenWidth, const float screenHeight) {
		doPan = true;
		panStartPosition[0] = x;
		panStartPosition[1] = y;
	}
	inline void panRecordMouseMotion(const float x, const float y, const float screenWidth, const float screenHeight) {
		if (!doPan) return;

		float dx = (panStartPosition[0] - x) / screenWidth * panSensitivity[0];
		float dy = (panStartPosition[1] - y) / screenHeight * panSensitivity[1];

		panCamera(dx, dy);

		panStartPosition[0] = x;
		panStartPosition[1] = y;
	}
	inline void panRecordMouseRelease(const float x, const float y, const float screenWidth, const float screenHeight) {
		doPan = false;
	}

	inline void rotateCameraAboutCenter(const float yaw, const float pitch) {
		Vector3f lookAt = getLookAtPos();
		lookAt += left * -yaw;
		lookAt += up * -pitch;
		lookDir = eye - lookAt;
		updateCamera();
	}

	inline void panCamera(const float h, const float v) {
		eye += left * h;
		eye += up * v;
		/*Vector3f lookAt = getLookAtPos();
		lookAt += left*h;
		lookAt += up*v;
		lookDir = eye-lookAt;*/
		updateCamera();
	}

	inline void dollyCamera(const float d) {
		if (orthographic) {
			float left = getLeftClipPlane();
			float right = getRightClipPlane();
			float bottom = getBottomClipPlane();
			float top = getTopClipPlane();

			float ratio = (right - left) / (top - bottom);

			left -= d;
			right += d;
			bottom -= d/ratio;
			top += d/ratio;
			
			if (left >= right || bottom >= top) return;

			setOrthographic(left, right, bottom, top, getNearClipPlane(), getFarClipPlane());
		}
		else {
			eye += lookDir.unit() * d;
		}
		updateCamera();
	}

	inline void lookUp(const float radians) {
		Matrix4x4f rotate;
		rotate.setAsRotate(left, radiansToDegrees(radians));
		up = rotate * up;
		lookDir = rotate * lookDir;
		updateCamera();
	}

	inline void lookLeft(const float radians) {
		Matrix4x4f rotate;
		rotate.setAsRotate(absoluteUp, radiansToDegrees(radians));
		up = rotate * up;
		left = rotate * left;
		lookDir = rotate * lookDir;
		updateCamera();
	}

	inline void rollCamera(const float radians) {
		Matrix4x4f rotate;
		rotate.setAsRotate(up, radiansToDegrees(radians));
		left = rotate * left;
		lookDir = rotate * lookDir;
		updateCamera();
	}

	inline Vector3f getPosition() const {
		return eye;
	}
	inline Vector3f getLookAtPos() const {
		return eye + lookDir.unit();
	}
	inline Vector3f getLookAtDirection() const {
		return lookDir.unit();
	}
	inline Vector3f getUpVector() const {
		return up.unit();
	}
	inline Vector3f getLeftVector() const {
		return left.unit();
	}


	inline void reversePanX() {
		panSensitivity[0] *= -1;
	}
	inline void reversePanY() {
		panSensitivity[1] *= -1;
	}
	inline void setPanSensitivityX(const float sensitivity) {
		panSensitivity[0] = (panSensitivity[0] < 0) ? -sensitivity : sensitivity;
	}
	inline void setPanSensitivityY(const float sensitivity) {
		panSensitivity[1] = (panSensitivity[1] < 0) ? -sensitivity : sensitivity;
	}

	inline void reverseDolly() {
		dollySensitivity *= -1;
	}
	inline void setDollySensitivity(const float sensitivity) {
		dollySensitivity = (dollySensitivity < 0) ? -sensitivity : sensitivity;
	}

	inline void reverseRotateX() {
		rotateSensitivity[0] *= -1;
	}
	inline void reverseRotateY() {
		rotateSensitivity[1] *= -1;
	}
	inline void setRotateSensitivityX(const float sensitivity) {
		rotateSensitivity[0] = (rotateSensitivity[0] < 0) ? -sensitivity : sensitivity;
	}
	inline void setRotateSensitivityY(const float sensitivity) {
		rotateSensitivity[1] = (rotateSensitivity[1] < 0) ? -sensitivity : sensitivity;
	}

	inline void resetPerspective() {
		setPerspective(60.0, 1, 0.1, 100.0);
	}

	inline void setOrthographic(const double left, const double right, const double bottom, const double top, const double nearClip, const double farClip) {
		orthographic = true;
		setPerspectiveProperties(left, right, bottom, top, nearClip, farClip);
	}

	inline void setPerspective(const double fovy_degrees, const double aspect, const double nearClip, const double farClip) {
		double rad = MathUtils::degreesToRadians(fovy_degrees);
		double top = nearClip * tanf(rad * 0.5);
		double bottom = -top;
		setPerspective(aspect*bottom, aspect*top, bottom, top, nearClip, farClip);
	}

	inline void setPerspective(const double left, const double right, const double bottom, const double top, const double nearClip, const double farClip) {
		orthographic = false;
		setPerspectiveProperties(left, right, bottom, top, nearClip, farClip);
	}

	inline double getLeftClipPlane() const {
		return _left;
	}
	inline double getRightClipPlane() const {
		return _right;
	}
	inline double getBottomClipPlane() const {
		return _bottom;
	}
	inline double getTopClipPlane() const {
		return _top;
	}
	inline double getNearClipPlane() const {
		return _near;
	}
	inline double getFarClipPlane() const {
		return _far;
	}
	inline double getFOV() const {
		return MathUtils::radiansToDegrees(2.0 * atan(_top / _near));
	}
	inline double getAspect() const {
		return _left / _bottom;
	}

	inline bool isOrthographic() const {
		return orthographic;
	}

	inline void setGLPerspective() const {
		if (orthographic) glOrtho(_left, _right, _bottom, _top, _near, _far);
		else glFrustum(_left, _right, _bottom, _top, _near, _far);
	}
	inline Matrix4x4f getPerspectiveMatrix() const {
		return perspectiveMat;
	}

	inline void getLeftFrustumPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		lowerLeft = leftClipPlane[0];
		topLeft = leftClipPlane[2];
		lowerRight = leftClipPlane[1];
		topRight = leftClipPlane[3];
	}

	inline void getRightFrustumPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		lowerLeft = rightClipPlane[0];
		topLeft = rightClipPlane[2];
		lowerRight = rightClipPlane[1];
		topRight = rightClipPlane[3];
	}

	inline void getBottomFrustumPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		lowerLeft = bottomClipPlane[0];
		topLeft = bottomClipPlane[2];
		lowerRight = bottomClipPlane[1];
		topRight = bottomClipPlane[3];
	}

	inline void getTopFrustumPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		lowerLeft = topClipPlane[0];
		topLeft = topClipPlane[2];
		lowerRight = topClipPlane[1];
		topRight = topClipPlane[3];
	}

	inline void getNearFrustumPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		lowerLeft = nearClipPlane[0];
		topLeft = nearClipPlane[2];
		lowerRight = nearClipPlane[1];
		topRight = nearClipPlane[3];
	}

	inline void getFarFrustumPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		lowerLeft = farClipPlane[0];
		topLeft = farClipPlane[2];
		lowerRight = farClipPlane[1];
		topRight = farClipPlane[3];
	}

	inline Vector3f getRayDirection(const Vector3f &pos) const {
		if (orthographic) return getLookAtDirection().unit();
		else return pos -getPosition();
	}

private:
	inline void updateCamera() {
		left.normalize();
		up.normalize();
		cameraMat.setAsCamera(getPosition(), getLookAtPos(), getUpVector());
		updateViewFrustum();
	}

	Vector3f eye, originalEye, lookDir, up, left, absoluteUp;
	Matrix4x4f cameraMat;

	bool doPan, doDolly, doRotate;

	float panStartPosition[2];
	float panSensitivity[2];

	float dollyStartPosition[2];
	float dollySensitivity;

	float rotateStartPosition[2];
	float rotateSensitivity[2];



	inline void setPerspectiveProperties(double left, double right, double bottom, double top, double nearClip, double farClip) {
		_left = left;
		_right = right;
		_bottom = bottom;
		_top = top;
		_near = nearClip;
		_far = farClip;
		updatePerspective();
		updateViewFrustum();
	}

	inline void updatePerspective() {
		if (orthographic) perspectiveMat.setAsOrthographic(_left, _right, _bottom, _top, _near, _far);
		else perspectiveMat.setAsPerspective(_left, _right, _bottom, _top, _near, _far);
	}

	bool orthographic;
	double _left, _right, _bottom, _top, _near, _far;
	Matrix4x4f perspectiveMat;


	inline void updateViewFrustum() {
		cpMatrix = cameraMat * perspectiveMat;
		invCpMatrix = cpMatrix.getInverse();

		Vector3f pos = getPosition();
		Vector3f lookAt = getLookAtDirection();
		Vector3f up = getUpVector();
		Vector3f left = getLeftVector();

		Vector3f nearCenter = pos + lookAt * _near;
		Vector3f farCenter = pos + lookAt * _far;

		double aspect = getAspect();
		double tang = _top / _near;

		double nearHeight = _near * tang;
		double nearWidth = nearHeight * aspect;

		double farHeight, farWidth;
		if (orthographic) {
			farHeight = nearHeight;
			farWidth = nearWidth;
		}
		else {
			farHeight = _far * tang;
			farWidth = farHeight * aspect;
		}

		// compute the 4 corners of the frustum on the near plane
		Vector3f ntl = nearCenter + up * nearHeight - left * nearWidth;
		Vector3f ntr = nearCenter + up * nearHeight + left * nearWidth;
		Vector3f nbl = nearCenter - up * nearHeight - left * nearWidth;
		Vector3f nbr = nearCenter - up * nearHeight + left * nearWidth;

		// compute the 4 corners of the frustum on the far plane
		Vector3f ftl = farCenter + up * farHeight - left * farWidth;
		Vector3f ftr = farCenter + up * farHeight + left * farWidth;
		Vector3f fbl = farCenter - up * farHeight - left * farWidth;
		Vector3f fbr = farCenter - up * farHeight + left * farWidth;

		leftClipPlane[0] = nbl;
		leftClipPlane[1] = ntl;
		leftClipPlane[2] = fbl;
		leftClipPlane[3] = ftl;
		rightClipPlane[0] = nbr;
		rightClipPlane[1] = ntr;
		rightClipPlane[2] = fbr;
		rightClipPlane[3] = ftr;

		bottomClipPlane[0] = nbl;
		bottomClipPlane[1] = nbr;
		bottomClipPlane[2] = fbl;
		bottomClipPlane[3] = fbr;
		topClipPlane[0] = ntl;
		topClipPlane[1] = ntr;
		topClipPlane[2] = ftl;
		topClipPlane[3] = ftr;

		nearClipPlane[0] = nbl;
		nearClipPlane[1] = nbr;
		nearClipPlane[2] = ntl;
		nearClipPlane[3] = ntr;
		farClipPlane[0] = fbl;
		farClipPlane[1] = fbr;
		farClipPlane[2] = ftl;
		farClipPlane[3] = ftr;
	}


	Matrix4x4f cpMatrix, invCpMatrix;
	Vector3f rightClipPlane[4], leftClipPlane[4], topClipPlane[4], bottomClipPlane[4], nearClipPlane[4], farClipPlane[4];
};

#endif








