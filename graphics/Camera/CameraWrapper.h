#ifndef CAMERA_WRAPPER_H
#define CAMERA_WRAPPER_H

#include "Camera.h"
#include "VirtualTrackball.h"
#include "graphics/GFXMath/VectorN.h"
#include "graphics/GFXMath/Matrix4x4.h"

#include <vector>
using namespace std;

// to do: region selection, object picking, saving a camera view, scroll wheel, recording movement path, get view frustum
//    for region selection pass in a vector for where points are recorded and set a flag when the button is released
//      set a flag when the mouse is released indicating the user should use those points for something
//    similarly with object picking, have the user pass in a set of models... put probably do not want to have #include "models" in this

// to do: push/pop matrix
//        CameraWrapper_Misc for creating a camera, glui controls
//        my own glUnproject (look at man pages for implentation)

// to do: RotateAboutLookAtCenter seems to only rotate about (0,0,0)
//           also does not rotate right if originally centered at (0,0,0) then you pan to another location
//          it needs to be translated somehow before doing rotate, maybe the camera matrix shouldn't have the translate
//           component added until after the rotate

// to do: be able to set the viewport end pixel, right now setting the start.x = x and width = windowWidth-2x
//         results in the end pixel being outside of the window because of the GLUT resize
//         may need to record the window size, which would also enable the ability to do resizing options
//          (such as expand start.x to move with the window resize)

class CameraWrapper {
public:
	CameraWrapper() {
		setLeftMouse_RotateAboutLookAtCenter();
		setRightMouse_Dolly();
		setMiddleMouse_Pan();
		reset();

		regionLineColor.set(1,1,1);
		zoomRegionPts.clear();
		selectRegionPts.clear();

		setRotationAngleStep(5);
		setRotationAxis(Vector3f(0,1,0));
		doNotModifyRotationAnimationByMouse();

		limitMouseMovementToViewport = false;
		limitMouseClickToViewport = true;
	}

	inline CameraWrapper& operator= (const CameraWrapper &cw) {
		reset();
		trackball = cw.trackball;
		camera = cw.camera;
		setCamera(cw.getCameraPos(), cw.getCameraLookAtPos(), cw.getCameraUpVector());
		return *this;
	}

	void reset() {
		trackball.reset();
		camera.reset();
		resetRotationAnimation();
		setCamera(Vector3f(0, 0, 10), Vector3f(0,0,0), Vector3f(0,1,0));
	}

	inline void setLeftMouse_RotateAboutLookAtCenter() { setButton_RotateAboutLookAtCenter(leftButtonAction); }
	inline void setLeftMouse_RotateAboutCameraCenter() { setButton_RotateAboutCameraCenter(leftButtonAction); }
	inline void setLeftMouse_Pan() { setButton_Pan(leftButtonAction); }
	inline void setLeftMouse_Dolly() { setButton_Dolly(leftButtonAction); }
	inline void setLeftMouse_SelectRegion() { setButton_SelectRegion(leftButtonAction); }
	inline void setLeftMouse_ZoomRegion() { setButton_ZoomRegion(leftButtonAction); }

	inline void setRightMouse_RotateAboutLookAtCenter() { setButton_RotateAboutLookAtCenter(rightButtonAction); }
	inline void setRightMouse_RotateAboutCameraCenter() { setButton_RotateAboutCameraCenter(rightButtonAction); }
	inline void setRightMouse_Pan() { setButton_Pan(rightButtonAction); }
	inline void setRightMouse_Dolly() { setButton_Dolly(rightButtonAction); }
	inline void setRightMouse_SelectRegion() { setButton_SelectRegion(rightButtonAction); }
	inline void setRightMouse_ZoomRegion() { setButton_ZoomRegion(rightButtonAction); }

	inline void setMiddleMouse_RotateAboutLookAtCenter() { setButton_RotateAboutLookAtCenter(middleButtonAction); }
	inline void setMiddleMouse_RotateAboutCameraCenter() { setButton_RotateAboutCameraCenter(middleButtonAction); }
	inline void setMiddleMouse_Pan() { setButton_Pan(middleButtonAction); }
	inline void setMiddleMouse_Dolly() { setButton_Dolly(middleButtonAction); }
	inline void setMiddleMouse_SelectRegion() { setButton_SelectRegion(middleButtonAction); }
	inline void setMiddleMouse_ZoomRegion() { setButton_ZoomRegion(middleButtonAction); }

	inline void recordLeftMouseDown(const int x, const int y) { recordButtonDown(leftButtonAction, x, y); }
	inline void recordLeftMouseUp(const int x, const int y) { recordButtonUp(leftButtonAction, x, y); }
	inline void recordRightMouseDown(const int x, const int y) { recordButtonDown(rightButtonAction, x, y); }
	inline void recordRightMouseUp(const int x, const int y) { recordButtonUp(rightButtonAction, x, y); }
	inline void recordMiddleMouseDown(const int x, const int y) { recordButtonDown(middleButtonAction, x, y); }
	inline void recordMiddleMouseUp(const int x, const int y) { recordButtonUp(middleButtonAction, x, y); }

	inline void reversePanX() { camera.reversePanX(); }
	inline void reversePanY() { camera.reversePanY(); }
	inline void reverseDolly() { camera.reverseDolly(); }
	inline void reverseRotateAboutCameraCenterX() { camera.reverseRotateX(); }
	inline void reverseRotateAboutCameraCenterY() { camera.reverseRotateY(); }
	inline void reverseRotateAboutLookAtCenterX() { trackball.reverseRotateX(); }
	inline void reverseRotateAboutLookAtCenterY() { trackball.reverseRotateY(); }

	inline void setPanSensitivityX(const float sensitivity) { camera.setPanSensitivityX(sensitivity); }
	inline void setPanSensitivityY(const float sensitivity) { camera.setPanSensitivityY(sensitivity); }
	inline void setDollySensitivity(const float sensitivity) { camera.setDollySensitivity(sensitivity); }
	inline void setRotateAboutCameraCenterSensitivityX(const float sensitivity) { camera.setRotateSensitivityX(sensitivity); }
	inline void setRotateAboutCameraCenterSensitivityY(const float sensitivity) { camera.setRotateSensitivityY(sensitivity); }
	inline void setRotateAboutLookAtCenterSensitivityX(const float sensitivity) { trackball.setRotateSensitivityX(sensitivity); }
	inline void setRotateAboutLookAtCenterSensitivityY(const float sensitivity) { trackball.setRotateSensitivityY(sensitivity); }

	inline void setCamera(const Vector3f &location, const Vector3f &lookAtPos, const Vector3f &upVector) {
		trackball.reset();
		camera.set(location, lookAtPos, upVector);
		rotatedCameraPos = camera.getPosition();
		rotatedLookAtPos = camera.getLookAtPos();
		rotatedUpVector = camera.getUpVector();
		update();
	}

	inline void scaleCameraToObject(const Vector3f &objectCenter, const Vector3f &objectBoundingBoxDim, const Vector3f &upVector,
									const float perspectiveFovyDegrees, const float perspectiveAspect, const float relativeSensitivity = 2.5) {
		double maxDim = objectBoundingBoxDim.maxValue();
		setCamera(objectCenter+Vector3f(0,0,maxDim*1.25), objectCenter, upVector);
		setPerspective(perspectiveFovyDegrees, perspectiveAspect, 0.001, maxDim*2.5 + objectCenter.z + 10);

		float sensitivityScale = maxDim*relativeSensitivity;
		setDollySensitivity(sensitivityScale);
		setPanSensitivityX(sensitivityScale);
		setPanSensitivityY(sensitivityScale);
	}

	inline void resizeViewport(const int newWidth, const int newHeight) {
		viewportWidth = newWidth;
		viewportHeight = newHeight;

		invViewportWidth = 1.0/(float)viewportWidth;
		invViewportHeight = 1.0/(float)viewportHeight;

		glViewport(viewportStartX, viewportStartY, viewportWidth, viewportHeight);
	}

	inline void setViewport(const int startX, const int startY, const int width, const int height) {
		viewportStartX = startX;
		viewportStartY = startY;
		resizeViewport(width, height);
	}

	Vector3f getCameraPos() const { return rotatedCameraPos; }
	Vector3f getCameraLookAtPos() const { return rotatedLookAtPos; }
	Vector3f getCameraLookAtDir() const { return (getCameraLookAtPos()-getCameraPos()).unit(); }
	Vector3f getCameraUpVector() const { return rotatedUpVector; }

	inline void getFrustrumBoundingBox(Vector3f &bbMin, Vector3f &bbMax) const {
		Vector3f pts[8];
		getFrustrum(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7]);
		bbMin = bbMax = pts[0];
		for (unsigned int i=1; i<8; i++) {
			bbMin = bbMin.minVector(pts[i]);
			bbMax = bbMax.maxVector(pts[i]);
		}
	}

	inline void getFrustrum(Vector3f &pt000, Vector3f &pt001, Vector3f &pt010, Vector3f &pt011,
							Vector3f &pt100, Vector3f &pt101, Vector3f &pt110, Vector3f &pt111) const {
		getNearPlane(pt000, pt100, pt010, pt110);
		getFarPlane(pt001, pt101, pt011, pt111);
	}

	inline void getLeftPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		camera.getLeftFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);
		rotatePointsByTrackballAndAnimator(lowerLeft, lowerRight, topLeft, topRight);
	}
	inline void getRightPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		camera.getRightFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);
		rotatePointsByTrackballAndAnimator(lowerLeft, lowerRight, topLeft, topRight);
	}
	inline void getBottomPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		camera.getBottomFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);
		rotatePointsByTrackballAndAnimator(lowerLeft, lowerRight, topLeft, topRight);
	}
	inline void getTopPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		camera.getTopFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);
		rotatePointsByTrackballAndAnimator(lowerLeft, lowerRight, topLeft, topRight);
	}
	inline void getNearPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		camera.getNearFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);
		rotatePointsByTrackballAndAnimator(lowerLeft, lowerRight, topLeft, topRight);
	}
	inline void getFarPlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const {
		camera.getFarFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);
		rotatePointsByTrackballAndAnimator(lowerLeft, lowerRight, topLeft, topRight);
	}
	inline void getCameraImagePlane(Vector3f &lowerLeft, Vector3f &lowerRight, Vector3f &topLeft, Vector3f &topRight) const
		{ getNearPlane(lowerLeft, lowerRight, topLeft, topRight); }

	inline Vector3f getPixelPosition(const float pixelX, const float pixelY) const {
		float tx = pixelX*invViewportWidth;
		float ty = pixelY*invViewportHeight;

		Vector3f lowerLeft, lowerRight, topLeft, topRight;
		getCameraImagePlane(lowerLeft, lowerRight, topLeft, topRight);

		return getRotatedPointByTrackballAndAnimator(MathUtils::bilinearLerp(lowerLeft, topLeft, lowerRight, topRight, tx, ty));
	}

	inline Vector3f getRayDirection(const float pixelX, const float pixelY) const {
		float tx = 1.0 - pixelX*invViewportWidth;
		float ty = 1.0 - pixelY*invViewportHeight;

		Vector3f lowerLeft, lowerRight, topLeft, topRight;
		camera.getNearFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);

		return getRotatedPointByTrackballAndAnimator(camera.getRayDirection(MathUtils::bilinearLerp(lowerLeft, topLeft, lowerRight, topRight, tx, ty)));
	}

	inline void getRay(const float pixelX, const float pixelY, Vector3f &rayOrigin, Vector3f &rayDirection) const {
		float tx = 1.0 - pixelX*invViewportWidth;
		float ty = 1.0 - pixelY*invViewportHeight;

		Vector3f lowerLeft, lowerRight, topLeft, topRight;
		camera.getNearFrustumPlane(lowerLeft, lowerRight, topLeft, topRight);

		rayOrigin = MathUtils::bilinearLerp(lowerLeft, topLeft, lowerRight, topRight, tx, ty);
		rayDirection = camera.getRayDirection(rayOrigin);

		rayOrigin = getRotatedPointByTrackballAndAnimator(rayOrigin);
		rayDirection = getRotatedPointByTrackballAndAnimator(rayDirection);
	}

	inline void getActualPixelCoordinates(int &x, int &y) const {
		bool invertPositionX = false;
		bool invertPositionY = false;

		#ifdef GLUT_API_VERSION
			invertPositionY = true;
		#endif

		if (invertPositionX) x = viewportWidth-x;
		if (invertPositionY) y = viewportHeight-y;
	}

	void recordMouseMotion(const int pixelX, const int pixelY) {
		int x = pixelX, y = pixelY;
		getActualPixelCoordinates(x, y);

		if (limitMouseMovementToViewport) {
			if (x < viewportStartX || x >= viewportStartX+viewportWidth || y < viewportStartY || y >= viewportStartY+viewportHeight) {
				std::cerr << "mouse cursor is out of viewport. ignoring\n";
				return;
			}
		}

		trackball.trackballRecordMouseMotion(x, y, viewportWidth, viewportHeight);
		camera.recordMouseMotion(x, y, viewportWidth, viewportHeight);

		if (zoomRegionPts.size() > 0 || selectRegionPts.size() > 0) {
			Vector3f unPos = unProject(x, y);
			unPos.z = 0;
			if (zoomRegionPts.size() > 0) {
				zoomRegionPts[2] = unPos;
				zoomRegionPts[1].set(zoomRegionPts[0].x, zoomRegionPts[2].y, 0);
				zoomRegionPts[3].set(zoomRegionPts[2].x, zoomRegionPts[0].y, 0);
			}
			else selectRegionPts.push_back(unPos);
		}

		update();
	}

	inline void rotateAboutLookAtCenter(const int x, const int y) { recordIncrementalMovement(MOUSE_ROTATE_ABOUT_LOOKAT, x, y); }
	//inline void rotateAboutCameraCenter(const int x, const int y) { recordIncrementalMovement(MOUSE_ROTATE_ABOUT_CAMERA, x, y); } doesn't seem to work
	inline void pan(const int x, const int y) { recordIncrementalMovement(MOUSE_PAN, x, y); }
	inline void dolly(const int x) { recordIncrementalMovement(MOUSE_DOLLY, x, x); }

	Vector3f unProject(const int x, const int y, bool useObjectPosition = false) const {
		GLint viewport[4];
		GLdouble modelview[16];
		GLdouble projection[16];

		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetIntegerv(GL_VIEWPORT, viewport);

		// get z position
		GLfloat z;
		if (useObjectPosition) glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
		else {
			GLdouble winX, winY, winZ;
			gluProject(0, 0, 0, modelview, projection, viewport, &winX, &winY, &winZ);
			z = winZ;
		}

		GLdouble posX, posY, posZ;
		gluUnProject(x, y, z, modelview, projection, viewport, &posX, &posY, &posZ);

		// gluUnProject can return bad values if at camera center (possibly other times too)
		if (posX != posX) posX = 0;
		if (posY != posY) posY = 0;
		if (posZ != posZ) posZ = 0;

		return Vector3f(posX, posY, posZ);
	}

	Vector3f project(const Vector3f &pt) const {
		Vector4f pixel = gluProjectMatrix*Vector4f(pt.x, pt.y, pt.z, 1);
		return Vector3f((float)viewportStartX +  (float)viewportWidth*((pixel.x/pixel.w)+1.0)*0.5,
						(float)viewportStartY + (float)viewportHeight*((pixel.y/pixel.w)+1.0)*0.5,
						((pixel.z/pixel.w)+1.0)*0.5);
	}

	void setProjectionMatrix() {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		camera.setGLPerspective();
	}

	void setModelviewMatrix() {
		Vector3f pos = getCameraPos();
		Vector3f lookAt = getCameraLookAtPos();
		Vector3f up = getCameraUpVector();

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(pos.x, pos.y, pos.z, lookAt.x, lookAt.y, lookAt.z, up.x, up.y, up.z);
	}

	Vector3f regionLineColor;
	std::vector<Vector3f> zoomRegionPts, selectRegionPts;

	void renderSelectRegion() {
		if (zoomRegionPts.size() == 0 && selectRegionPts.size() == 0) return;

		glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);
		glLineWidth(2.0);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glColor3fv(regionLineColor.data);
		glBegin(GL_LINE_LOOP);
		for (unsigned int i=0; i<zoomRegionPts.size(); i++) glVertex3fv(zoomRegionPts[i].data);
		glEnd();

		glColor3fv(regionLineColor.data);
		glBegin(GL_LINE_LOOP);
		for (unsigned int i=0; i<selectRegionPts.size(); i++) glVertex3fv(selectRegionPts[i].data);
		glEnd();

		glPopAttrib();
	}

	inline void resetRotationAnimation() { rotationAngle = 0; }
	inline void setRotationAngleStep(const float degrees) { rotationAngleStep = degrees; }
	inline void setRotationAxis(Vector3f axis) { rotationAxis = axis; }
	inline void modifyRotationAnimationByMouse() { modifyAxisByTrackball = true; }
	inline void doNotModifyRotationAnimationByMouse() { modifyAxisByTrackball = false; }
	inline void stepRotationAnimation() {
		rotationAngle += rotationAngleStep;
		if (rotationAngle > 360) rotationAngle -= 360;
		update();
	}

	inline void setOrthographic(const double left, const double right, const double bottom, const double top, const double nearClip, const double farClip)
		{ camera.setOrthographic(left, right, bottom, top, nearClip, farClip);  update(); }
	inline void setPerspective(const double left, const double right, const double bottom, const double top, const double nearClip, const double farClip)
		{ camera.setPerspective(left, right, bottom, top, nearClip, farClip);  update(); }
	inline void setPerspective(const double fovy_degrees, const double aspect, const double nearClip, const double farClip)
		{ camera.setPerspective(fovy_degrees, aspect, nearClip, farClip);  update(); }

	inline double getLeftClipPlane() const { return camera.getLeftClipPlane(); }
	inline double getRightClipPlane() const { return camera.getRightClipPlane(); }
	inline double getBottomClipPlane() const { return camera.getBottomClipPlane(); }
	inline double getTopClipPlane() const { return camera.getTopClipPlane(); }
	inline double getNearZ() const { return camera.getNearClipPlane(); }
	inline double getFarZ() const { return camera.getFarClipPlane(); }

	inline Matrix4x4f getInverseCompositeRotationMatrix() { return invCompositeRotationMatrix; }

	inline void markChanged() { cameraHasChanged = true; }
	inline void markUnchanged() { cameraHasChanged = false; }
	inline bool hasCameraChanged() { return cameraHasChanged; }

private:
	int leftButtonAction, rightButtonAction, middleButtonAction;
	enum { MOUSE_ROTATE_ABOUT_LOOKAT, MOUSE_ROTATE_ABOUT_CAMERA, MOUSE_PAN, MOUSE_DOLLY, MOUSE_REGION_SELECT, MOUSE_REGION_ZOOM };

	inline void setButton_RotateAboutLookAtCenter(int &button) { button = MOUSE_ROTATE_ABOUT_LOOKAT; }
	inline void setButton_RotateAboutCameraCenter(int &button) { button = MOUSE_ROTATE_ABOUT_CAMERA; }
	inline void setButton_Pan(int &button) { button = MOUSE_PAN; }
	inline void setButton_Dolly(int &button) { button = MOUSE_DOLLY; }
	inline void setButton_SelectRegion(int &button) { button = MOUSE_REGION_SELECT; }
	inline void setButton_ZoomRegion(int &button) { button = MOUSE_REGION_ZOOM; }

	void recordButtonDown(const int button, const int pixelX, const int pixelY) {
		int x = pixelX, y = pixelY;
		getActualPixelCoordinates(x, y);

		if (limitMouseClickToViewport) {
			if (x < viewportStartX || x >= viewportStartX+viewportWidth || y < viewportStartY || y >= viewportStartY+viewportHeight)
				return;
		}

		if (button == MOUSE_ROTATE_ABOUT_LOOKAT) trackball.trackballRecordMouseDown(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_ROTATE_ABOUT_CAMERA) camera.rotateRecordMouseDown(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_PAN) camera.panRecordMouseDown(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_DOLLY) camera.dollyRecordMouseDown(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_REGION_ZOOM) {
			Vector3f unPos = unProject(x, y);
			for (int i=0; i<4; i++) zoomRegionPts.push_back(Vector3f(unPos.x, unPos.y, 0));
		}
		if (button == MOUSE_REGION_SELECT) selectRegionPts.push_back(unProject(x, y));
	}
	void recordButtonUp(const int button, const int pixelX, const int pixelY) {
		int x = pixelX, y = pixelY;
		getActualPixelCoordinates(x, y);

		if (button == MOUSE_ROTATE_ABOUT_LOOKAT) trackball.trackballRecordMouseRelease(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_ROTATE_ABOUT_CAMERA) camera.rotateRecordMouseRelease(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_PAN) camera.panRecordMouseRelease(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_DOLLY) camera.dollyRecordMouseRelease(x, y, viewportWidth, viewportHeight);
		if (button == MOUSE_REGION_ZOOM) {
			Vector3f unPos = unProject(x, y);

			Vector3f diff = unPos-zoomRegionPts[0];
			diff.z = 0;
			diff.absoluteValue();

			if (diff.magnitudeSquared() > 0.01) {
				float max = diff.x>diff.y ? diff.x : diff.y;
				Vector3f average = (unPos+zoomRegionPts[0]) * 0.5;
				average.z = 0;

				Vector3f prevPos = camera.getPosition();

				printf("to do:  do this with the camera and not the trackball\n");
				float oldX, oldY, oldZ, scale;
				trackball.getTranslate(oldX, oldY, oldZ);
				trackball.getScale(scale);

				trackball.setTranslate(oldX-(average.x/scale), oldY-(average.y/scale), 0);

				Vector3f unPos1 = unProject(0, 0);
				Vector3f unPos2 = unProject(viewportWidth, viewportHeight);

				diff = unPos1-unPos2;
				diff.absoluteValue();
				float minDiff = diff.x<diff.y ? diff.x : diff.y;

				trackball.setScale(scale * minDiff / max);
			}

			zoomRegionPts.clear();
		}
		if (button == MOUSE_REGION_SELECT) {
			selectRegionPts.clear();
		}
	}

	inline void recordIncrementalMovement(const int button, const int x, const int y) {
		if (-x > viewportWidth/2 || -y > viewportHeight/2)
			printf("Warning: incremental movement is too much (CameraWrapper.recordIncrementalMovement)\n");

		int startX = viewportWidth/2, startY = viewportHeight/2;
		int endX   = startX + x,      endY   = startY + y;

		getActualPixelCoordinates(startX, startY);
		getActualPixelCoordinates(endX, endY);

		recordButtonDown(button, startX, startY);
		recordMouseMotion(endX, endY);
		recordButtonUp(button, endX, endY);
	}

	inline void update() {
		// update the camera's position, lookAt, etc. using the trackball's transformations and animators
		updateCompositeRotationMatrix();

		rotatedCameraPos = camera.getPosition();
		rotatedLookAtPos = camera.getLookAtPos();
		rotatedUpVector = camera.getUpVector();

		Matrix4x4f cameraMat;
		cameraMat.setAsCamera(rotatedCameraPos, rotatedLookAtPos, rotatedUpVector);
		cameraMat = cameraMat * compositeRotationMatrix;

		cameraMat.getCameraProperties(rotatedCameraPos, rotatedLookAtPos, rotatedUpVector);

		gluProjectMatrix = camera.getPerspectiveMatrix()*cameraMat;

		markChanged();
	}

	inline void updateCompositeRotationMatrix() {
		float angle, axis[3];
		trackball.getRotationAxisAndAngle(angle, axis);

		Matrix4x4f rotMat, animateMat;
		rotMat.setAsRotate(axis[0], axis[1], axis[2], angle);
		animateMat.setAsRotate(rotationAxis, rotationAngle);

		Vector3f pos = camera.getLookAtPos();
		Matrix4x4f translateMat;
		translateMat.setAsTranslate(pos.x, pos.y, pos.z);

		if (modifyAxisByTrackball) compositeRotationMatrix = rotMat * animateMat;
		else compositeRotationMatrix = animateMat * rotMat;

		invCompositeRotationMatrix = compositeRotationMatrix.getInverse();
	}

	inline void rotatePointsByTrackballAndAnimator(Vector3f &pt1, Vector3f &pt2, Vector3f &pt3, Vector3f &pt4) const {
		pt1 = getRotatedPointByTrackballAndAnimator(pt1);		pt2 = getRotatedPointByTrackballAndAnimator(pt2);
		pt3 = getRotatedPointByTrackballAndAnimator(pt3);		pt4 = getRotatedPointByTrackballAndAnimator(pt4);
	}

	inline Vector3f getRotatedPointByTrackballAndAnimator(const Vector3f &pt) const { return invCompositeRotationMatrix * pt; }

	Camera camera;
	VirtualTrackball trackball;
	bool cameraHasChanged;

	Vector3f rotatedCameraPos;
	Vector3f rotatedLookAtPos;
	Vector3f rotatedUpVector;

	Vector3f rotationAxis;
	float rotationAngle, rotationAngleStep;
	bool modifyAxisByTrackball;

	Matrix4x4f compositeRotationMatrix, invCompositeRotationMatrix, gluProjectMatrix;

	int viewportStartX, viewportStartY, viewportWidth, viewportHeight;
	float invViewportWidth, invViewportHeight;
	bool limitMouseMovementToViewport, limitMouseClickToViewport;
};

// to do:  add this old code to the camera
	//else if (selectRegion || deselectRegion) {
	//			y = height - y - 1;
	//			double unX, unY, unZ;
	//			unProject(x, y, 0, unX, unY, unZ);
	//
	//			Vector3f newPt = Vector3f(unX, unY, 0);
	//			Vector3f lastPt = regionPts[regionPts.size()-1];

	//			float dx = fabs(newPt.x-lastPt.x);
	//			float dy = fabs(newPt.y-lastPt.y);

	//			if (dx > 0.01 || dy > 0.01) regionPts.push_back(newPt);

	//			// figure out what points in the data set are in the region
	//			int xDim = visManagers[dataSetToModify].getXDim();
	//			int yDim = visManagers[dataSetToModify].getYDim();
	//
	//			Vector3f minPt(unX, unY, 0), maxPt(unX, unY, 0);
	//			for (unsigned int i=0; i<regionPts.size(); i++) {
	//				if (minPt.x > regionPts[i].x) minPt.x = regionPts[i].x;
	//				if (minPt.y > regionPts[i].y) minPt.y = regionPts[i].y;
	//				if (maxPt.x < regionPts[i].x) maxPt.x = regionPts[i].x;
	//				if (maxPt.y < regionPts[i].y) maxPt.y = regionPts[i].y;
	//			}

	//			float tx, ty, tz, scale;
	//			trackball.getTranslate(tx, ty, tz);
	//			trackball.getScale(scale);

	//			minPt = Vector3f(-tx, -ty, 0) + minPt/scale;
	//			maxPt = Vector3f(-tx, -ty, 0) + maxPt/scale;

	//			bool valueToChangeTo = true;
	//			if (deselectRegion) valueToChangeTo = false;

	//			dx = (visManagers[dataSetToModify].getXPosition(0) - visManagers[dataSetToModify].getXPosition(xDim-1))/(float)((xDim-1) * 2.0);
	//			dy = (visManagers[dataSetToModify].getYPosition(0) - visManagers[dataSetToModify].getYPosition(yDim-1))/(float)((yDim-1) * 2.0);

	//			for (int i=0; i<xDim; i++) {
	//				float xPos1 = visManagers[dataSetToModify].getXPosition(i);
	//				float xPos2 = xPos1+dx;

	//				if ((xPos1 >= minPt.x && xPos1 <= maxPt.x) || (xPos2 >= minPt.x && xPos2 <= maxPt.x) ) {
	//					if (xPos1 > xPos2) {
	//						float temp = xPos1;
	//						xPos1 = xPos2;
	//						xPos2 = temp;
	//					}

	//					for (int j=0; j<yDim; j++) {
	//						float yPos1 = visManagers[dataSetToModify].getYPosition(j);
	//						float yPos2 = yPos1+dy;

	//						if ((yPos1 >= minPt.y && yPos1 <= maxPt.y) || (yPos2 >= minPt.y && yPos2 <= maxPt.y)) {
	//							if (yPos1 > yPos2) {
	//								float temp = yPos1;
	//								yPos1 = yPos2;
	//								yPos2 = temp;
	//							}

	//							bool isInRegion = false;

	//							// point in polygon test
	//							int npol = (int)regionPts.size();
	//							int idxI, idxJ, c = 0;

	//							float x = (xPos1+xPos2)*0.5;
	//							float y = (yPos1+yPos2)*0.5;

	//							for (idxI = 0, idxJ = npol-1; idxI < npol; idxJ = idxI++) {
	//								Vector3f pos1 = Vector3f(-tx, -ty, 0) + regionPts[idxI]/scale;
	//								Vector3f pos2 = Vector3f(-tx, -ty, 0) + regionPts[idxJ]/scale;

	//								float x1 = pos1.x, y1 = pos1.y;
	//								float x2 = pos2.x, y2 = pos2.y;

	//								if ((((y1 <= y) && (y < y2)) || ((y2 <= y) && (y < y1))) &&
	//									(x < (x2 - x1) * (y - y1) / (y2 - y1) + x1))
	//										c = !c;
	//							}

	//							if (c!=0) isInRegion = true;

	//							if (isInRegion) visManagers[dataSetToModify].data->volume->pointInMap[i][j] = valueToChangeTo;
	//						}
	//					}
	//				}
	//			}

#endif
