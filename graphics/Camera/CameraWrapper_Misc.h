#ifndef CAMERA_WRAPPER_MISC_H
#define CAMERA_WRAPPER_MISC_H

#include "CameraWrapper.h"

namespace CameraWrapper_GUI {
	void glutMouseCallback(CameraWrapper *camera, int button, int state, int x, int y) {
		if (state == GLUT_DOWN) {
			if (button == GLUT_LEFT_BUTTON) camera->recordLeftMouseDown(x, y);
			if (button == GLUT_RIGHT_BUTTON) camera->recordRightMouseDown(x, y);
			if (button == GLUT_MIDDLE_BUTTON) camera->recordMiddleMouseDown(x, y);
		}

		else if (state == GLUT_UP) {
			if (button == GLUT_LEFT_BUTTON) camera->recordLeftMouseUp(x, y);
			if (button == GLUT_RIGHT_BUTTON) camera->recordRightMouseUp(x, y);
			if (button == GLUT_MIDDLE_BUTTON) camera->recordMiddleMouseUp(x, y);
		}
	}

	void glutMouseMotionCallback(CameraWrapper *camera, int x, int y) {
		camera->recordMouseMotion(x, y);
	}

	void glutResizeCallback(CameraWrapper *camera, int width, int height) {
		camera->resizeViewport(width, height);
	}
}

#endif