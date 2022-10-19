#pragma once
#ifndef CUDA_CAMERA_H
#define CUDA_CAMERA_H

#include "../common code/CUDA/CudaRay.cu"
#include "../common code/CUDA/CudaUtils.cu"

__constant__ float3 imagePlaneLowerLeft, imagePlaneLowerRight, imagePlaneTopLeft, imagePlaneTopRight, cameraPos;
__constant__ float cameraNearZ, cameraFarZ;

extern "C"
void setCameraVariables(float* nearCameraPlaneLowerLeft, float* nearCameraPlaneLowerRight, float* nearCameraPlaneTopLeft,
						float* nearCameraPlaneTopRight, float* cameraPosition, float zNear, float zFar) {
	cutilSafeCall(cudaMemcpyToSymbol(imagePlaneLowerLeft, nearCameraPlaneLowerLeft, 3*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(imagePlaneLowerRight, nearCameraPlaneLowerRight, 3*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(imagePlaneTopLeft, nearCameraPlaneTopLeft, 3*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(imagePlaneTopRight, nearCameraPlaneTopRight, 3*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(cameraPos, cameraPosition, 3*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(cameraNearZ, &zNear, sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(cameraFarZ, &zFar, sizeof(float)));
}

__constant__ unsigned int cameraImageWidth, cameraImageHeight;

// previously was also pre-recording invRayDirection, slope classification and slope values
//  but CUDA did not like that much because of all the malloc/free being used

cudaArray *cameraRayPositions = 0, *cameraRayDirections = 0, *cameraRayIntersectsVolume = 0;
texture<float4, 2, cudaReadModeElementType> cameraRayPositionsTex, cameraRayDirectionsTex;
texture<int, 2, cudaReadModeElementType> cameraRayIntersectsVolumeTex;

extern "C"
void freeComputedCameraRays() {
	if (cameraRayPositions != NULL) cutilSafeCall(cudaFreeArray(cameraRayPositions));
	if (cameraRayDirections != NULL) cutilSafeCall(cudaFreeArray(cameraRayDirections));
	if (cameraRayIntersectsVolume != NULL) cutilSafeCall(cudaFreeArray(cameraRayIntersectsVolume));
}

__device__ int computeRayFromCamera(float x, float y, const float3 bbMin, const float3 bbMax,
									float3 &rayOrigin, float3 &rayDirection, float &tNear, float &tFar) {
	float tx = x / (float)cameraImageWidth;
	float ty = y / (float)cameraImageHeight;

	float3 yLerp1 = lerp(imagePlaneLowerLeft, imagePlaneTopLeft, ty);
	float3 yLerp2 = lerp(imagePlaneLowerRight, imagePlaneTopRight, ty);

	rayOrigin = lerp(yLerp1, yLerp2, tx);

	rayDirection = normalize(rayOrigin-cameraPos);

	//float3 invRayDirection;
	//int slopeClassification;
	//float slopes[12];
	//computeInverseRayDirection(rayOrigin, rayDirection, &invRayDirection, &slopeClassification, slopes);

	const float eps = 0.001;

	// find intersection with data bounding box
	//int hit = intersectBox(rayOrigin, invRayDirection, slopeClassification, slopes, dataBoundingBoxMin, dataBoundingBoxMax, &tNear, &tFar);
	int hit = intersectBox(rayOrigin, rayDirection, bbMin, bbMax, &tNear, &tFar);
	if (!hit || tNear+2.0*eps > tFar) return 0;
	else if (tNear < 0.0f) tNear = 0.0f;  // clamp to near plane
	else tNear += eps; // little bump due to floating point errors

	return 1;
}

__global__ void computeRaysFromCameraKernel(const float3 dataBoundingBoxMin, const float3 dataBoundingBoxMax, float4 *rayPositions,
											float4 *rayDirections, int *rayIntersectsVolume) {
	unsigned int x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	unsigned int y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;

	if (x >= cameraImageWidth || y >= cameraImageHeight) return;

	unsigned int idx = x + y*cameraImageWidth;

	float3 rayOrigin, rayDirection;
	float tNear, tFar;
	if (computeRayFromCamera(x, y, dataBoundingBoxMin, dataBoundingBoxMax, rayOrigin, rayDirection, tNear, tFar)) {
		rayPositions[idx] = make_float4(rayOrigin, tNear);
		rayDirections[idx] = make_float4(rayDirection, tFar);
		rayIntersectsVolume[idx] = 1;
	}
	else rayIntersectsVolume[idx] = 0;
}

extern "C"
void cudaLoadCameraImageDimensions(const unsigned int imageW, const unsigned int imageH) {
	cutilSafeCall(cudaMemcpyToSymbol(cameraImageWidth, &imageW, sizeof(unsigned int)));
	cutilSafeCall(cudaMemcpyToSymbol(cameraImageHeight, &imageH, sizeof(unsigned int)));
}

extern "C"
void computeRaysFromCamera(const unsigned int imageW, const unsigned int imageH, const float3 dataBoundingBoxMin, const float3 dataBoundingBoxMax) {
	unsigned int lastWidth, lastHeight;
	cutilSafeCall(cudaMemcpyFromSymbol((void*)&lastWidth, cameraImageWidth, sizeof(unsigned int)));
	cutilSafeCall(cudaMemcpyFromSymbol((void*)&lastHeight, cameraImageHeight, sizeof(unsigned int)));

	cudaLoadCameraImageDimensions(imageW, imageH);

	dim3 blockSize = dim3(8, 8, 1);
	dim3 gridSize = dim3((imageW % blockSize.x != 0) ? (imageW / blockSize.x + 1) : (imageW / blockSize.x),
						 (imageH % blockSize.y != 0) ? (imageH / blockSize.y + 1) : (imageH / blockSize.y), 1);

	float4 *rayPositions, *rayDirections;
	int *rayIntersectsVolume;

	unsigned int numberOfRays = imageW*imageH;

	cudaMalloc((void**)&rayPositions, sizeof(float4)*numberOfRays);
	cudaMalloc((void**)&rayDirections, sizeof(float4)*numberOfRays);
	cudaMalloc((void**)&rayIntersectsVolume, sizeof(int)*numberOfRays);

	computeRaysFromCameraKernel<<<gridSize, blockSize>>>(dataBoundingBoxMin, dataBoundingBoxMax, rayPositions,
														 rayDirections, rayIntersectsVolume);

	if (lastWidth != imageW || lastHeight != imageH) {
		freeComputedCameraRays();

		loadAreaToTexture(rayPositions, imageW, imageH, &cameraRayPositions, &cameraRayPositionsTex, cudaMemcpyDeviceToDevice);
		loadAreaToTexture(rayDirections, imageW, imageH, &cameraRayDirections, &cameraRayDirectionsTex, cudaMemcpyDeviceToDevice);
		loadAreaToTexture(rayIntersectsVolume, imageW, imageH, &cameraRayIntersectsVolume, &cameraRayIntersectsVolumeTex, cudaMemcpyDeviceToDevice);
	}
	else {
		cutilSafeCall(cudaMemcpyToArray(cameraRayPositions, 0, 0, (void*)rayPositions, numberOfRays*sizeof(float4), cudaMemcpyDeviceToDevice));
		cutilSafeCall(cudaMemcpyToArray(cameraRayDirections, 0, 0, (void*)rayDirections, numberOfRays*sizeof(float4), cudaMemcpyDeviceToDevice));
		cutilSafeCall(cudaMemcpyToArray(cameraRayIntersectsVolume, 0, 0, (void*)rayIntersectsVolume, numberOfRays*sizeof(int), cudaMemcpyDeviceToDevice));
	}

	cutilSafeCall(cudaFree(rayPositions));  rayPositions = NULL;
	cutilSafeCall(cudaFree(rayDirections));  rayDirections = NULL;
	cutilSafeCall(cudaFree(rayIntersectsVolume));  rayIntersectsVolume = NULL;
}

#endif