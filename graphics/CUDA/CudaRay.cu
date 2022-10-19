#ifndef CUDA_RAY_H
#define CUDA_RAY_H

#include <cutil_inline.h>
#include <cutil_math.h>

#ifndef GPU_RAY_CLASSIFICATION
#define GPU_RAY_CLASSIFICATION
enum GPURayClassification
{	GPU_MMM = 0, GPU_MMP = 1, GPU_MMO = 2, GPU_MPM = 4, GPU_MPP = 5, GPU_MPO = 6, GPU_MOM = 8, GPU_MOP = 9, GPU_MOO = 10,
	GPU_PMM = 16, GPU_PMP = 17, GPU_PMO = 18, GPU_PPM = 20, GPU_PPP = 21, GPU_PPO = 22, GPU_POM = 24, GPU_POP = 25, GPU_POO = 26,
	GPU_OMM = 32, GPU_OMP = 33, GPU_OMO = 34, GPU_OPM = 36, GPU_OPP = 37, GPU_OPO = 38, GPU_OOM = 40, GPU_OOP = 41 };
#endif

// "Fast Ray/Axis-Aligned Bounding Box Overlap Tests using Ray Slopes"
__device__ void computeInverseRayDirection(float3 rayOrigin, float3 rayDirection, float3 *invRayDirection, int *raySlopeClassification, float *raySlopes) {
	invRayDirection->x = 1.0f/rayDirection.x;
	invRayDirection->y = 1.0f/rayDirection.y;
	invRayDirection->z = 1.0f/rayDirection.z;

	//ray slope
	raySlopes[0] = rayDirection.x * invRayDirection->y;
	raySlopes[1] = rayDirection.y * invRayDirection->x;
	raySlopes[2] = rayDirection.y * invRayDirection->z;
	raySlopes[3] = rayDirection.z * invRayDirection->y;
	raySlopes[4] = rayDirection.x * invRayDirection->z;
	raySlopes[5] = rayDirection.z * invRayDirection->x;

	raySlopes[6]  = rayOrigin.y - raySlopes[1] * rayOrigin.x;
	raySlopes[7]  = rayOrigin.z - raySlopes[5] * rayOrigin.x;
	raySlopes[8]  = rayOrigin.x - raySlopes[0] * rayOrigin.y;
	raySlopes[9]  = rayOrigin.z - raySlopes[3] * rayOrigin.y;
	raySlopes[10] = rayOrigin.x - raySlopes[4] * rayOrigin.z;
	raySlopes[11] = rayOrigin.y - raySlopes[2] * rayOrigin.z;	

	//ray slope classification
	*raySlopeClassification = 0;
	if (rayDirection.z > 0) *raySlopeClassification += 1;
	else if (rayDirection.z == 0) *raySlopeClassification += 2;
	if (rayDirection.y > 0) *raySlopeClassification += 4;
	else if (rayDirection.y == 0) *raySlopeClassification += 8;
	if (rayDirection.x > 0) *raySlopeClassification += 16;
	else if (rayDirection.x == 0) *raySlopeClassification += 32;
}

__device__ int getEndOfBoxIntersection(const float3 rayOrigin, const float3 invRayDirection, const int raySlopeClassification,
									   const float *raySlopes, const float3 boxMin, float3 boxMax, float *tFar) {
	switch (raySlopeClassification) {
		case GPU_MMM: {
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_MMP: {		
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_MPM: {		
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) 
				|| (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_MPP: {
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) 
				|| (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0) 
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PMM: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PMP: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PPM: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0)
				|| (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PPP: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0)
				|| (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_OMM: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)) return false;

			*tFar = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_OMP: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)) return false;

			*tFar = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_OPM: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)) return false;

			*tFar = (boxMax.y - rayOrigin.y) * invRayDirection.y;		
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_OPP: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0)) return false;
			
			*tFar = (boxMax.y - rayOrigin.y) * invRayDirection.y;		
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MOM: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.z < boxMin.z) 
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MOP: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.z > boxMax.z) 
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_POM: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_POP: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x > boxMax.x) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MMO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y)  
				|| (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0) || (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0)) return false;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MPO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) 
				|| (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) || (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0)) return false;
			
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_PMO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) 
				|| (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0) || (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0)) return false;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_PPO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) 
				|| (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0) || (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0)) return false;
		
			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MOO: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			return true;
		}
		case GPU_POO: {
			if((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			return true;
		}
		case GPU_OMO: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;
			
			*tFar = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			return true;
		}
		case GPU_OPO: {
			if((rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tFar = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			return true;
		}
		case GPU_OOM: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y))
				return false;

			*tFar = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			return true;
		}
		case GPU_OOP: {
			if((rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y))
				return false;

			*tFar = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			return true;
		}
	}
	return false;
}

__device__ int intersectBox(const float3 rayOrigin, const float3 invRayDirection, const int raySlopeClassification,
							const float *raySlopes, const float3 boxMin, const float3 boxMax, float *tNear) {
	switch (raySlopeClassification) {
		case GPU_MMM: {
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_MMP: {		
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_MPM: {		
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) 
				|| (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_MPP: {
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) 
				|| (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0) 
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_PMM: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_PMP: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_PPM: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0)
				|| (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_PPP: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0)
				|| (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;
			return true;
		}
		case GPU_OMM: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)) return false;

			*tNear = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_OMP: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)) return false;

			*tNear = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_OPM: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)) return false;

			*tNear = (boxMin.y - rayOrigin.y) * invRayDirection.y;		
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_OPP: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0)) return false;
			
			*tNear = (boxMin.y - rayOrigin.y) * invRayDirection.y;		
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_MOM: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.z < boxMin.z) 
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_MOP: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.z > boxMax.z) 
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;

			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_POM: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_POP: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x > boxMax.x) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_MMO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y)  
				|| (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0) || (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0)) return false;

			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_MPO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) 
				|| (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) || (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_PMO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) 
				|| (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0) || (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0)) return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_PPO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) 
				|| (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0) || (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0)) return false;
		
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			return true;
		}
		case GPU_MOO: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			return true;
		}
		case GPU_POO: {
			if((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			return true;
		}
		case GPU_OMO: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;
			
			*tNear = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			return true;
		}
		case GPU_OPO: {
			if((rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tNear = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			return true;
		}
		case GPU_OOM: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y))
				return false;

			*tNear = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			return true;
		}
		case GPU_OOP: {
			if((rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y))
				return false;

			*tNear = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			return true;
		}
	}
	return false;
}

__device__ int intersectBox(const float3 rayOrigin, const float3 invRayDirection, const int raySlopeClassification,
							const float *raySlopes, const float3 boxMin, const float3 boxMax, float *tNear, float *tFar) {
	switch (raySlopeClassification) {
		case GPU_MMM: {
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_MMP: {		
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_MPM: {		
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) 
				|| (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_MPP: {
			if ((rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) 
				|| (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0) || (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0) 
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PMM: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PMP: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0)
				|| (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PPM: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0)
				|| (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_PPP: {
			if ((rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z) || (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0)
				|| (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0) || (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;
			float t2 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t2 > *tNear) *tNear = t2;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			t2 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t2 < *tFar) *tFar = t2;
			return true;
		}
		case GPU_OMM: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[2] * boxMin.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMax.z + raySlopes[9] > 0)) return false;

			*tNear = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_OMP: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[2] * boxMax.z - boxMax.y + raySlopes[11] > 0) || (raySlopes[3] * boxMin.y - boxMin.z + raySlopes[9] < 0)) return false;

			*tNear = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_OPM: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[2] * boxMin.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMax.z + raySlopes[9] > 0)) return false;

			*tNear = (boxMin.y - rayOrigin.y) * invRayDirection.y;		
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMax.y - rayOrigin.y) * invRayDirection.y;		
			t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_OPP: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[2] * boxMax.z - boxMin.y + raySlopes[11] < 0) || (raySlopes[3] * boxMax.y - boxMin.z + raySlopes[9] < 0)) return false;
			
			*tNear = (boxMin.y - rayOrigin.y) * invRayDirection.y;		
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMax.y - rayOrigin.y) * invRayDirection.y;		
			t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MOM: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.z < boxMin.z) 
				|| (raySlopes[5] * boxMin.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMax.x + raySlopes[10] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MOP: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.z > boxMax.z) 
				|| (raySlopes[5] * boxMin.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMax.x + raySlopes[10] > 0)) return false;

			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_POM: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z)
				|| (raySlopes[5] * boxMax.x - boxMax.z + raySlopes[7] > 0) || (raySlopes[4] * boxMin.z - boxMin.x + raySlopes[10] < 0)) return false;
			
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_POP: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.x > boxMax.x) || (rayOrigin.z > boxMax.z)
				|| (raySlopes[5] * boxMax.x - boxMin.z + raySlopes[7] < 0) || (raySlopes[4] * boxMax.z - boxMin.x + raySlopes[10] < 0)) return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MMO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y)  
				|| (raySlopes[1] * boxMin.x - boxMax.y + raySlopes[6] > 0) || (raySlopes[0] * boxMin.y - boxMax.x + raySlopes[8] > 0)) return false;

			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MPO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.y > boxMax.y) 
				|| (raySlopes[1] * boxMin.x - boxMin.y + raySlopes[6] < 0) || (raySlopes[0] * boxMax.y - boxMax.x + raySlopes[8] > 0)) return false;
			
			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_PMO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) 
				|| (raySlopes[1] * boxMax.x - boxMax.y + raySlopes[6] > 0) || (raySlopes[0] * boxMin.y - boxMin.x + raySlopes[8] < 0)) return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_PPO: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z) || (rayOrigin.x > boxMax.x) || (rayOrigin.y > boxMax.y) 
				|| (raySlopes[1] * boxMax.x - boxMin.y + raySlopes[6] < 0) || (raySlopes[0] * boxMax.y - boxMin.x + raySlopes[8] < 0)) return false;
		
			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			float t1 = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			if (t1 > *tNear) *tNear = t1;

			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			t1 = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			if (t1 < *tFar) *tFar = t1;
			return true;
		}
		case GPU_MOO: {
			if((rayOrigin.x < boxMin.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tNear = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			*tFar = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			return true;
		}
		case GPU_POO: {
			if((rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tNear = (boxMin.x - rayOrigin.x) * invRayDirection.x;
			*tFar = (boxMax.x - rayOrigin.x) * invRayDirection.x;
			return true;
		}
		case GPU_OMO: {
			if((rayOrigin.y < boxMin.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;
			
			*tNear = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			*tFar = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			return true;
		}
		case GPU_OPO: {
			if((rayOrigin.y > boxMax.y) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.z < boxMin.z) || (rayOrigin.z > boxMax.z))
				return false;

			*tNear = (boxMin.y - rayOrigin.y) * invRayDirection.y;
			*tFar = (boxMax.y - rayOrigin.y) * invRayDirection.y;
			return true;
		}
		case GPU_OOM: {
			if((rayOrigin.z < boxMin.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y))
				return false;

			*tNear = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			*tFar = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			return true;
		}
		case GPU_OOP: {
			if((rayOrigin.z > boxMax.z) || (rayOrigin.x < boxMin.x) || (rayOrigin.x > boxMax.x) || (rayOrigin.y < boxMin.y) || (rayOrigin.y > boxMax.y))
				return false;

			*tNear = (boxMin.z - rayOrigin.z) * invRayDirection.z;
			*tFar = (boxMax.z - rayOrigin.z) * invRayDirection.z;
			return true;
		}
	}
	return false;
}

// "An Efficient and Robust Ray–Box Intersection Algorithm"
//   the only thing different from "Efficiency issues for ray tracing" Brian Smits is that you test the invDir at each step
__device__ int intersectBox(const float3 rayOrigin, const float3 rayDirection, const float3 boxMin, const float3 boxMax,
							float *tNear, float *tFar) {
	float tMin, tMax;
	float invRayDirection = 1.0f / rayDirection.x;
	if (invRayDirection >= 0) {
		tMin = (boxMin.x - rayOrigin.x) * invRayDirection;
		tMax = (boxMax.x - rayOrigin.x) * invRayDirection;
	}
	else {
		tMin = (boxMax.x - rayOrigin.x) * invRayDirection;
		tMax = (boxMin.x - rayOrigin.x) * invRayDirection;
	}

	float tMin2, tMax2;
	invRayDirection = 1.0f / rayDirection.y;
	if (invRayDirection >= 0) {
		tMin2 = (boxMin.y - rayOrigin.y) * invRayDirection;
		tMax2 = (boxMax.y - rayOrigin.y) * invRayDirection;
	}
	else {
		tMin2 = (boxMax.y - rayOrigin.y) * invRayDirection;
		tMax2 = (boxMin.y - rayOrigin.y) * invRayDirection;
	}

	if (tMin > tMax2 || tMin2 > tMax) return false;
	if (tMin2 > tMin) tMin = tMin2;
	if (tMax2 < tMax) tMax = tMax2;

	invRayDirection = 1.0f / rayDirection.z;
	if (invRayDirection >= 0) {
		tMin2 = (boxMin.z - rayOrigin.z) * invRayDirection;
		tMax2 = (boxMax.z - rayOrigin.z) * invRayDirection;
	}
	else {
		tMin2 = (boxMax.z - rayOrigin.z) * invRayDirection;
		tMax2 = (boxMin.z - rayOrigin.z) * invRayDirection;
	}

	if (tMin > tMax2 || tMin2 > tMax) return false;
	if (tMin2 > tMin) tMin = tMin2;
	if (tMax2 < tMax) tMax = tMax2;

	*tNear = tMin;
	*tFar = tMax;

	return tMax > tMin;
}

__device__ int getEndOfBoxIntersection(const float3 rayOrigin, const float3 rayDirection, const float3 boxMin, const float3 boxMax,
									   float *tFar) {
	float tMin, tMax;
	float invRayDirection = 1.0f / rayDirection.x;
	if (invRayDirection >= 0) {
		tMin = (boxMin.x - rayOrigin.x) * invRayDirection;
		tMax = (boxMax.x - rayOrigin.x) * invRayDirection;
	}
	else {
		tMin = (boxMax.x - rayOrigin.x) * invRayDirection;
		tMax = (boxMin.x - rayOrigin.x) * invRayDirection;
	}

	float tMin2, tMax2;
	invRayDirection = 1.0f / rayDirection.y;
	if (invRayDirection >= 0) {
		tMin2 = (boxMin.y - rayOrigin.y) * invRayDirection;
		tMax2 = (boxMax.y - rayOrigin.y) * invRayDirection;
	}
	else {
		tMin2 = (boxMax.y - rayOrigin.y) * invRayDirection;
		tMax2 = (boxMin.y - rayOrigin.y) * invRayDirection;
	}

	if (tMin > tMax2 || tMin2 > tMax) return false;
	if (tMin2 > tMin) tMin = tMin2;
	if (tMax2 < tMax) tMax = tMax2;

	invRayDirection = 1.0f / rayDirection.z;
	if (invRayDirection >= 0) {
		tMin2 = (boxMin.z - rayOrigin.z) * invRayDirection;
		tMax2 = (boxMax.z - rayOrigin.z) * invRayDirection;
	}
	else {
		tMin2 = (boxMax.z - rayOrigin.z) * invRayDirection;
		tMax2 = (boxMin.z - rayOrigin.z) * invRayDirection;
	}

	if (tMin > tMax2 || tMin2 > tMax) return false;
	if (tMin2 > tMin) tMin = tMin2;
	if (tMax2 < tMax) tMax = tMax2;

	*tFar = tMax;

	return tMax > tMin;
}

__device__ int rayIntersectsSphere(const float3 rayOrigin, const float3 rayDirection, const float3 sphereCenter,
								   const float sphereRadiusSqr, float &tNear, float &tFar) {
	float3 diff = rayOrigin - sphereCenter;
	float A = dot(diff, diff) - sphereRadiusSqr;
	float B = dot(diff, rayDirection);
	if (A <= 0) { // inside sphere
		float C = B*B - A;
		tNear = 0;
		tFar = -B + sqrt(C);
		return 1;
	}
	else if (B >= 0) return 0;
	
	float C = B*B - A;
	if (C < 0) return 0;
	else {
		float root = sqrt(C);
		tNear = -B - root;
		tFar = -B + root;
		return 1;
	}
}

#endif