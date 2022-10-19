#ifndef CUDA_CELL_3_D_H
#define CUDA_CELL_3_D_H

#include <cutil_math.h>

enum CudaCellTypes { CUDA_CELL_UNKNOWN, CUDA_CELL_TETRAHEDRON, CUDA_CELL_PRISM, CUDA_CELL_PYRAMID, CUDA_CELL_HEXAHEDRON };

__device__ float4 tetrahedronGetBarycentricCoordinates(const float3 pt, const float3 v0, const float3 v1, const float3 v2, const float3 v3) {
	float3 a1 = v0-v3, a2 = v1-v3, a3 = v2-v3;
	float3 B = pt-v3;
	float3 A = cross(a2,a3);

	float invDenominator = 1.0f / dot(A,a1);

	float4 bary;
	bary.x = dot(A,B) * invDenominator;
	bary.y = -dot(cross(a1,a3), B) * invDenominator;
	bary.z = dot(cross(a1,a2), B) * invDenominator;
	bary.w = 1.0f - bary.x - bary.y - bary.z;
	return bary;
}

__device__ int tetrahedronIsPointInside(const float3 pt, const float3 v0, const float3 v1, const float3 v2, const float3 v3) {
	float3 a1 = v0-v3, a2 = v1-v3, a3 = v2-v3;
	float3 B = pt-v3;
	float3 A = cross(a2,a3);

	float denominator = dot(A,a1), bary = dot(A,B);
	float alpha = denominator*1.0e-6, dMinusB = denominator-bary;

	if (denominator > 0 ? (bary>alpha && dMinusB>alpha) : (bary<alpha && dMinusB<alpha)) {
		bary = -dot(cross(a1,a3), B); // = bary.y * denominator
		dMinusB -= bary;

		if (denominator > 0 ? (bary>alpha && dMinusB>alpha) : (bary<alpha && dMinusB<alpha)) {
			bary = dot(cross(a1,a2), B); // = bary.z * denominator
			dMinusB -= bary;

			if (denominator > 0 ? (bary>alpha && dMinusB>alpha) : (bary<alpha && dMinusB<alpha)) return 1;
		}
	}
	return 0;
}

__device__ float tetrahedronInterpolateValue(const float3 pt, const float3 v0, const float3 v1, const float3 v2, const float3 v3,
											 const float val0, const float val1, const float val2, const float val3) {
	float4 bary = tetrahedronGetBarycentricCoordinates(pt, v0, v1, v2, v3);
	return val0*bary.x + val1*bary.y + val2*bary.z + val3*bary.w;
}

#endif