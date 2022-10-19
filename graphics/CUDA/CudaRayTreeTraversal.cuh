#ifndef CUDA_RAY_TREE_TRAVERSAL_CUH
#define CUDA_RAY_TREE_TRAVERSAL_CUH

#include "CudaUtils.cu"
#include "CudaRay.cu"
#include "CudaMatrixTree.cuh"
#include "CudaRayCasterConstants.cuh"

#define CUDA_RAY_TREE_TRAVERSAL_EPS 0.0001

__device__ int mtGetRayBounds(const float3 rayOrigin, const float3 rayDirection, const MTNode &node, float &tMin, float &tMax) {
	float3 nodeMin = matricesGetCellMinPos(node.index, node.level, NODE_MATRIX_ID);
	float3 nodeMax = nodeMin + matricesCellSizes[NODE_MATRIX_ID][node.level];
	return intersectBox(rayOrigin, rayDirection, nodeMin, nodeMax, &tMin, &tMax);
}

__device__ int mtMoveRayOverChildNode(const float3 rayOrigin, const float3 rayDirection, float &rayT,
									  const float &rayTEnd, MTNode &leaf, const uint3 childIdx) {
	// the current node exists and is a leaf (for the given point, but not at bottom level of tree), so skip empty space of child
	leaf.index = childIdx;
	leaf.level++;
	
	float3 nodeMin = matricesGetCellMinPos(leaf.index, leaf.level, NODE_MATRIX_ID);
	float3 nodeMax = nodeMin + matricesCellSizes[NODE_MATRIX_ID][leaf.level];
	
	float tFar;
	int intersects = getEndOfBoxIntersection(rayOrigin, rayDirection, nodeMin, nodeMax, &tFar);
	//if (tFar < 0) return 0; // should not get here

	// move ray past child node
	rayT = tFar + CUDA_RAY_TREE_TRAVERSAL_EPS;

	// ray is now outside tree's bounding box
	if (rayTEnd-rayT < CUDA_RAY_TREE_TRAVERSAL_EPS) return 0;

	// update the node we are in at the current level
	float3 pos = rayOrigin + rayDirection*rayT;
	leaf.index = mtGetNodeIndex(pos, leaf.level);

	// old way, does not seem to make much difference either way, but above does not have branching so using that...
	/*uint3 numCells = matricesNumberOfCells[NODE_MATRIX_ID][leaf.level];

	if		(pos.x <= nodeMin.x) { if (leaf.index.x==0)			   return 0;  else leaf.index.x--; }
	else if (pos.x >= nodeMax.x) { if (leaf.index.x==numCells.x-1) return 0;  else leaf.index.x++; }
	if		(pos.y <= nodeMin.y) { if (leaf.index.y==0)			   return 0;  else leaf.index.y--; }
	else if (pos.y >= nodeMax.y) { if (leaf.index.y==numCells.y-1) return 0;  else leaf.index.y++; }
	if		(pos.z <= nodeMin.z) { if (leaf.index.z==0)		       return 0;  else leaf.index.z--; }
	else if (pos.z >= nodeMax.z) { if (leaf.index.z==numCells.z-1) return 0;  else leaf.index.z++; }*/


	return 1;
}

__device__ int mtGetNextLeaf(const float3 rayOrigin, const float3 rayDirection, float &rayT,
							 const float &rayTEnd, MTNode &leaf, const bool getNodesOnlyAtBottom) {
	#pragma unroll 8
	while (1) {
		uint3 childIdx;

		if (mtNodeExists(leaf)) {
			if (leaf.level >= matricesLevels[NODE_MATRIX_ID]-1) break; // the found node is the leaf we are looking for

			// keep going down the tree until the leaf is found, linear search seems to perform better
			float3 pos = rayOrigin + rayDirection*rayT;
			mtSearchDownTreeToFindLeaf_Linear(pos, leaf, childIdx);

			if (!getNodesOnlyAtBottom || leaf.level >= matricesLevels[NODE_MATRIX_ID]-1) break; // the found node is the leaf we are looking for
		}
		else {
			// go up the tree until the leaf is found, binary search seems to perform better
			mtBinarySearchAlongPath(0, leaf, leaf, childIdx);
			//mtSearchUpTreeToFindLeaf_Linear(leaf, childIdx);

			if (!getNodesOnlyAtBottom) break;
		}

		if (!mtMoveRayOverChildNode(rayOrigin, rayDirection, rayT, rayTEnd, leaf, childIdx)) return 0;
	}

	return 1;
}

__device__ int mtFindInitialLeaf(const float3 rayOrigin, const float3 rayDirection, float &rayT,
								 const float &rayTEnd, MTNode &leaf, const bool getNodesOnlyAtBottom) {
	float3 pos = rayOrigin + rayDirection*rayT;

	uint3 childIdx;
	int foundLeaf = mtFindLeaf(pos, leaf, childIdx);
	//if (!foundLeaf) return 0; // should not get here

	if (!getNodesOnlyAtBottom || leaf.level >= matricesLevels[NODE_MATRIX_ID]-1) return 1;

	if (!mtMoveRayOverChildNode(rayOrigin, rayDirection, rayT, rayTEnd, leaf, childIdx)) return 0;
	else return mtGetNextLeaf(rayOrigin, rayDirection, rayT, rayTEnd, leaf, getNodesOnlyAtBottom);
}

#endif