#ifndef CUDA_MATRIX_SPATIAL_STRUCTURE_CUH
#define CUDA_MATRIX_SPATIAL_STRUCTURE_CUH

#include "CudaUtils.cu"
#include "CudaMatrices.cuh"
#include "CudaConstants.cuh"



struct MTNode {
	unsigned int level;
	uint3 index;
};

__device__ void mtPrintNode(const MTNode &node) { cudaPrint("%d %d %d %d\n", node.index.x, node.index.y, node.index.z, node.level); }

__device__ int mtNodeIsEqual(const MTNode &node1, const MTNode &node2)
	{ return node1.index.x == node2.index.x && node1.index.y == node2.index.y && node1.index.z == node2.index.z && node1.level == node2.level; }
__device__ int mtNodeIsEqual(const MTNode &node, const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int level) {
	MTNode n;  n.index.x = x;  n.index.y = y;  n.index.z = z;  n.level = level;
	return mtNodeIsEqual(node, n);
}



__device__ void mtGetCornerIndexAtLeafLevel(MTNode &node, const int cornerID) {
	uint3 factor = (matricesNumberOfCells[CORNER_MATRIX_ID][matricesLevels[CORNER_MATRIX_ID]-1]-make_uint3(1)) / (matricesNumberOfCells[NODE_MATRIX_ID][node.level]);

	node.level = matricesLevels[CORNER_MATRIX_ID]-1;
	switch (cornerID) {
		case CORNER_LBB:			node.index.x = __umul24(factor.x, node.index.x);   node.index.y = __umul24(factor.y, node.index.y);   node.index.z = __umul24(factor.z, node.index.z);   break;
		case CORNER_LBF:			node.index.x = __umul24(factor.x, node.index.x);   node.index.y = __umul24(factor.y, node.index.y);   node.index.z = __umul24(factor.z, node.index.z+1); break;
		case CORNER_LTB:			node.index.x = __umul24(factor.x, node.index.x);   node.index.y = __umul24(factor.y, node.index.y+1); node.index.z = __umul24(factor.z, node.index.z);   break;
		case CORNER_LTF:			node.index.x = __umul24(factor.x, node.index.x);   node.index.y = __umul24(factor.y, node.index.y+1); node.index.z = __umul24(factor.z, node.index.z+1); break;
		case CORNER_RBB:			node.index.x = __umul24(factor.x, node.index.x+1); node.index.y = __umul24(factor.y, node.index.y);   node.index.z = __umul24(factor.z, node.index.z);   break;
		case CORNER_RBF:			node.index.x = __umul24(factor.x, node.index.x+1); node.index.y = __umul24(factor.y, node.index.y);   node.index.z = __umul24(factor.z, node.index.z+1); break;
		case CORNER_RTB:			node.index.x = __umul24(factor.x, node.index.x+1); node.index.y = __umul24(factor.y, node.index.y+1); node.index.z = __umul24(factor.z, node.index.z);   break;
		default:/*CORNER_RTF*/	node.index.x = __umul24(factor.x, node.index.x+1); node.index.y = __umul24(factor.y, node.index.y+1); node.index.z = __umul24(factor.z, node.index.z+1); break;
	}
}

#ifdef CUDA_MSS_COMPILE_AS_OCTREE

	__device__ uint3 mtGetParentIndices(const MTNode child) {
		return make_uint3(child.index.x >> 1, child.index.y >> 1, child.index.z >> 1);
	}

	__device__ uint3 mtGetParentIndices(const unsigned int levelDifference, const MTNode child) {
		return make_uint3(child.index.x >> levelDifference, child.index.y >> levelDifference, child.index.z >> levelDifference);
	}

	/*__device__ uint3 mtGetParentIndices(const unsigned int parentLevel, const MTNode child) {
		unsigned int factor = mtNodeNumberOfCells[child.level] / mtNodeNumberOfCells[parentLevel];
		return make_uint3(child.index.x / factor, child.index.y / factor, child.index.z / factor);
	}*/

	__device__ uint3 mtGetChildIndices(const MTNode parent, const unsigned int childID) {
		switch (childID) {
			case CORNER_LBB:		return make_uint3( parent.index.x<<1,     parent.index.y<<1,     parent.index.z<<1);
			case CORNER_LBF:		return make_uint3( parent.index.x<<1,     parent.index.y<<1,    (parent.index.z<<1)+1);
			case CORNER_LTB:		return make_uint3( parent.index.x<<1,    (parent.index.y<<1)+1,  parent.index.z<<1);
			case CORNER_LTF:		return make_uint3( parent.index.x<<1,    (parent.index.y<<1)+1, (parent.index.z<<1)+1);
			case CORNER_RBB:		return make_uint3((parent.index.x<<1)+1,  parent.index.y<<1,     parent.index.z<<1);
			case CORNER_RBF:		return make_uint3((parent.index.x<<1)+1,  parent.index.y<<1,    (parent.index.z<<1)+1);
			case CORNER_RTB:		return make_uint3((parent.index.x<<1)+1, (parent.index.y<<1)+1,  parent.index.z<<1);
			default:/*CORNER_RTF*/	return make_uint3((parent.index.x<<1)+1, (parent.index.y<<1)+1, (parent.index.z<<1)+1);
		}
	}

	__device__ int mtGetChildThatPosLiesIn(const float3 pos, const MTNode &parent) {
		float3 nodeCenter = matricesGetCellCenter(parent.index, parent.level, NODE_MATRIX_ID);
		if (pos.x <= nodeCenter.x) {
			if (pos.y <= nodeCenter.y) {
				if (pos.z <= nodeCenter.z) return CORNER_LBB;
				else return CORNER_LBF;
			}
			else {
				if (pos.z <= nodeCenter.z) return CORNER_LTB;
				else return CORNER_LTF;
			}
		}
		else {
			if (pos.y <= nodeCenter.y) {
				if (pos.z <= nodeCenter.z) return CORNER_RBB;
				else return CORNER_RBF;
			}
			else {
				if (pos.z <= nodeCenter.z) return CORNER_RTB;
				else return CORNER_RTF;
			}
		}
	}

	__device__ void mtGetCornerIndex(MTNode &node, const int cornerID) {
		// "base" case that will result in an infinite loop if we don't check for it
		if (node.index.x == 0 && node.index.y == 0 && node.index.z == 0 && cornerID == CORNER_LBB) {
			node.level = 0;
			return;
		}

		// get index of the corner at lowest level
		mtGetCornerIndexAtLeafLevel(node, cornerID);

		// check to see if the corner data was flattened to a single matrix
		if (matricesLevels[CORNER_MATRIX_ID] == 1) return;

		// if all the indices are even, that means we need to go a level higher to find the correct place
		//while (node.index.x%2==0 && node.index.y%2==0 && node.index.z%2==0) {  // below uses the more efficient mod
		#pragma unroll 16
		while ((node.index.x&1)==0 && (node.index.y&1)==0 && (node.index.z&1)==0) {
			node.level--;
			node.index.x = node.index.x>>1;
			node.index.y = node.index.y>>1;
			node.index.z = node.index.z>>1;
		}
	}

#else // KD-Tree
	__constant__ int mtKDTreeSplitDirections[CUDA_MAX_MATRIX_LEVEL-1];
	__constant__ uint3 mtKDTreeSplitSums[CUDA_MAX_MATRIX_LEVEL];

	extern "C"
	void mtLoadKDTreeSplitDirections(unsigned int numSplits, int *splitDirections) {
		if (numSplits > CUDA_MAX_MATRIX_LEVEL-1) {
			printf("Height too large, increase constant variable CUDA_MAX_MATRIX_LEVEL\n");
			return;
		}
		
		uint3 *splitSums;  cutilSafeCall(cudaMallocHost((void**)&splitSums, (numSplits+1)*sizeof(uint3)));
		splitSums[0] = make_uint3(0,0,0);
		for (unsigned int i=0; i<numSplits; i++) {
			splitSums[i+1] = splitSums[i];
			if (splitDirections[i] == DIRECTION_X) splitSums[i+1].x++;
			else if (splitDirections[i] == DIRECTION_Y) splitSums[i+1].y++;
			else if (splitDirections[i] == DIRECTION_Z) splitSums[i+1].z++;
			else printf("Unknown split direction in mtLoadKDTreeSplitDirections()\n");
		}

		cutilSafeCall(cudaMemcpyToSymbol(mtKDTreeSplitDirections, splitDirections, numSplits*sizeof(int)));
		cutilSafeCall(cudaMemcpyToSymbol(mtKDTreeSplitSums, splitSums, (numSplits+1)*sizeof(uint3)));

		cutilSafeCall(cudaFreeHost(splitSums));  splitSums = NULL;
	}

	__device__ int mtGetSplitDirection(const unsigned int level) { return mtKDTreeSplitDirections[level]; }	


	__device__ uint3 mtGetParentIndices(const MTNode child) {
		switch (mtGetSplitDirection(child.level-1)) {
			case DIRECTION_X: { return make_uint3(child.index.x >> 1, child.index.y     , child.index.z); }
			case DIRECTION_Y: { return make_uint3(child.index.x     , child.index.y >> 1, child.index.z); }
			default: /*Z*/      { return make_uint3(child.index.x     , child.index.y     , child.index.z >> 1); }
		}
	}

	__device__ uint3 mtGetParentIndices(const unsigned int levelDifference, const MTNode child) {
		uint3 diff = mtKDTreeSplitSums[child.level] - mtKDTreeSplitSums[child.level-levelDifference];
		return make_uint3(child.index.x >> diff.x, child.index.y >> diff.y, child.index.z >> diff.z);
	}

	__device__ uint3 mtGetChildIndices(const MTNode parent, const unsigned int childID) {
		uint3 child = parent.index;
		switch (mtGetSplitDirection(parent.level)) {
			case DIRECTION_X: { child.x <<= 1;  if (childID) child.x++;  break; }
			case DIRECTION_Y: { child.y <<= 1;  if (childID) child.y++;  break; }
			default: /*Z*/      { child.z <<= 1;  if (childID) child.z++;  break; }
		}
		return child;
	}

	__device__ int mtGetChildThatPosLiesIn(const float3 pos, const MTNode &parent) {
		float3 nodeCenter = mtGetNodeCenter(parent);
		switch (mtGetSplitDirection(parent.level)) {
			case DIRECTION_X: { 
				if (pos.x <= nodeCenter.x) return DIRECTION_NEGATIVE;
				else return DIRECTION_POSITIVE;
			}
			case DIRECTION_Y: { 
				if (pos.y <= nodeCenter.y) return DIRECTION_NEGATIVE;
				else return DIRECTION_POSITIVE;
			}
			default: /*Z*/ {
				if (pos.z <= nodeCenter.z) return DIRECTION_NEGATIVE;
				else return DIRECTION_POSITIVE;
			}
		}
	}

	__device__ void mtGetCornerIndex(MTNode &node, const int cornerID) {
		// "base" case that will result in an infinite loop if we don't check for it
		if (node.index.x == 0 && node.index.y == 0 && node.index.z == 0 && cornerID == CORNER_LBB) {
			node.level = 0;
			return;
		}

		// get index of the corner at lowest level
		mtGetCornerIndexAtLeafLevel(node, cornerID);

		// check to see if the corner data was flattened to a single matrix
		if (mtCornerLevels == 1) return;

		// if all the indices are even, that means we need to go a level higher to find the correct place
		#pragma unroll 16
		while (node.level > 0) {
			switch (mtGetSplitDirection(node.level-1)) {
				case DIRECTION_X: {
					if ((node.index.x&1)==0) return;
					else { node.index.x >>= 1;  break; }
				}
				case DIRECTION_Y: {
					if ((node.index.y&1)==0) return;
					else { node.index.y >>= 1;  break; }
				}
				default: /*Z*/ {
					if ((node.index.z&1)==0) return;
					else { node.index.z >>= 1;  break; }
				}
			}

			node.level--;
		}
	}

#endif


__device__ float2 mtGetMinMax(const unsigned int dataIndex) { return cudaTexFetch(matricesGetFloat2DataTexture(NODE_MATRIX_ID), dataIndex); }
__device__ float2 mtGetMinMax(const MTNode &node) { return mtGetMinMax(blockGetDataIndex(node.index, node.level, NODE_MATRIX_ID)); }



//__device__ unsigned int mtGetNodeData(const MTNode &node) {
//	int type = matricesMatrixTypes[NODE_MATRIX_ID][node.level];
//	if (type == 0) { // GRID_REGULAR
//		unsigned int idx = fullMatrixGetIndex(node.index, matricesNumberOfCells[NODE_MATRIX_ID][node.level]);
//		return cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), idx + matricesDataOffsets[NODE_MATRIX_ID][node.level]);
//	}
//	else { // GRID_REGULAR_CSR
//		int relativeDataIndex = csrGetIndex(node.index, node.level, NODE_MATRIX_ID);
//		int dataOffset = matricesDataOffsets[NODE_MATRIX_ID][node.level];
//		return cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), dataOffset+relativeDataIndex);
//	}
//}
//

__device__ uint2 mtGetNodeDataAndNextInTexture(const MTNode &node) {
	//int type = matricesMatrixTypes[NODE_MATRIX_ID][node.level];
	//uint2 result;
	//if (type == 0) { // GRID_REGULAR
	//	unsigned int idx = fullMatrixGetIndex(node.index, matricesNumberOfCells[NODE_MATRIX_ID][node.level]);
	//	cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), idx + matricesDataOffsets[NODE_MATRIX_ID][node.level], result.x, result.y);
	//}
	//else { // GRID_REGULAR_CSR
	//	int relativeDataIndex = csrGetIndex(node.index, node.level, NODE_MATRIX_ID);
	//	int dataOffset = matricesDataOffsets[NODE_MATRIX_ID][node.level];
	//	cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), dataOffset+relativeDataIndex, result.x, result.y);
	//}
	//return result;

	unsigned int idx = blockGetDataIndex(node.index, node.level, NODE_MATRIX_ID);
	uint2 result;
	cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), idx, result.x, result.y);
	return result;
}

//
//__device__ int mtNodeExists(const MTNode &node, uint2 &vectorOffsets) {
//	int type = matricesMatrixTypes[NODE_MATRIX_ID][node.level];
//	if (type == 0) { // type = 0: GRID_REGULAR
//		unsigned int idx = fullMatrixGetIndex(node.index, matricesNumberOfCells[NODE_MATRIX_ID][node.level]);
//		unsigned int val = cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), idx + matricesDataOffsets[NODE_MATRIX_ID][node.level]);
//		if (val != matricesNullValue_uint[NODE_MATRIX_ID]) {
//			cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), idx + matricesDataOffsets[NODE_MATRIX_ID][node.level], vectorOffsets.x, vectorOffsets.y);
//			return 1;
//		}
//		else return 0;
//	}
//	else {
//		// when calling this function, assume grid is of type GRID_REGULAR_CSR and not GRID_REGULAR_SPARSE_NODATA
//		int idx = csrGetIndex(node.index, node.level, NODE_MATRIX_ID);
//
//		if (idx >= 0) {
//			cudaTexFetch(matricesGetUintDataTexture(NODE_MATRIX_ID), idx + matricesDataOffsets[NODE_MATRIX_ID][node.level], vectorOffsets.x, vectorOffsets.y);
//			return 1;
//		}
//		else return 0;
//	}
//}

__device__ int mtNodeExists(const MTNode &node) {
	return blockCellExists(node.index, node.level, NODE_EXISTS_MATRIX_ID, matricesGetCharDataTexture(NODE_EXISTS_MATRIX_ID), (char)0);
}

template <class T> __device__ int mtNodeExists(const MTNode &node, texture<T, 2, cudaReadModeElementType> tex, const T nullVal) {
	return blockCellExists(node.index, node.level, NODE_EXISTS_MATRIX_ID, tex, nullVal);
}

template <class T> __device__ void mtBinarySearchAlongPath(const unsigned int topLevel, const MTNode &bottomNode, MTNode &leafNode, texture<T, 2, cudaReadModeElementType> tex, const T nullVal) {
	unsigned int minLevel = topLevel, maxLevel = bottomNode.level;

	MTNode parentNode;
	//#pragma unroll 4 // only need to unroll lg(MAX_OCTREE_SIZE)
	while(maxLevel-minLevel > 1) {
		parentNode.level = (maxLevel+minLevel)>>1;

		parentNode.index = mtGetParentIndices(bottomNode.level-parentNode.level, bottomNode);

		if (mtNodeExists(parentNode, tex, nullVal)) minLevel = parentNode.level;
		else maxLevel = parentNode.level;
	}

	leafNode.index = mtGetParentIndices(bottomNode.level-minLevel, bottomNode);
	leafNode.level = minLevel;
}

//__device__ void mtBinarySearchAlongPath(const unsigned int topLevel, const MTNode &bottomNode,
//										 MTNode &leafNode, uint2 &vectorOffsets) {
//	unsigned int minLevel = topLevel, maxLevel = bottomNode.level;
//
//	MTNode parentNode;
//	//#pragma unroll 4 // only need to unroll lg(MAX_OCTREE_SIZE)
//	while(maxLevel-minLevel > 1) {
//		parentNode.level = (maxLevel+minLevel)>>1;
//
//		parentNode.index = mtGetParentIndices(bottomNode.level-parentNode.level, bottomNode);
//
//		uint2 offsets;
//		if (mtNodeExists(parentNode, offsets)) { minLevel = parentNode.level;  vectorOffsets = offsets; }
//		else maxLevel = parentNode.level;
//	}
//
//	leafNode.index = mtGetParentIndices(bottomNode.level-minLevel, bottomNode);
//	leafNode.level = minLevel;
//
//	if (topLevel-bottomNode.level <= 1) vectorOffsets = mtGetNodeDataAndNextInTexture(leafNode);
//}

__device__ void mtBinarySearchAlongPath(const unsigned int topLevel, const MTNode &bottomNode,
										 MTNode &leafNode, uint3 &childIdx) {
	unsigned int minLevel = topLevel, maxLevel = bottomNode.level;

	MTNode parentNode;
	while(maxLevel-minLevel > 1) {
		parentNode.level = (maxLevel+minLevel)>>1;

		parentNode.index = mtGetParentIndices(bottomNode.level-parentNode.level, bottomNode);

		//printf("b: "); mtPrintNode(parentNode);
		//if (mtNodeExists(parentNode)) printf("\t\ttrue\n");
		//else printf("\t\tfalse\n");

		if (mtNodeExists(parentNode)) minLevel = parentNode.level;
		else maxLevel = parentNode.level;
	}

	childIdx = mtGetParentIndices(bottomNode.level-maxLevel, bottomNode);

	leafNode.index = mtGetParentIndices(bottomNode.level-minLevel, bottomNode);
	leafNode.level = minLevel;
}

__device__ int mtBinarySearchAlongPath(const unsigned int topLevel, const MTNode &bottomNode, const float searchValue,
										MTNode &leafNode, uint3 &childIdx) {
	unsigned int minLevel = topLevel, maxLevel = bottomNode.level;

	MTNode parentNode;
	while(maxLevel-minLevel > 1) {
		parentNode.level = (maxLevel+minLevel)>>1;

		parentNode.index = mtGetParentIndices(bottomNode.level-parentNode.level, bottomNode);

		int nodeIsValid = mtNodeExists(parentNode);
		if (nodeIsValid) {
			float2 minMax = mtGetMinMax(parentNode);
			if (minMax.x /*min*/ > searchValue || minMax.y /*max*/ < searchValue) nodeIsValid = 0;
		}

		if (nodeIsValid) minLevel = parentNode.level;
		else maxLevel = parentNode.level;
	}

	childIdx = mtGetParentIndices(bottomNode.level-maxLevel, bottomNode);

	leafNode.index = mtGetParentIndices(bottomNode.level-minLevel, bottomNode);
	leafNode.level = minLevel;

	float2 minMax = mtGetMinMax(leafNode);

	if (minMax.x /*min*/ <= searchValue && minMax.y /*max*/ >= searchValue) return 1;
	else return 0;
}

__device__ uint3 mtGetNodeIndex(const float3 pos, const unsigned int level) {
	float3 invLeafNodeSize = matricesInvCellSizes[NODE_MATRIX_ID][level];
	uint3 index;
	index.x = __fmul_rz(pos.x - mtMinPos.x, invLeafNodeSize.x);
	index.y = __fmul_rz(pos.y - mtMinPos.y, invLeafNodeSize.y);
	index.z = __fmul_rz(pos.z - mtMinPos.z, invLeafNodeSize.z);
	return index;
}

template <class T>
__device__ int mtFindLeaf(float3 pos, MTNode &leafNode, texture<T, 2, cudaReadModeElementType> tex, const T nullVal) {
	if (!isPointInBox(pos, mtMinPos, mtMaxPos)) return 0;

	leafNode.level = matricesLevels[NODE_MATRIX_ID]-1;
	leafNode.index = mtGetNodeIndex(pos, leafNode.level);

	if (mtNodeExists(leafNode, tex, nullVal)) return 1;

	mtBinarySearchAlongPath(0, leafNode, leafNode, tex, nullVal);

	return 1;
}

//__device__ int mtFindLeaf(float3 pos, MTNode &leafNode, uint2 &vectorOffsets) {
//	if (!isPointInBox(pos, mtMinPos, mtMaxPos)) return 0;
//
//	leafNode.level = matricesLevels[NODE_MATRIX_ID]-1;
//	leafNode.index = mtGetNodeIndex(pos, leafNode.level);
//
//	if (mtNodeExists(leafNode, vectorOffsets)) return 1;
//
//	mtBinarySearchAlongPath(0, leafNode, leafNode, vectorOffsets);
//
//	return 1;
//}

__device__ int mtFindLeaf(float3 pos, MTNode &leafNode, uint3 &childIdx) {
	if (!isPointInBox(pos, mtMinPos, mtMaxPos)) return 0;

	leafNode.level = matricesLevels[NODE_MATRIX_ID]-1;
	leafNode.index = mtGetNodeIndex(pos, leafNode.level);

	if (mtNodeExists(leafNode)) return 1;

	mtBinarySearchAlongPath(0, leafNode, leafNode, childIdx);

	return 1;
}

__device__ int mtFindLeaf(const float3 pos, const float searchValue, MTNode &leafNode, uint3 &childIdx) {
	if (!isPointInBox(pos, mtMinPos, mtMaxPos)) return 0;

	leafNode.level = matricesLevels[NODE_MATRIX_ID]-1;
	leafNode.index = mtGetNodeIndex(pos, leafNode.level);

	if (mtNodeExists(leafNode)) {
		float2 minMax = mtGetMinMax(leafNode);
		if (minMax.x /*min*/ <= searchValue && minMax.y /*max*/ >= searchValue) return 1;
	}

	return mtBinarySearchAlongPath(0, leafNode, searchValue, leafNode, childIdx);
}


// returns whether the point lies inside a neighbor
//   0 = same cell, 1 = neighbor cell, -1 = outside structure
// computes all 27 touching cells, not just the 8 ones sharing a face 
__device__ int mtGetNeighborRelativeToPos(MTNode &node, const float3 worldPos) {
	float3 nodeMin(matricesGetCellMinPos(node.index, node.level, NODE_MATRIX_ID));
	float3 nodeMax(nodeMin + matricesCellSizes[NODE_MATRIX_ID][node.level]);
	float3 dim = matricesCellSizes[NODE_MATRIX_ID][node.level];

	uint3 numCells = matricesNumberOfCells[NODE_MATRIX_ID][node.level]-make_uint3(1);

	int moved = 0;
	if (worldPos.x < nodeMin.x) {
		do { if (node.index.x == 0)			 return -1;  else node.index.x--;  nodeMin.x -= dim.x; } while (worldPos.x < nodeMin.x);
		moved = 1;
	}
	else if (worldPos.x > nodeMax.x) {
		do { if (node.index.x == numCells.x) return -1;  else node.index.x++;  nodeMax.x += dim.x; } while (worldPos.x > nodeMax.x);
		moved = 1;
	}

	if (worldPos.y < nodeMin.y) {
		do { if (node.index.y == 0)			 return -1;  else node.index.y--;  nodeMin.y -= dim.y; } while (worldPos.y < nodeMin.y);
		moved = 1;
	}
	else if (worldPos.y > nodeMax.y) {
		do { if (node.index.y == numCells.y) return -1;  else node.index.y++;  nodeMax.y += dim.y; } while (worldPos.y > nodeMax.y);
		moved = 1;
	}

	if (worldPos.z < nodeMin.z) {
		do { if (node.index.z == 0)			 return -1;  else node.index.z--;  nodeMin.z -= dim.z; } while (worldPos.z < nodeMin.z);
		moved = 1;
	}
	else if (worldPos.z > nodeMax.z) {
		do { if (node.index.z == numCells.z) return -1;  else node.index.z++;  nodeMax.z += dim.z; } while (worldPos.z > nodeMax.z);
		moved = 1;
	}

	return moved;
}

// provided node should exist, pt is used to determine the path
__device__ void mtSearchDownTreeToFindLeaf_Linear(const float3 pt, MTNode &leaf, uint3 &childIdx) {
	MTNode child;
	
	unsigned int leafLevel = matricesLevels[NODE_MATRIX_ID]-1;
	//#pragma unroll 16
	while (leaf.level < leafLevel) {
		int childID = mtGetChildThatPosLiesIn(pt, leaf);
		child.index = mtGetChildIndices(leaf, childID);
		child.level = leaf.level+1;

		if (!mtNodeExists(child)) break;

		leaf.index = child.index;
		leaf.level = child.level;
	}

	childIdx = child.index;
}

__device__ void mtSearchDownTreeToFindLeaf_Linear(const float3 pt, const float searchValue, MTNode &leaf, uint3 &childIdx) {
	MTNode child;
	
	unsigned int leafLevel = matricesLevels[NODE_MATRIX_ID]-1;
	//#pragma unroll 16
	while (leaf.level < leafLevel) {
		int childID = mtGetChildThatPosLiesIn(pt, leaf);
		child.index = mtGetChildIndices(leaf, childID);
		child.level = leaf.level+1;

		if (!mtNodeExists(child)) break;

		float2 minMax = mtGetMinMax(child);
		if (minMax.x /*min*/ > searchValue || minMax.y /*max*/ < searchValue) break;

		leaf.index = child.index;
		leaf.level = child.level;
	}

	childIdx = child.index;
}

__device__ void mtSearchDownTreeToFindMax_Linear(const float3 pt, const float minValue, MTNode &leaf, uint3 &childIdx) {
	MTNode child;
	
	unsigned int leafLevel = matricesLevels[NODE_MATRIX_ID]-1;
	while (leaf.level < leafLevel) {
		int childID = mtGetChildThatPosLiesIn(pt, leaf);
		child.index = mtGetChildIndices(leaf, childID);
		child.level = leaf.level+1;

		if (!mtNodeExists(child)) break;

		float2 minMax = mtGetMinMax(child);
		if (minMax.y /*max*/ <= minValue) break;

		leaf.index = child.index;
		leaf.level = child.level;
	}

	childIdx = child.index;
}

__device__ void mtSearchUpTreeToFindLeaf_Linear(MTNode &node, uint3 &childIdx) {
	childIdx = node.index;

	node.index = mtGetParentIndices(node);
	node.level--;

	#pragma unroll 16
	while (!mtNodeExists(node)) {
		childIdx = node.index;
		
		node.index = mtGetParentIndices(node);
		node.level--;
	}
}

// to do: break up into seperate functions like with what I did for MIP
__device__ void mtSearchUpTreeToFindLeaf_Linear(MTNode &node, const float searchValue, uint3 &childIdx) {
	childIdx = node.index;

	node.index = mtGetParentIndices(node);
	node.level--;

	// find first actual node
	#pragma unroll 16
	while (!mtNodeExists(node)) {
		childIdx = node.index;
		
		node.index = mtGetParentIndices(node);
		node.level--;
	}

	// now make sure that node is within the min/max range
	float2 minMax = mtGetMinMax(node);
	if (minMax.x /*min*/ <= searchValue && minMax.y /*max*/ >= searchValue) return;

	while (node.level > 0) {
		childIdx = node.index;

		node.index = mtGetParentIndices(node);
		node.level--;

		float2 minMax = mtGetMinMax(node);
		if (minMax.x /*min*/ <= searchValue && minMax.y /*max*/ >= searchValue) break;
	}
}

// passed leaf must exist and leaf.max <= minValue
__device__ void searchUpTreeToFindMax_Linear(MTNode &node, const float minValue, uint3 &childIdx) {
	while (node.level > 0) {
		childIdx = node.index;

		node.index = mtGetParentIndices(node);
		node.level--;

		float2 minMax = mtGetMinMax(node);
		if (minMax.y /*max*/ > minValue) break;
	}
}


__device__ float mtGetActualCornerData(const MTNode &cornerNode) {
	unsigned int idx = blockGetDataIndex(cornerNode.index, cornerNode.level, CORNER_MATRIX_ID);
	return cudaTexFetch(matricesGetFloatDataTexture(CORNER_MATRIX_ID), idx);
}

__device__ float mtGetCornerData(const MTNode &node, const int cornerID) {
	MTNode cornerNode;
	cornerNode.level = node.level;
	cornerNode.index.x = node.index.x;
	cornerNode.index.y = node.index.y;
	cornerNode.index.z = node.index.z;
	mtGetCornerIndex(cornerNode, cornerID);

	return mtGetActualCornerData(cornerNode);
}

__device__ void mtGetCornerIndex(MTNode &node) {
	// "base" case that will result in an infinite loop if we don't check for it
	if (node.index.x == 0 && node.index.y == 0 && node.index.z == 0) { node.level = 0;  return; }

	// if all the indices are even, that means we need to go a level higher to find the correct place
	#pragma unroll 16
	while ((node.index.x&1)==0 && (node.index.y&1)==0 && (node.index.z&1)==0) {
		node.level--;
		node.index.x = node.index.x>>1;
		node.index.y = node.index.y>>1;
		node.index.z = node.index.z>>1;
	}
}

__device__ float mtGetCornerData(const unsigned int x, const unsigned int y, const unsigned int z) {
	MTNode cornerNode;
	cornerNode.level = matricesLevels[CORNER_MATRIX_ID]-1;
	cornerNode.index.x = x;
	cornerNode.index.y = y;
	cornerNode.index.z = z;

	return mtGetActualCornerData(cornerNode);
}

// assumes that all values exist
__device__ void mtGet2x2x2CornerBlock(const uint3 lowIndex,
									  float &val000, float &val001, float &val010, float &val011,
									  float &val100, float &val101, float &val110, float &val111) {
	texture<float, 2, cudaReadModeElementType> tex = matricesGetFloatDataTexture(CORNER_MATRIX_ID);
	blockGet2x2x2DataBlock(lowIndex, 0, CORNER_MATRIX_ID, tex, val000, val001, val010, val011, val100, val101, val110, val111);
}

__device__ float3 mtGetCellSpacePosForCornerValueInterpolation(const float3 worldPos) {
	unsigned int level = matricesLevels[NODE_MATRIX_ID]-1;
	float3 mtNodeGridMinPos = mtMinPos - matricesCellSizes[NODE_MATRIX_ID][level]*0.5f;
	return ((worldPos-mtNodeGridMinPos)*matricesInvCellSizes[NODE_MATRIX_ID][level])-0.5f;
}

__device__ float3 mtGetTValueForCornerValueInterpolation(const float3 worldPos, const uint3 index) {
	float3 cellSpacePos = mtGetCellSpacePosForCornerValueInterpolation(worldPos);
	return make_float3(index.x+1, index.y+1, index.z+1) - cellSpacePos;
}


//__device__ float getLinearFromCornerData(const float3 pos) {
//	unsigned int level = matricesLevels[NODE_MATRIX_ID]-1;
//	float3 mtNodeGridMinPos = mtMinPos - matricesCellSizes[NODE_MATRIX_ID][level]*0.5f;
//	float3 cellSpacePos = ((pos-mtNodeGridMinPos)*matricesInvCellSizes[NODE_MATRIX_ID][level])-0.5f;
//
//	// cellSpacePos will be a float [0, numberCells-1]
//	// the integers that bound that number will be the corner positions
//	// the result seems to be the same if they are computed this way as opposed using the getBoundingIntegers() function
//
//	int i0 = truncf(cellSpacePos.x), j0 = truncf(cellSpacePos.y), k0 = truncf(cellSpacePos.z);
//	int i1 = i0+1, j1 = j0+1, k1 = k0+1;
//
//	//if (mtCornerLevels == 1) {
//		// since everything is on the same level, we can exploit contiguous memory for faster lookups
//		float val000, val001, val010, val011, val100, val101, val110, val111;
//		mtGet2x2CornerBlock(i0, j0, k0, val000, val001, val010, val011, val100, val101, val110, val111);
//
//		return trilinearLerp(val000, val001, val010, val011, val100, val101, val110, val111,
//							 i1-cellSpacePos.x, j1-cellSpacePos.y, k1-cellSpacePos.z);
//	//}
//	//else {
//	//	// using a tolerance and checking to see if i0==i1, j0==j1, k0==k1 actually makes this perform slower due to the branching
//	//	float val000 = mtGetCornerData(i0,j0,k0);
//	//	float val001 = mtGetCornerData(i0,j0,k1);
//	//	float val010 = mtGetCornerData(i0,j1,k0);
//	//	float val011 = mtGetCornerData(i0,j1,k1);
//	//	float val100 = mtGetCornerData(i1,j0,k0);
//	//	float val101 = mtGetCornerData(i1,j0,k1);
//	//	float val110 = mtGetCornerData(i1,j1,k0);
//	//	float val111 = mtGetCornerData(i1,j1,k1);
//
//	//	return trilinearLerp(val000, val001, val010, val011, val100, val101, val110, val111,
//	//						 i1-cellSpacePos.x, j1-cellSpacePos.y, k1-cellSpacePos.z);
//	//}
//}




__device__ int verifyTexture(const texture<float2,2,cudaReadModeElementType> tex, const unsigned int textureIdx, const float2 expectedValue) {
	float2 val = cudaTexFetch(tex, textureIdx);
	if (val.x == expectedValue.x && val.y == expectedValue.y) return 1;
	else return 0;
}

__global__ void mtVerifyMinMaxDataKernel(const uint4 *nodes, const unsigned int numberNodes, const float2 *expectedValues,
										  const unsigned int *expectedIndex, int *dataVerified) {
	unsigned int idx = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

	if (idx >= numberNodes) return;

	dataVerified[idx] = 0;

	MTNode node;
	node.index.x = nodes[idx].x;
	node.index.y = nodes[idx].y;
	node.index.z = nodes[idx].z;
	node.level = nodes[idx].w;

	int dataIndex = blockGetDataIndex(node.index, node.level, NODE_MATRIX_ID);
	if (dataIndex != expectedIndex[idx]) {
		cudaPrint("%d %d\n", dataIndex, expectedIndex[idx]);
		dataVerified[idx] += 1;
	}

	float2 minMax = mtGetMinMax(dataIndex);
	if (minMax.x != expectedValues[idx].x || minMax.y != expectedValues[idx].y) {
		cudaPrint("%d %d %d %d\n", node.index.x, node.index.y, node.index.z, node.level);
		cudaPrint("\t%d = %f %f\n", dataIndex, minMax.x, minMax.y);
		cudaPrint("\t%d = %f %f\n", expectedIndex[idx], expectedValues[idx].x, expectedValues[idx].y);
		dataVerified[idx] += 2;
	}
}

extern "C"
void mtVerifyMinMaxData(const uint4 *nodes, const unsigned int numberNodes, const float2 *expectedValues, const unsigned int *expectedIndex) {
	dim3 blockSize = dim3(64, 1, 1);
	dim3 gridSize = dim3((numberNodes % blockSize.x != 0) ? (numberNodes / blockSize.x + 1) : (numberNodes / blockSize.x));

	uint4 *deviceNodes;  cutilSafeCall(cudaMalloc((void**)&deviceNodes, numberNodes*sizeof(uint4)));
	cutilSafeCall(cudaMemcpy(deviceNodes, nodes, numberNodes*sizeof(uint4), cudaMemcpyHostToDevice));

	float2 *deviceExpectedValues;  cutilSafeCall(cudaMalloc((void**)&deviceExpectedValues, numberNodes*sizeof(float2)));
	cutilSafeCall(cudaMemcpy(deviceExpectedValues, expectedValues, numberNodes*sizeof(float2), cudaMemcpyHostToDevice));

	unsigned int *deviceExpectedIndex;  cutilSafeCall(cudaMalloc((void**)&deviceExpectedIndex, numberNodes*sizeof(unsigned int)));
	cutilSafeCall(cudaMemcpy(deviceExpectedIndex, expectedIndex, numberNodes*sizeof(unsigned int), cudaMemcpyHostToDevice));

	int *deviceOutput;  cutilSafeCall(cudaMalloc((void**)&deviceOutput, numberNodes*sizeof(int)));

	mtVerifyMinMaxDataKernel<<<gridSize, blockSize>>>(nodes, numberNodes, expectedValues, expectedIndex, deviceOutput);

	int *output;  cutilSafeCall(cudaMallocHost((void**)&output, numberNodes*sizeof(int)));
	cutilSafeCall(cudaMemcpy(output, deviceOutput, numberNodes*sizeof(int), cudaMemcpyDeviceToHost));

	for (unsigned int i=0; i<numberNodes; i++) {
		if (output[i] != 0) printf("here\n");
	}

	cutilSafeCall(cudaFree(deviceNodes));
	cutilSafeCall(cudaFree(deviceExpectedValues));
	cutilSafeCall(cudaFree(deviceExpectedIndex));
	cutilSafeCall(cudaFree(deviceOutput));
}

#endif