#ifndef CUDA_MATRIX_CSR_CUH
#define CUDA_MATRIX_CSR_CUH

#include "CudaMatrices.cuh"

// to do: some of the csr matrix type stuff has not been tested since moved from the matrix tree code and into here

__constant__ unsigned int matricesRowColIndexOffsets[CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL+1];
__constant__ unsigned int matricesDepthIndexOffsets [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL+1];

cudaArray                                           *matricesRowColIndices0 = 0, *matricesDepthIndices0 = 0;
texture<unsigned int, 2, cudaReadModeElementType> matricesRowColIndicesTex0,   matricesDepthIndicesTex0;

cudaArray                                           *matricesRowColIndices1 = 0, *matricesDepthIndices1 = 0;
texture<unsigned int, 2, cudaReadModeElementType> matricesRowColIndicesTex1,   matricesDepthIndicesTex1;

cudaArray                                           *matricesRowColIndices2 = 0, *matricesDepthIndices2 = 0;
texture<unsigned int, 2, cudaReadModeElementType> matricesRowColIndicesTex2,   matricesDepthIndicesTex2;


__device__ __host__ t_uint_2 csrGetRowColIndexTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesRowColIndicesTex, matrixID) }
__device__ __host__ t_uint_2* csrGetRowColIndexTexturePtr(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesRowColIndicesTex, matrixID) }
__device__ __host__ cudaArray* csrGetRowColIndexArray(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesRowColIndices, matrixID) }

__device__ __host__ t_uint_2 csrGetDepthIndexTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDepthIndicesTex, matrixID) }
__device__ __host__ t_uint_2* csrGetDepthIndexTexturePtr(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesDepthIndicesTex, matrixID) }
__device__ __host__ cudaArray* csrGetDepthIndexArray(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDepthIndices, matrixID) }


extern "C"
void matricesLoadMatrixIndicesToCuda(const CudaMatrixIDs matrixID, unsigned int levels, unsigned int *rowColIndexOffsets,
									 unsigned int *depthIndexOffsets, unsigned int *rowColIndices, unsigned int *depthIndices) {
	if (levels >= CUDA_MAX_MATRIX_LEVEL) printf("Number of matrices too large, increase #define CUDA_MAX_MATRIX_LEVEL\n");
	if (matrixID >= CUDA_MAX_MATRIX_GROUPS) printf("Matrix ID too large, increase #define CUDA_MAX_MATRIX_GROUPS\n");

	cutilSafeCall(cudaMemcpyToSymbol(matricesRowColIndexOffsets, rowColIndexOffsets, (levels+1)*sizeof(unsigned int), matrixID*(CUDA_MAX_MATRIX_LEVEL+1)*sizeof(unsigned int)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesDepthIndexOffsets, depthIndexOffsets, (levels+1)*sizeof(unsigned int), matrixID*(CUDA_MAX_MATRIX_LEVEL+1)*sizeof(unsigned int)));

	cudaArray *rowColArray = csrGetRowColIndexArray(matrixID);
	cudaArray *depthArray = csrGetDepthIndexArray(matrixID);

	texture<unsigned int, 2, cudaReadModeElementType> *rowColTex = csrGetRowColIndexTexturePtr(matrixID);
	texture<unsigned int, 2, cudaReadModeElementType> *depthTex = csrGetDepthIndexTexturePtr(matrixID);

	if (rowColIndexOffsets[levels] > 0) loadVectorToTexture2D(rowColIndices, rowColIndexOffsets[levels], &rowColArray, rowColTex);
	if (depthIndexOffsets[levels] > 0) loadVectorToTexture2D(depthIndices, depthIndexOffsets[levels], &depthArray, depthTex);
}

__device__ unsigned int csrMapToRowColIndex(const unsigned int x, const unsigned int y, const uint3 numCells) {
	return __umul24(x,numCells.y) + y;
	//return interleaveBits(x,y);
}
__device__ unsigned int csrGetRightRowColIndex(const unsigned int idx, const uint3 numCells) { return idx+numCells.y; }
__device__ unsigned int csrGetLeftRowColIndex(const unsigned int idx, const uint3 numCells) { return idx-numCells.y; }
__device__ unsigned int csrGetUpRowColIndex(const unsigned int idx) { return idx+1; }
__device__ unsigned int csrGetDownRowColIndex(const unsigned int idx) { return idx-1; }


// returns the index relative to the start of the offset for that level, -1 if not found
__device__ int csrGetIndex(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	unsigned int xyIndex = csrMapToRowColIndex(index.x, index.y, matricesNumberOfCells[matrixID][level]);
	unsigned int rowColOffset = matricesRowColIndexOffsets[matrixID][level];
	unsigned int depthOffset = matricesDepthIndexOffsets[matrixID][level];
	
	// the supposedly more efficient texture fetch that gets multiple neighboring values actually makes this slower
	texture<unsigned int, 2, cudaReadModeElementType> rowColTex = csrGetRowColIndexTexture(matrixID);
	unsigned int rowColStart = cudaTexFetch(rowColTex, rowColOffset+xyIndex)+depthOffset;
	unsigned int rowColEnd = cudaTexFetch(rowColTex, rowColOffset+xyIndex+1)+depthOffset;

	int idx = binarySearch_sortedLowToHigh(csrGetDepthIndexTexture(matrixID), index.z, rowColStart, rowColEnd);
	
	// and again, for whatever reason, doing this if-else test results in better performance than just returning idx-depthOffset
	if (idx >= 0) return idx-depthOffset;
	else return -1;
}

__device__ int csrGetDataIndex(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	uint3 numCells = matricesNumberOfCells[matrixID][level];
	unsigned int xyIndex = csrMapToRowColIndex(index.x, index.y, numCells);
	unsigned int rowColOffset = matricesRowColIndexOffsets[matrixID][level];
	unsigned int depthOffset = matricesDepthIndexOffsets[matrixID][level];
	
	// the supposedly more efficient texture fetch that gets multiple neighboring values actually makes this slower
	texture<unsigned int, 2, cudaReadModeElementType> rowColTex = csrGetRowColIndexTexture(matrixID);
	unsigned int rowColStart = cudaTexFetch(rowColTex, rowColOffset+xyIndex)+depthOffset;
	unsigned int rowColEnd = cudaTexFetch(rowColTex, rowColOffset+xyIndex+1)+depthOffset;

	// this finds the index but does not do the last test at the end, the other will return -1 if it fails which we do not want
	int idx = binarySearch_sortedLowToHigh_guaranteedToBeFound(csrGetDepthIndexTexture(matrixID), index.z, rowColStart, rowColEnd);

	idx -= depthOffset;

	int type = matricesMatrixTypes[matrixID][level];
	if (type == 1) return idx + matricesDataOffsets[matrixID][level]; // GRID_REGULAR_CSR or GRID_REGULAR_SPARSE_NODATA, normal mapping
	else return numCells.z*xyIndex + index.z - idx + matricesDataOffsets[matrixID][level]; // non-normal mapping
}

template <class T> __device__ int csrCellExists(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	int type = matricesMatrixTypes[matrixID][level];
	int idx = csrGetIndex(index, level, matrixID);
	
	// type = 1: GRID_REGULAR_CSR or GRID_REGULAR_SPARSE_NODATA, normal mapping
	//    so found index means that the node exists
	// type = 2: GRID_REGULAR_SPARSE_NODATA, inverse mapping, which means a found index indicates a non-existent node
	if (idx >= 0) {
		if (type == 1) return 1;
		else return 0;
	}
	else {
		if (type == 1) return 0;
		else return 1;
	}
}

template <class T> __device__ void csrGet2x2x2DataBlock(const uint3 lowIndex, const unsigned int level, const CudaMatrixIDs matrixID,
														texture<T, 2, cudaReadModeElementType> tex,
														float &val000, float &val001, float &val010, float &val011,
														float &val100, float &val101, float &val110, float &val111) {
	uint3 numCells = matricesNumberOfCells[matrixID][level];
	unsigned int xyIndex = csrMapToRowColIndex(lowIndex.x, lowIndex.y, numCells);
	
	texture<unsigned int, 2, cudaReadModeElementType> rowColTex = csrGetRowColIndexTexture(matrixID);
	unsigned int rowColStart = cudaTexFetch(rowColTex, xyIndex);
	unsigned int rowColEnd = cudaTexFetch(rowColTex, xyIndex+1);

	texture<unsigned int, 2, cudaReadModeElementType> depthTex = csrGetDepthIndexTexture(matrixID);
	unsigned int depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, rowColEnd);
	val000 = cudaTexFetch(tex, depth);
	val001 = cudaTexFetch(tex, depth+1);  // we can assume that the next value in the texture is the one we want
	
	// chances are the position of what we just found will be close to the other rows, so use it to speed up the search
	//   worst case, this will basically increase the numbers of iterations of the binary search by 1
	unsigned int offset = depth-rowColStart;

	xyIndex = csrGetUpRowColIndex(xyIndex);
	rowColStart = cudaTexFetch(rowColTex, xyIndex);
	rowColEnd = cudaTexFetch(rowColTex, xyIndex+1);
	unsigned int dim = rowColEnd-rowColStart;


	if (offset < dim) {
		depth = rowColStart+offset;
		unsigned int initialGuess = cudaTexFetch(depthTex, depth);

		if (lowIndex.z < initialGuess) depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, depth);
		else if (lowIndex.z > initialGuess) depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, depth, rowColEnd);
	}
	else depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, rowColEnd);

	val010 = cudaTexFetch(tex, depth);
	val011 = cudaTexFetch(tex, depth+1);


	offset = depth-rowColStart;
	xyIndex = csrGetRightRowColIndex(xyIndex, numCells);
	rowColStart = cudaTexFetch(rowColTex, xyIndex);
	rowColEnd = cudaTexFetch(rowColTex, xyIndex+1);
	dim = rowColEnd-rowColStart;

	if (offset < dim) {
		depth = rowColStart+offset;
		unsigned int initialGuess = cudaTexFetch(depthTex, depth);

		if (lowIndex.z < initialGuess) depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, depth);
		else if (lowIndex.z > initialGuess) depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, depth, rowColEnd);
	}
	else depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, rowColEnd);
	
	val110 = cudaTexFetch(tex, depth);
	val111 = cudaTexFetch(tex, depth+1);


	offset = depth-rowColStart;
	xyIndex = csrGetDownRowColIndex(xyIndex);
	rowColStart = cudaTexFetch(rowColTex, xyIndex);
	rowColEnd = cudaTexFetch(rowColTex, xyIndex+1);
	dim = rowColEnd-rowColStart;

	if (offset < dim) {
		depth = rowColStart+offset;
		unsigned int initialGuess = cudaTexFetch(depthTex, depth);

		if (lowIndex.z < initialGuess) depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, depth);
		else if (lowIndex.z > initialGuess) depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, depth, rowColEnd);
	}
	else depth = binarySearch_sortedLowToHigh_guaranteedToBeFound(depthTex, lowIndex.z, rowColStart, rowColEnd);
	
	val100 = cudaTexFetch(tex, depth);
	val101 = cudaTexFetch(tex, depth+1);
}

extern "C"
void cudaMatrixCSRFreeData() {
	for (unsigned int i=0; i<CUDA_MAX_MATRIX_GROUPS; i++) {
		CudaMatrixIDs matrixID = matrixConvertIntToIDEnum(i);

		cudaArray *rowColArray = csrGetRowColIndexArray(matrixID);
		cudaArray *depthArray = csrGetDepthIndexArray(matrixID);

		if (rowColArray != NULL) cutilSafeCall(cudaFreeArray(rowColArray));
		if (depthArray != NULL) cutilSafeCall(cudaFreeArray(depthArray));
	}
}

#endif