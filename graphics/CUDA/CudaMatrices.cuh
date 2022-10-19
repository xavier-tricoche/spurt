#ifndef CUDA_MATRICES_CUH
#define CUDA_MATRICES_CUH

#include "CudaConstants.cuh"

// to do: some of the full and csr matrix type stuff has not been tested since moved from the matrix tree code and into here

#define CUDA_MAX_MATRIX_LEVEL 32
#define CUDA_MAX_MATRIX_GROUPS 3

#define matrixArrayFetch(_texbasename, _texnum)\
switch(_texnum) {\
	case 0: return _texbasename##0;\
	case 1: return _texbasename##1;\
	default: return _texbasename##2;\
}

#define matrixArrayPtrFetch(_texbasename, _texnum)\
switch(_texnum) {\
	case 0: return &_texbasename##0;\
	case 1: return &_texbasename##1;\
	default: return &_texbasename##2;\
}

#define t_char_2 texture<char, 2, cudaReadModeElementType>
#define t_uint_2 texture<unsigned int, 2, cudaReadModeElementType>
#define t_float_2 texture<float, 2, cudaReadModeElementType>
#define t_float2_2 texture<float2, 2, cudaReadModeElementType>

__device__ __host__ CudaMatrixIDs matrixConvertIntToIDEnum(const int matrixID) { matrixArrayFetch(CUDA_MATRIX_, matrixID); }


// ----------------------- Position and Index Info ------------

__constant__ float3 mtMinPos, mtMaxPos, mtInvDim;

extern "C"
void mtLoadPositionInfo(float* dataMinPos, float* dataMaxPos) {
	cutilSafeCall(cudaMemcpyToSymbol(mtMinPos, dataMinPos, sizeof(float3)));
	cutilSafeCall(cudaMemcpyToSymbol(mtMaxPos, dataMaxPos, sizeof(float3)));

	float invDim[3] = {1.0f/(dataMaxPos[0]-dataMinPos[0]), 1.0f/(dataMaxPos[1]-dataMinPos[1]), 1.0f/(dataMaxPos[2]-dataMinPos[2])};
	cutilSafeCall(cudaMemcpyToSymbol(mtInvDim, invDim, sizeof(float3)));
}

__constant__ unsigned int matricesLevels            [CUDA_MAX_MATRIX_GROUPS];
__constant__ int		  matricesMatrixTypes       [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL];

__constant__ unsigned int matricesDataOffsets       [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL+1];

__constant__ uint3		  matricesNumberOfCells     [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL];
__constant__ float3		  matricesCellSizes         [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL];
__constant__ float3		  matricesInvCellSizes      [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL];

__constant__ unsigned int matricesBlockIndexOffsets [CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL+1];
__constant__ uint3        matricesBlockNumberOfCells[CUDA_MAX_MATRIX_GROUPS];
__constant__ float3       matricesBlockInvNumberOfCells[CUDA_MAX_MATRIX_GROUPS];
__constant__ uint3        matricesBlockOverlapNumberOfCells[CUDA_MAX_MATRIX_GROUPS];
__constant__ unsigned int matricesBlockMaxCells     [CUDA_MAX_MATRIX_GROUPS];
__constant__ uint3        matricesBlockRedirectionNumberOfCells[CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL];
__constant__ unsigned int matricesBlockRedirectionMaxCells[CUDA_MAX_MATRIX_GROUPS][CUDA_MAX_MATRIX_LEVEL];


// unfortunately have to list out the texture instead of doing an array right now...
cudaArray                                           *matricesBlockIndices0 = 0, *matricesBlockIndices1 = 0, *matricesBlockIndices2 = 0;
texture<unsigned int, 2, cudaReadModeElementType> matricesBlockIndicesTex0,   matricesBlockIndicesTex1,   matricesBlockIndicesTex2;


__device__ __host__ t_uint_2 blockGetIndexTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesBlockIndicesTex, matrixID) }
t_uint_2* blockGetIndexTexturePtr(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesBlockIndicesTex, matrixID) }
cudaArray* blockGetIndexArray(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesBlockIndices, matrixID) }

extern "C"
void matricesLoadMatrixIndicesToCuda(const CudaMatrixIDs matrixID, unsigned int levels, int *matrixTypes, unsigned int *blockIndexOffsets,
									 unsigned int *blockIndices, unsigned int *numberOfCells, float *cellSizes) {
	if (levels >= CUDA_MAX_MATRIX_LEVEL) printf("Number of matrices too large, increase #define CUDA_MAX_MATRIX_LEVEL\n");
	if (matrixID >= CUDA_MAX_MATRIX_GROUPS) printf("Matrix ID too large, increase #define CUDA_MAX_MATRIX_GROUPS\n");

	cutilSafeCall(cudaMemcpyToSymbol(matricesLevels, &levels, sizeof(unsigned int), matrixID*sizeof(unsigned int)));

	cutilSafeCall(cudaMemcpyToSymbol(matricesMatrixTypes, matrixTypes, levels*sizeof(int), matrixID*CUDA_MAX_MATRIX_LEVEL*sizeof(int)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockIndexOffsets, blockIndexOffsets, (levels+1)*sizeof(unsigned int), matrixID*(CUDA_MAX_MATRIX_LEVEL+1)*sizeof(unsigned int)));
	
	float *invCellSizes = (float*)malloc(levels*sizeof(float)*3);
	for (unsigned int i=0; i<levels; i++) {
		invCellSizes[i*3 + 0] = 1.0f/cellSizes[i*3 + 0]/* - 0.000001f*/;  // the little extra added is to fix 
		invCellSizes[i*3 + 1] = 1.0f/cellSizes[i*3 + 1]/* - 0.000001f*/;  //   some floating point errors that 
		invCellSizes[i*3 + 2] = 1.0f/cellSizes[i*3 + 2]/* - 0.000001f*/;  //   occur later on
	}

	cutilSafeCall(cudaMemcpyToSymbol(matricesNumberOfCells, numberOfCells, levels*sizeof(uint3), matrixID*CUDA_MAX_MATRIX_LEVEL*sizeof(uint3)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesCellSizes, cellSizes, levels*sizeof(float3), matrixID*CUDA_MAX_MATRIX_LEVEL*sizeof(float3)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesInvCellSizes, invCellSizes, levels*sizeof(float3), matrixID*CUDA_MAX_MATRIX_LEVEL*sizeof(float3)));

	free(invCellSizes);		invCellSizes = NULL;

	cudaArray *blockArray = blockGetIndexArray(matrixID);
	texture<unsigned int, 2, cudaReadModeElementType> *blockTex = blockGetIndexTexturePtr(matrixID);

	if (blockIndexOffsets[levels] > 0) loadVectorToTexture2D(blockIndices, blockIndexOffsets[levels], &blockArray, blockTex);



	uint3 blockSize = make_uint3(8,8,8);
	float3 invBlockSize = make_float3(1.0f/(float)blockSize.x, 1.0f/(float)blockSize.y, 1.0f/(float)blockSize.z);
	uint3 overlapBlockSize = make_uint3(blockSize.x+1,blockSize.y+1,blockSize.z+1);
	if (matrixID != CORNER_MATRIX_ID) overlapBlockSize = blockSize;

	unsigned int maxCellsPerOverlapBlock = overlapBlockSize.x*overlapBlockSize.y*overlapBlockSize.z;

	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockNumberOfCells, &blockSize, sizeof(uint3), matrixID*sizeof(uint3)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockInvNumberOfCells, &invBlockSize, sizeof(float3), matrixID*sizeof(float3)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockOverlapNumberOfCells, &overlapBlockSize, sizeof(uint3), matrixID*sizeof(uint3)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockMaxCells, &maxCellsPerOverlapBlock, sizeof(unsigned int), matrixID*sizeof(unsigned int)));

	uint3 *numRedirectionCells = (uint3*)malloc(levels*sizeof(uint3));
	unsigned int *maxRedirectionCells = (unsigned int*)malloc(levels*sizeof(unsigned int));
	for (unsigned int i=0; i<levels; i++) {
		float3 numCells = make_float3((float)numberOfCells[i*3+0], (float)numberOfCells[i*3+1], (float)numberOfCells[i*3+2]);
		float3 size = make_float3((float)blockSize.x, (float)blockSize.y, (float)blockSize.z);

		numRedirectionCells[i] = make_uint3((unsigned int)ceil(numCells.x/size.x), (unsigned int)ceil(numCells.y/size.y), (unsigned int)ceil(numCells.z/size.z));
		maxRedirectionCells[i] = numRedirectionCells[i].x*numRedirectionCells[i].y*numRedirectionCells[i].z;
	}

	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockRedirectionNumberOfCells, numRedirectionCells, levels*sizeof(uint3), matrixID*CUDA_MAX_MATRIX_LEVEL*sizeof(uint3)));
	cutilSafeCall(cudaMemcpyToSymbol(matricesBlockRedirectionMaxCells, maxRedirectionCells, levels*sizeof(unsigned int), matrixID*CUDA_MAX_MATRIX_LEVEL*sizeof(unsigned int)));

	free(numRedirectionCells);		numRedirectionCells = NULL;
	free(maxRedirectionCells);		maxRedirectionCells = NULL;
}



// ---------------------- Data ------------------------

bool mtLoadDataLevels(unsigned int levels, unsigned int *dataOffsets, const CudaMatrixIDs matrixID) {
	cutilSafeCall(cudaMemcpyToSymbol(matricesDataOffsets, dataOffsets, (levels+1)*sizeof(unsigned int), matrixID*(CUDA_MAX_MATRIX_LEVEL+1)*sizeof(unsigned int)));
	return dataOffsets[levels] > 0;
}


cudaArray                                        *matricesData0 = 0,             *matricesData1 = 0,             *matricesData2 = 0;
cudaArray                                        *matricesVectorData0 = 0,       *matricesVectorData1 = 0,       *matricesVectorData2 = 0;

texture<char, 2, cudaReadModeElementType>         matricesDataTex_char_0,         matricesDataTex_char_1,         matricesDataTex_char_2;
texture<char, 2, cudaReadModeElementType>         matricesVectorDataTex_char_0,   matricesVectorDataTex_char_1,   matricesVectorDataTex_char_2;

texture<unsigned int, 2, cudaReadModeElementType> matricesDataTex_uint_0,         matricesDataTex_uint_1,         matricesDataTex_uint_2;
texture<unsigned int, 2, cudaReadModeElementType> matricesVectorDataTex_uint_0,   matricesVectorDataTex_uint_1,   matricesVectorDataTex_uint_2;

texture<float, 2, cudaReadModeElementType>        matricesDataTex_float_0,        matricesDataTex_float_1,        matricesDataTex_float_2;
texture<float, 2, cudaReadModeElementType>        matricesVectorDataTex_float_0,  matricesVectorDataTex_float_1,  matricesVectorDataTex_float_2;

texture<float2, 2, cudaReadModeElementType>       matricesDataTex_float2_0,       matricesDataTex_float2_1,       matricesDataTex_float2_2;
texture<float2, 2, cudaReadModeElementType>       matricesVectorDataTex_float2_0, matricesVectorDataTex_float2_1, matricesVectorDataTex_float2_2;

__constant__ char matricesNullValue_char[CUDA_MAX_MATRIX_GROUPS];
__constant__ unsigned int matricesNullValue_uint[CUDA_MAX_MATRIX_GROUPS];
__constant__ float matricesNullValue_float[CUDA_MAX_MATRIX_GROUPS];
__constant__ float2 matricesNullValue_float2[CUDA_MAX_MATRIX_GROUPS];


cudaArray* matricesGetDataArrayHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesData, matrixID) }
cudaArray* matricesGetVectorDataArrayHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorData, matrixID) }

cudaArray** matricesGetDataArrayPtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesData, matrixID) }
cudaArray** matricesGetVectorDataArrayPtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesVectorData, matrixID) }


__device__ t_char_2 matricesGetCharDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_char_, matrixID) }
t_char_2 matricesGetCharDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_char_, matrixID) }
t_char_2* matricesGetCharDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesDataTex_char_, matrixID) }

__device__ t_char_2 matricesGetCharVectorDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_char_, matrixID) }
t_char_2 matricesGetCharVectorDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_char_, matrixID) }
t_char_2* matricesGetCharVectorDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesVectorDataTex_char_, matrixID) }

__device__ t_uint_2 matricesGetUintDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_uint_, matrixID) }
t_uint_2 matricesGetUintDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_uint_, matrixID) }
t_uint_2* matricesGetUintDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesDataTex_uint_, matrixID) }

__device__ t_uint_2 matricesGetUintVectorDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_uint_, matrixID) }
t_uint_2 matricesGetUintVectorDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_uint_, matrixID) }
t_uint_2* matricesGetUintVectorDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesVectorDataTex_uint_, matrixID) }

__device__ t_float_2 matricesGetFloatDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_float_, matrixID) }
t_float_2 matricesGetFloatDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_float_, matrixID) }
t_float_2* matricesGetFloatDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesDataTex_float_, matrixID) }

__device__ t_float_2 matricesGetFloatVectorDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_float_, matrixID) }
t_float_2 matricesGetFloatVectorDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_float_, matrixID) }
t_float_2* matricesGetFloatVectorDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesVectorDataTex_float_, matrixID) }

__device__ t_float2_2 matricesGetFloat2DataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_float2_, matrixID) }
t_float2_2 matricesGetFloat2DataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesDataTex_float2_, matrixID) }
t_float2_2* matricesGetFloat2DataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesDataTex_float2_, matrixID) }

__device__ t_float2_2 matricesGetFloat2VectorDataTexture(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_float2_, matrixID) }
t_float2_2 matricesGetFloat2VectorDataTextureHost(const CudaMatrixIDs matrixID) { matrixArrayFetch(matricesVectorDataTex_float2_, matrixID) }
t_float2_2* matricesGetFloat2VectorDataTexturePtrHost(const CudaMatrixIDs matrixID) { matrixArrayPtrFetch(matricesVectorDataTex_float2_, matrixID) }


extern "C"
void mtInitializeAllTextures() {
	for (unsigned int i=0; i<CUDA_MAX_MATRIX_GROUPS; i++) {
		CudaMatrixIDs id = matrixConvertIntToIDEnum(i);

		t_char_2 *cdt = matricesGetCharDataTexturePtrHost(id), *cvdt = matricesGetCharVectorDataTexturePtrHost(id);
		t_uint_2 *uidt = matricesGetUintDataTexturePtrHost(id), *uivdt = matricesGetUintVectorDataTexturePtrHost(id);
		t_float_2 *fdt = matricesGetFloatDataTexturePtrHost(id), *fvdt = matricesGetFloatVectorDataTexturePtrHost(id);
		t_float2_2 *f2dt = matricesGetFloat2DataTexturePtrHost(id), *f2vdt = matricesGetFloat2VectorDataTexturePtrHost(id);


		cdt->filterMode = cvdt->filterMode = uidt->filterMode = uivdt->filterMode =
			fdt->filterMode = fvdt->filterMode = f2dt->filterMode = f2vdt->filterMode = cudaFilterModePoint;
	}
}

// "error C2733: second C linkage of overloaded function 'matricesLoadDataToCuda' not allowed"
//   so basically cannot overload and cannot template

// load data
extern "C" void matricesLoadDataToCuda_char(unsigned int levels, unsigned int *dataOffsets, char *data, const CudaMatrixIDs matrixID)
	{ if (mtLoadDataLevels(levels, dataOffsets, matrixID)) loadVectorToTexture2D(data, dataOffsets[levels], matricesGetDataArrayPtrHost(matrixID), matricesGetCharDataTexturePtrHost(matrixID)); }
extern "C" void matricesLoadDataToCuda_uint(unsigned int levels, unsigned int *dataOffsets, unsigned int *data, const CudaMatrixIDs matrixID)
	{ if (mtLoadDataLevels(levels, dataOffsets, matrixID)) loadVectorToTexture2D(data, dataOffsets[levels], matricesGetDataArrayPtrHost(matrixID), matricesGetUintDataTexturePtrHost(matrixID)); }
extern "C" void matricesLoadDataToCuda_float(unsigned int levels, unsigned int *dataOffsets, float *data, const CudaMatrixIDs matrixID)
	{ if (mtLoadDataLevels(levels, dataOffsets, matrixID)) loadVectorToTexture2D(data, dataOffsets[levels], matricesGetDataArrayPtrHost(matrixID), matricesGetFloatDataTexturePtrHost(matrixID)); }
extern "C" void matricesLoadDataToCuda_float2(unsigned int levels, unsigned int *dataOffsets, float2 *data, const CudaMatrixIDs matrixID)
	{ if (mtLoadDataLevels(levels, dataOffsets, matrixID)) loadVectorToTexture2D(data, dataOffsets[levels], matricesGetDataArrayPtrHost(matrixID), matricesGetFloat2DataTexturePtrHost(matrixID)); }

// load null values
extern "C" void matricesSetNullValue_char(const char nullVal, const CudaMatrixIDs matrixID) { cutilSafeCall(cudaMemcpyToSymbol(matricesNullValue_char, &nullVal, sizeof(char), matrixID*sizeof(char))); }
extern "C" void matricesSetNullValue_uint(const unsigned int nullVal, const CudaMatrixIDs matrixID) { cutilSafeCall(cudaMemcpyToSymbol(matricesNullValue_uint, &nullVal, sizeof(unsigned int), matrixID*sizeof(unsigned int))); }
extern "C" void matricesSetNullValue_float(const float nullVal, const CudaMatrixIDs matrixID) { cutilSafeCall(cudaMemcpyToSymbol(matricesNullValue_float, &nullVal, sizeof(float), matrixID*sizeof(float))); }
extern "C" void matricesSetNullValue_float2(const float2 nullVal, const CudaMatrixIDs matrixID) { cutilSafeCall(cudaMemcpyToSymbol(matricesNullValue_float2, &nullVal, sizeof(float2), matrixID*sizeof(float2))); }

// load vector data
extern "C" void matricesLoadVectorDataToCuda_char(unsigned int vectorDataSize, char *vectorData, const CudaMatrixIDs matrixID)
	{ if (vectorDataSize > 0) loadVectorToTexture2D(vectorData, vectorDataSize, matricesGetVectorDataArrayPtrHost(matrixID), matricesGetCharVectorDataTexturePtrHost(matrixID)); }
extern "C" void matricesLoadVectorDataToCuda_uint(unsigned int vectorDataSize, unsigned int *vectorData, const CudaMatrixIDs matrixID)
	{ if (vectorDataSize > 0) loadVectorToTexture2D(vectorData, vectorDataSize, matricesGetVectorDataArrayPtrHost(matrixID), matricesGetUintVectorDataTexturePtrHost(matrixID)); }
extern "C" void matricesLoadVectorDataToCuda_float(unsigned int vectorDataSize, float *vectorData, const CudaMatrixIDs matrixID)
	{ if (vectorDataSize > 0) loadVectorToTexture2D(vectorData, vectorDataSize, matricesGetVectorDataArrayPtrHost(matrixID), matricesGetFloatVectorDataTexturePtrHost(matrixID)); }




// ----------------- Position Info -------------------------

__device__ float3 matricesGetCellCenter(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	return make_float3(((float)index.x + 0.5f)*matricesCellSizes[matrixID][level].x + mtMinPos.x,
					   ((float)index.y + 0.5f)*matricesCellSizes[matrixID][level].y + mtMinPos.y,
					   ((float)index.z + 0.5f)*matricesCellSizes[matrixID][level].z + mtMinPos.z);
}

__device__ float3 matricesGetCellMinPos(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	return make_float3(((float)index.x)*matricesCellSizes[matrixID][level].x + mtMinPos.x,
					   ((float)index.y)*matricesCellSizes[matrixID][level].y + mtMinPos.y,
					   ((float)index.z)*matricesCellSizes[matrixID][level].z + mtMinPos.z);
}


// ------------------ Normal (Full) -------------------------
__device__ unsigned int fullMatrixGetIndex(const uint3 index, const uint3 numCells) {
	return __umul24(numCells.z, (__umul24(index.x, numCells.y) + index.y)) + index.z;
	//return interleaveBits(index.x, index.y, index.z);
}
__device__ unsigned int fullMatrixGetRightIndex(const unsigned int idx, const uint3 numCells) { return idx+__umul24(numCells.y, numCells.z); }
__device__ unsigned int fullMatrixGetLeftIndex(const unsigned int idx, const uint3 numCells) { return idx-__umul24(numCells.y, numCells.z); }
__device__ unsigned int fullMatrixGetUpIndex(const unsigned int idx, const uint3 numCells) { return idx+numCells.z; }
__device__ unsigned int fullMatrixGetDownIndex(const unsigned int idx, const uint3 numCells) { return idx-numCells.z; }
__device__ unsigned int fullMatrixGetFrontIndex(const unsigned int idx) { return idx+1; }
__device__ unsigned int fullMatrixGetBackIndex(const unsigned int idx) { return idx-1; }

__device__ int fullMatrixGetDataIndex(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	return fullMatrixGetIndex(index, matricesNumberOfCells[matrixID][level]) + matricesDataOffsets[matrixID][level];
}

template <class T> __device__ int fullMatrixCellExists(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID,
													   texture<T, 2, cudaReadModeElementType> tex, const T nullVal) {
	unsigned int idx = fullMatrixGetDataIndex(index, level, matrixID);
	if (cudaTexFetch(tex, idx) != nullVal) return 1;
	else return 0;
}

template <class T> __device__ void fullMatrixGet2x2x2DataBlock(const uint3 lowIndex, const unsigned int level, const CudaMatrixIDs matrixID,
															   texture<T, 2, cudaReadModeElementType> tex,
															   float &val000, float &val001, float &val010, float &val011,
															   float &val100, float &val101, float &val110, float &val111) {
	uint3 numberOfCells = matricesNumberOfCells[matrixID][level];

	unsigned int idx = fullMatrixGetIndex(lowIndex, numberOfCells);
	val000 = cudaTexFetch(tex, idx);
	val001 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));

	idx = fullMatrixGetUpIndex(idx, numberOfCells);
	val010 = cudaTexFetch(tex, idx);
	val011 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));

	idx = fullMatrixGetRightIndex(idx, numberOfCells);
	val110 = cudaTexFetch(tex, idx);
	val111 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));

	idx = fullMatrixGetDownIndex(idx, numberOfCells);
	val100 = cudaTexFetch(tex, idx);
	val101 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));
}




// ------------------ Block ------------------------

__device__ uint3 blockGetBlockNumber(const uint3 index, const float3 invBlockSize)
	{ return make_uint3(index.x*invBlockSize.x, index.y*invBlockSize.y, index.z*invBlockSize.z); }
	//{ return make_uint3(__fmul_rz(index.x, invBlockSize.x), __fmul_rz(index.y, invBlockSize.y), __fmul_rz(index.z, invBlockSize.z)); }
__device__ uint3 blockGetBlockNumber(const uint3 index, const CudaMatrixIDs matrixID) { return blockGetBlockNumber(index, matricesBlockInvNumberOfCells[matrixID]); }


__device__ uint3 blockGetIndexWithinBlock(const uint3 index, const uint3 block, const uint3 blockSize)
	//{ return index - block*blockSize; }
	{ return index - make_uint3(__umul24(block.x, blockSize.x), __umul24(block.y, blockSize.y), __umul24(block.z, blockSize.z)); }
__device__ uint3 blockGetIndexWithinBlock(const uint3 index, const uint3 block, const CudaMatrixIDs matrixID) { return blockGetIndexWithinBlock(index, block, matricesBlockNumberOfCells[matrixID]); }



__device__ unsigned int blockGetBlockPointer(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	unsigned int idx = fullMatrixGetIndex(index, matricesBlockRedirectionNumberOfCells[matrixID][level]) + matricesBlockIndexOffsets[matrixID][level];	
	
	texture<unsigned int, 2, cudaReadModeElementType> tex = blockGetIndexTexture(matrixID);
	return cudaTexFetch(tex, idx);
}

__device__ int blockGetDataIndex(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID) {
	uint3 blockIndex = blockGetBlockNumber(index, matrixID);

	unsigned int redirectIndex = blockGetBlockPointer(blockIndex, level, matrixID);

	return fullMatrixGetIndex(blockGetIndexWithinBlock(index, blockIndex, matrixID), matricesBlockOverlapNumberOfCells[matrixID]) +
		   __umul24(redirectIndex, matricesBlockMaxCells[matrixID]) + matricesDataOffsets[matrixID][level];
}

template <class T> __device__ int blockCellExists(const uint3 index, const unsigned int level, const CudaMatrixIDs matrixID,
												  texture<T, 2, cudaReadModeElementType> tex, const T nullVal) {
	uint3 blockIndex = blockGetBlockNumber(index, matrixID);

	unsigned int redirectIndex = blockGetBlockPointer(blockIndex, level, matrixID);

	if (redirectIndex >= matricesBlockRedirectionMaxCells[matrixID][level]) return 0;

	unsigned int idx = fullMatrixGetIndex(blockGetIndexWithinBlock(index, blockIndex, matrixID), matricesBlockOverlapNumberOfCells[matrixID]) +
					   __umul24(redirectIndex, matricesBlockMaxCells[matrixID]) + matricesDataOffsets[matrixID][level];

	if (cudaTexFetch(tex, idx) != nullVal) return 1;
	else return 0;
}

template <class T> __device__ void blockGet2x2x2DataBlock(const uint3 lowIndex, const unsigned int level, const CudaMatrixIDs matrixID,
														  texture<T, 2, cudaReadModeElementType> tex,
														  float &val000, float &val001, float &val010, float &val011,
														  float &val100, float &val101, float &val110, float &val111) {
	unsigned int idx = blockGetDataIndex(lowIndex, level, matrixID);
	uint3 overlapBlockSize = matricesBlockOverlapNumberOfCells[matrixID];

	val000 = cudaTexFetch(tex, idx);
	val001 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));

	idx = fullMatrixGetUpIndex(idx, overlapBlockSize);
	val010 = cudaTexFetch(tex, idx);
	val011 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));

	idx = fullMatrixGetRightIndex(idx, overlapBlockSize);
	val110 = cudaTexFetch(tex, idx);
	val111 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));

	idx = fullMatrixGetDownIndex(idx, overlapBlockSize);
	val100 = cudaTexFetch(tex, idx);
	val101 = cudaTexFetch(tex, fullMatrixGetFrontIndex(idx));
}


extern "C"
void matricesFreeData() {
	for (unsigned int i=0; i<CUDA_MAX_MATRIX_GROUPS; i++) {
		CudaMatrixIDs matrixID = matrixConvertIntToIDEnum(i);

		cudaArray *blockArray = blockGetIndexArray(matrixID);
		if (blockArray != NULL) cutilSafeCall(cudaFreeArray(blockArray));
	}
}



#endif