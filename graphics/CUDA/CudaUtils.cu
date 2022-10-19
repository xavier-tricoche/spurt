#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include "cutil_inline.h"
#include "cutil_math.h"

#include "CudaConstants.cuh"

#ifdef __DEVICE_EMULATION__
    #define cudaPrint(...) printf(__VA_ARGS__)
	#define cudaTex1Dfetch(...) tex1D(__VA_ARGS__) /* As of version 2.2, tex1Dfetch gives errors in emulation mode */
#else
    #define cudaPrint(...) /* Nothing */
	#define cudaTex1Dfetch(...) tex1Dfetch(__VA_ARGS__)
#endif

// max size for the 2D array is 65536 x 32768 (as of release 2.2), in bytes
template <class T> unsigned int getMaxTextureWidth(texture<T,2,cudaReadModeElementType> *tex) {
	if (sizeof(T) <= 4) return MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE;
	else if (sizeof(T) == 8) return MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE;
	else /*if (sizeof(T) == 16)*/ return MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE;
}

template <class T>
__device__ T cudaTexFetch(const texture<T,1,cudaReadModeElementType> tex, const unsigned int index1D) {
	return cudaTex1Dfetch(tex, index1D);
}
template <class T>
__device__ void cudaTexFetch(const texture<T,1,cudaReadModeElementType> tex, const unsigned int index1D, T &val1, T &val2) {
	val1 = cudaTex1Dfetch(tex, index1D);
	val2 = cudaTex1Dfetch(tex, index1D+1);
}

template <class T>
__device__ T cudaTexFetch(const texture<T,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_LG; // index1D / MAX_WIDTH
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_MINUS_1; //index1D % MAX_WIDTH;
	return tex2D(tex, x, y);
}

template <> __device__ float2 cudaTexFetch(const texture<float2,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_LG;
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_MINUS_1;
	return tex2D(tex, x, y);
}
template <> __device__ int2 cudaTexFetch(const texture<int2,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_LG;
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_MINUS_1;
	return tex2D(tex, x, y);
}
template <> __device__ uint2 cudaTexFetch(const texture<uint2,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_LG;
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_8_BYTE_VALUE_MINUS_1;
	return tex2D(tex, x, y);
}

template <> __device__ float4 cudaTexFetch(const texture<float4,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_LG;
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_MINUS_1;
	return tex2D(tex, x, y);
}
template <> __device__ int4 cudaTexFetch(const texture<int4,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_LG;
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_MINUS_1;
	return tex2D(tex, x, y);
}
template <> __device__ uint4 cudaTexFetch(const texture<uint4,2,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_LG;
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_16_BYTE_VALUE_MINUS_1;
	return tex2D(tex, x, y);
}

template <class T>
__device__ void cudaTexFetch(const texture<T,2,cudaReadModeElementType> tex, const unsigned int index1D, T &val1, T &val2) {
	unsigned int y = index1D >> MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_LG; // index1D / MAX_WIDTH
	unsigned int x = index1D & MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_MINUS_1; //index1D % MAX_WIDTH;

	val1 = tex2D(tex, x, y);
	//val2 = tex2D(tex, x+1, y);

	// to do: not sure if this will create a performance hit, but it is more correct as I think the x value could be clamped
	if (x != MAX_TEXTURE_2D_WIDTH_4_BYTE_VALUE_MINUS_1) val2 = tex2D(tex, x+1, y);
	else val2 = tex2D(tex, 0, y+1);
}

template <class T>
__device__ T cudaTexFetch(const texture<T,3,cudaReadModeElementType> tex, const unsigned int index1D) {
	// to do: there may be problems when accessing a texture where z != 0
	//        its either from the index or from loading to the texture
	//        I think this may be fixed since I found a possible bug in the loading part and it may be fixed
	//           the bug is related to allocating memory on the host using malloc instead of cudaMallocHost
	//
	//        but also very possible that its that max indices below are in bytes and not 4-byte blocks (see how 2D is done)

	// 2048 = max 3D texture dim
	unsigned int z = index1D >> 22; // assumes width*height = 2^22 = 2048*2048
	unsigned int y = index1D >> 11; // assumes width = 2^11 = 2048
	unsigned int x = index1D & (2047); //index1D % 2048;
	return tex3D(tex, x, y, z);
}
template <class T>
__device__ void cudaTexFetch(const texture<T,3,cudaReadModeElementType> tex, const unsigned int index1D, T &val1, T &val2) {
	// 2048 = max 3D texture dim
	unsigned int z = index1D >> 22; // assumes width*height = 2^22 = 2048*2048
	unsigned int y = index1D >> 11; // assumes width = 2^11 = 2048
	unsigned int x = index1D & (2047); //index1D % 2048;

	val1 = tex3D(tex, x, y, z);
	//val2 = tex3D(tex, x+1, y, z);

	// to do: not sure if this will create a performance hit, but it is more correct as I think the x value could be clamped
	if (x != 2047) val2 = tex3D(tex, x+1, y, z);
	else if (y != 2047) val2 = tex3D(tex, 0, y+1, z);
	else tex3D(tex, 0, 0, z+1);
}


// to do: not sure if this works correctly with how its being loaded in the load function below
__device__ float cudaTex1DfetchFrom3D4(const texture<float4,3,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int idx = index1D >> 4;
	unsigned int z = idx >> 22, y = idx >> 11, x = idx & (2047);
	
	// divide by an extra 4 to get the index we want
	/*unsigned int z = index1D >> 24, y = index1D >> 13, x = index1D & (511);*/
	switch (index1D&3) {  // mod 4
		case 0:  return tex3D(tex, x, y, z).x;
		case 1:  return tex3D(tex, x, y, z).y;
		case 2:  return tex3D(tex, x, y, z).z;
		default: return tex3D(tex, x, y, z).w;
	}
}

__device__ unsigned int cudaTex1DfetchFrom3D4(texture<uint4,3,cudaReadModeElementType> tex, const unsigned int index1D) {
	unsigned int z = index1D >> 24, y = index1D >> 13, x = index1D & (511);
	switch (index1D&3) {  // mod 4
		case 0:  return tex3D(tex, x, y, z).x;
		case 1:  return tex3D(tex, x, y, z).y;
		case 2:  return tex3D(tex, x, y, z).z;
		default: return tex3D(tex, x, y, z).w;
	}
}


template<class T> // faster version of mod:  i%n   ==   i&(n-1)
__device__ T modulo(const T leftOfMod, const T rightOfMod) { return leftOfMod & (rightOfMod-1); }

__device__ unsigned int getMostSignificantBit(const unsigned int x) { return 32-__clz(x); }

__device__ unsigned int getNumberOfTrailingZeros(const int x) {
	const int MultiplyDeBruijnBitPosition[32] = 
		{	0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
			31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9	};
	return MultiplyDeBruijnBitPosition[((x & -x) * 0x077CB531U) >> 27];
	//return __ffs(x);
}

__device__ int getBit(const int x, const int bitIdx) { if (x & 1<<bitIdx) return 1;  else return 0; }

__device__ unsigned int interleaveBits(const unsigned int x, const unsigned int y, const unsigned int z) {
	unsigned int val = z&512;  // 20 shifts, 29 or's, 30 and's
	val =  (val << 2) | ((z&256) | (x&512)) | ((y&512) << 1);
	val =  (val << 2) | ((z&128) | (x&256)) | ((y&256) << 1);
	val =  (val << 2) | ((z&64)  | (x&128)) | ((y&128) << 1);
	val =  (val << 2) | ((z&32)  | (x&64) ) | ((y&64)  << 1);
	val =  (val << 2) | ((z&16)  | (x&32) ) | ((y&32)  << 1);
	val =  (val << 2) | ((z&8)   | (x&16) ) | ((y&16)  << 1);
	val =  (val << 2) | ((z&4)   | (x&8)  ) | ((y&8)   << 1);
	val =  (val << 2) | ((z&2)   | (x&4)  ) | ((y&4)   << 1);
	val =  (val << 2) | ((z&1)   | (x&2)  ) | ((y&2)   << 1);
	return (val << 2) |            (x&1)    | ((y&1)   << 1);
}

__device__ unsigned int interleaveBits(const unsigned int x, const unsigned int y) {
	unsigned int val = y&32768;  // 16 shifts, 32 or's, 32 and's
	val =  (val << 1) | (x&32768) | (y&16384);
	val =  (val << 1) | (x&16384) | (y&8192);
	val =  (val << 1) | (x&8192)  | (y&4096);
	val =  (val << 1) | (x&4096)  | (y&2048);
	val =  (val << 1) | (x&2048)  | (y&1024);
	val =  (val << 1) | (x&1024)  | (y&512);
	val =  (val << 1) | (x&512)   | (y&256);
	val =  (val << 1) | (x&256)   | (y&128);
	val =  (val << 1) | (x&128)   | (y&64);
	val =  (val << 1) | (x&64)    | (y&32);
	val =  (val << 1) | (x&32)    | (y&16);
	val =  (val << 1) | (x&16)    | (y&8);
	val =  (val << 1) | (x&8)     | (y&4);
	val =  (val << 1) | (x&4)     | (y&2);
	val =  (val << 1) | (x&2)     | (y&1);
	return (val << 1) | (x&1);
}

__device__ unsigned int rgbaFloatToInt(float4 rgba) {
	rgba.x = __saturatef(rgba.x);  // clamp to [0.0, 1.0]
	rgba.y = __saturatef(rgba.y);
	rgba.z = __saturatef(rgba.z);
	rgba.w = __saturatef(rgba.w);
	return (unsigned int(rgba.w*255)<<24) | (unsigned int(rgba.z*255)<<16) | (unsigned int(rgba.y*255)<<8) | unsigned int(rgba.x*255);
}

__device__ int isPointInBox(const float3 pt, const float3 boxMin, const float3 boxMax) {
	if (pt.x >= boxMin.x && pt.y >= boxMin.y && pt.z >= boxMin.z &&
		pt.x <= boxMax.x && pt.y <= boxMax.y && pt.z <= boxMax.z) return 1;
	else return 0;
}

template <class T> __device__ T bilinearLerp(const T val00, const T val01, const T val10, const T val11, const float tx, const float ty) {
	float ty2 = 1.0-ty;
	T i1 = ty*val00 + ty2*val01;
	T i2 = ty*val10 + ty2*val11;
	return tx*i1 + (1.0-tx)*i2;
}

template <class T>
__device__ T trilinearLerp(const T val000, const T val001, const T val010, const T val011, 
						   const T val100, const T val101, const T val110, const T val111,
						   const float tx, const float ty, const float tz) {
	float tz2 = 1.0-tz;	
	T i1 = tz*val000 + tz2*val001;
	T i2 = tz*val100 + tz2*val101;
	T j1 = tz*val010 + tz2*val011;
	T j2 = tz*val110 + tz2*val111;

	float ty2 = 1.0-ty;
	T w1 = ty*i1 + ty2*j1;
	T w2 = ty*i2 + ty2*j2;

	return tx*w1 + (1.0-tx)*w2;
}

__device__ void getBoundingIntegers(const float val, int &lower, int &upper, const float tolerance = 1e-6) {
	lower = floorf(val);
	if (fabs(lower-val) < tolerance) upper = lower; // val = x.00001 or -x.99999
	else {
		upper = lower+1;
		if (fabs(upper-val) < tolerance) lower = upper; // val = x.99999 or -x.00001
	}
}

// returns the index if it is found, or returns -1 if not
// searches between [lowIdx, highIdxPlus1-1], hence the plus 1 part
template <class T>
__device__ int binarySearch_sortedLowToHigh(T *vec, const T target, const unsigned int lowIdx, const unsigned int highIdxPlus1) {
	unsigned int low = lowIdx, high = highIdxPlus1;
	while (low < high) {
		unsigned int mid = (low+high)/2;

		// delays the == comparison until after the while loop
		if (vec[mid] < target) low = mid+1;
		else high = mid;
	}
	// high == low here
	if (low < highIdxPlus1 && vec[low] == target) return low;
	else return -1;
}

template <class T, int N>
__device__ int binarySearch_sortedLowToHigh(texture<T,N,cudaReadModeElementType> pseudo1DTex, const T target, const unsigned int lowIdx, const unsigned int highIdxPlus1) {
	unsigned int low = lowIdx, high = highIdxPlus1;

	// note: it seems that unrolling makes this slower
	//       and the only way to do the unroll with "#pragma unroll" is to keep track of the iterations of the while loop
	//         else the compiler (CUDA v2.3) will give a warning that the unroll is ignored
	while (low < high) {
		unsigned int mid = (low+high)>>1;
		if (cudaTexFetch(pseudo1DTex, mid) < target) low = mid+1;
		else high = mid;
	}
	if (low < highIdxPlus1 && cudaTexFetch(pseudo1DTex, low) == target) return low;
	else return -1;
}

template <class T, int N>
__device__ unsigned int binarySearch_sortedLowToHigh_guaranteedToBeFound(texture<T,N,cudaReadModeElementType> pseudo1DTex, const T target, const unsigned int lowIdx, const unsigned int highIdxPlus1) {
	unsigned int low = lowIdx, high = highIdxPlus1;

	while (low < high) {
		unsigned int mid = (low+high)>>1;
		if (cudaTexFetch(pseudo1DTex, mid) < target) low = mid+1;
		else high = mid;
	}
	return low;
}



template <class T>
void loadVolumeToTexture(T *sourceVolumeData, cudaExtent volumeSize, cudaArray **outputArray, texture<T,3,cudaReadModeElementType> *tex) {
	// create 3D array
	cudaChannelFormatDesc channelDesc = tex->channelDesc;
	cutilSafeCall(cudaMalloc3DArray(outputArray, &channelDesc, volumeSize));

	// copy data to 3D array
	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr   = make_cudaPitchedPtr((void*)sourceVolumeData, volumeSize.width*sizeof(T), volumeSize.width, volumeSize.height);
	copyParams.dstArray = *outputArray;
	copyParams.extent   = volumeSize;
	copyParams.kind     = cudaMemcpyHostToDevice;
	cutilSafeCall(cudaMemcpy3D(&copyParams));

	// bind array to 3D texture
	cutilSafeCall(cudaBindTextureToArray(*tex, *outputArray, channelDesc));
}

template <class T>
void loadVectorToArray(T *sourceVec, unsigned int numberEntries, cudaArray **destArray) {
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
	cutilSafeCall(cudaMallocArray(destArray, &channelDesc, numberEntries));
	cutilSafeCall(cudaMemcpyToArray(*destArray, 0, 0, (void*)sourceVec, numberEntries*sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void loadVectorToTexture(T *sourceVec, unsigned int numberEntries, cudaArray **destArray, texture<T,1,cudaReadModeElementType> *tex) {
	if (numberEntries >= 65536) printf("the max size for the 1D array/texture is 64k (as of release 2.2)\n");

	cudaChannelFormatDesc channelDesc = tex->channelDesc;
	cutilSafeCall(cudaMallocArray(destArray, &channelDesc, numberEntries));
	cutilSafeCall(cudaMemcpyToArray(*destArray, 0, 0, (void*)sourceVec, numberEntries*sizeof(T), cudaMemcpyHostToDevice));
	cutilSafeCall(cudaBindTextureToArray(*tex, destArray, channelDesc));
}


template <class T>
void loadAreaToTexture(T *sourceAreaData, unsigned int width, unsigned int height, cudaArray **outputArray,
					   texture<T,2,cudaReadModeElementType> *tex, cudaMemcpyKind memcpyType = cudaMemcpyHostToDevice) {
	cudaChannelFormatDesc channelDesc = tex->channelDesc;
	cutilSafeCall(cudaMallocArray(outputArray, &channelDesc, width, height));
	cutilSafeCall(cudaMemcpyToArray(*outputArray, 0, 0, (void*)sourceAreaData, width*height*sizeof(T), memcpyType));

	cutilSafeCall(cudaBindTextureToArray(*tex, *outputArray, channelDesc));
}

template <class T>
void loadVectorToTexture2D(T *sourceVec, const unsigned int numberEntries, cudaArray **destArray, texture<T,2,cudaReadModeElementType> *tex) {
	unsigned int width = getMaxTextureWidth(tex);
	unsigned int height = 1;
	if (width >= numberEntries) width = numberEntries;
	else height = (unsigned int)ceil((double)numberEntries/(double)width);
	
	printf("%d <= %d x %d = %d\n", numberEntries, width, height, width*height);
	
	size_t totalSize = width*height;
	if (totalSize > (size_t)numberEntries) {
		T *biggerSourceVec;
		cutilSafeCall(cudaMallocHost((void**)&biggerSourceVec, totalSize*sizeof(T)));

		// don't care about values in range [numberEntires, totalSize) because in theory they should not be accessed
		for (unsigned int i=0; i<numberEntries; i++) biggerSourceVec[i] = sourceVec[i];

		loadAreaToTexture(biggerSourceVec, width, height, destArray, tex);

		cutilSafeCall(cudaFreeHost((void*)biggerSourceVec));
	}
	else loadAreaToTexture(sourceVec, width, height, destArray, tex);
}


// don't really want to template this, but for whatever reason the compiler gives a link error
// "fatal error LNK1169: one or more multiply defined symbols found"
template <class T>
cudaExtent computeBrickSize(const T numberEntries) {
	// for now, just allocated blocks of say 2k x 2k x N until the vector fits in memory
	// to do: find a block size where the first 2 dimensions are powers of 2 and there is the least amount of wasted space
	
	unsigned int width = 2048;
	unsigned int height = 1, depth = 1;
	if (width >= numberEntries) width = numberEntries;
	else {
		if ((unsigned int)ceil((double)numberEntries/(double)width) > 2048) {
			height = 2048;
			depth = (unsigned int)ceil((double)numberEntries/(double)(width*height));
		}
		else height = (unsigned int)ceil((double)numberEntries/(double)width);
	}
	
	printf("%d <= %d x %d x %d = %d\n", numberEntries, width, height, depth, width*height*depth);

	return make_cudaExtent(width, height, depth);
}

// max for 3D array = 2048*2048*2048 (as of release 2.2)
template <class T>
void loadVectorToTexture3D(T *sourceVec, const unsigned int numberEntries, cudaArray **destArray, texture<T,3,cudaReadModeElementType> *tex) {
	cudaExtent volumeSize = computeBrickSize(numberEntries);

	// to do: is there a way to copy over the information without needing to allocate more memory?
	//           or at least try to limit it be doing a slice by slice copy using an offset parameter
	size_t totalSize = volumeSize.width*volumeSize.height*volumeSize.depth;
	if (totalSize > (size_t)numberEntries) {
		T *biggerSourceVec;
		cutilSafeCall(cudaMallocHost((void**)&biggerSourceVec, totalSize*sizeof(T)));

		// don't care about values in range [numberEntires, totalSize) because in theory they should not be accessed
		for (unsigned int i=0; i<numberEntries; i++) biggerSourceVec[i] = sourceVec[i];

		loadVolumeToTexture(biggerSourceVec, volumeSize, destArray, tex);

		cutilSafeCall(cudaFreeHost((void*)biggerSourceVec));
	}
	else loadVolumeToTexture(sourceVec, volumeSize, destArray, tex);
}

template <class T, class T4>
void loadVectorToTexture3D4(T *sourceVec, const unsigned int numberEntries, cudaArray **destArray, texture<T4,3,cudaReadModeElementType> *tex) {
	unsigned int actualSize = (unsigned int)ceil(numberEntries*0.25f);
	cudaExtent volumeSize = computeBrickSize(actualSize);
	size_t totalSize = volumeSize.width*volumeSize.height*volumeSize.depth;

	T4 *biggerSourceVec = (T4*)malloc(totalSize*sizeof(T4));

	unsigned int i, count=0;
	for (i=0; i<numberEntries; i+=4) {
		biggerSourceVec[count].x = sourceVec[ i ];
		biggerSourceVec[count].y = sourceVec[i+1];
		biggerSourceVec[count].z = sourceVec[i+2];
		biggerSourceVec[count].w = sourceVec[i+3];
		count++;
	}
	if (i<numberEntries) { biggerSourceVec[count].x = sourceVec[i];  i++; }
	if (i<numberEntries) { biggerSourceVec[count].y = sourceVec[i];  i++; }
	if (i<numberEntries) { biggerSourceVec[count].z = sourceVec[i];  i++; }

	loadVolumeToTexture(biggerSourceVec, volumeSize, destArray, tex);

	free(biggerSourceVec);  biggerSourceVec = NULL;
}

#endif