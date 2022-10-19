#ifndef CUDA_TRANSFER_FUNCTION_CUH
#define CUDA_TRANSFER_FUNCTION_CUH

#define MAX_TRANSFER_FUNCTIONS 3

// unfortunately have to list out the texture instead of doing an array right now...
texture<float4, 1, cudaReadModeElementType> transferTex0, transferTex1, transferTex2; // 1D transfer function texture
__constant__ float minimumDataValue[MAX_TRANSFER_FUNCTIONS], maximumDataValue[MAX_TRANSFER_FUNCTIONS], invDataValueDim[MAX_TRANSFER_FUNCTIONS];

#define transferFunctionArrayFetch(_texbasename, _texnum)\
switch(_texnum) {\
	case 0: return _texbasename##0;\
	case 1: return _texbasename##1;\
	default: return _texbasename##2;\
}

#define transferFunctionArrayPtrFetch(_texbasename, _texnum)\
switch(_texnum) {\
	case 0: return &_texbasename##0;\
	case 1: return &_texbasename##1;\
	default: return &_texbasename##2;\
}

__device__ __host__ texture<float4, 1, cudaReadModeElementType> transferFunctionGetTex(const unsigned int id) { transferFunctionArrayFetch(transferTex, id) }
texture<float4, 1, cudaReadModeElementType>* transferFunctionGetTexPtr(const unsigned int id) { transferFunctionArrayPtrFetch(transferTex, id) }

extern "C"
void cudaLoadTransferFunction(float *transferFunction, const unsigned int transferFunctionDim, const float minimumValue,
							  const float maximumValue, const unsigned int transferFunctionID) {
	if (transferFunctionID >= MAX_TRANSFER_FUNCTIONS) {
		printf("Bad ID passed to cudaLoadTransferFunction, may need to increase MAX_TRANSFER_FUNCTIONS\n");
		return;
	}

	// create transfer function texture
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	cudaArray *transferFuncArray;
	cutilSafeCall(cudaMallocArray(&transferFuncArray, &channelDesc, transferFunctionDim, 1)); 
	cutilSafeCall(cudaMemcpyToArray(transferFuncArray, 0, 0, transferFunction, transferFunctionDim*4*sizeof(float), cudaMemcpyHostToDevice));

	texture<float4, 1, cudaReadModeElementType> *transferTexPtr;

	transferTexPtr = transferFunctionGetTexPtr(transferFunctionID);

	//transferTexPtr->filterMode = cudaFilterModeLinear;
	transferTexPtr->normalized = true;    // access with normalized texture coordinates
	transferTexPtr->addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

	// bind the array to the texture
	cutilSafeCall(cudaBindTextureToArray(*transferTexPtr, transferFuncArray, channelDesc));
	transferTexPtr->filterMode = cudaFilterModeLinear;


	//float invDim = (float)1.0/(maximumValue-minimumValue);
	// the factor needs to be a little bigger than 1.0 to match CPU rendering
	//   below is still not correct though
	float invDim = (float)((float)transferFunctionDim / (float)(transferFunctionDim-3.0)) /(maximumValue-minimumValue);

	cutilSafeCall(cudaMemcpyToSymbol(minimumDataValue, &minimumValue, sizeof(float), transferFunctionID*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(maximumDataValue, &maximumValue, sizeof(float), transferFunctionID*sizeof(float)));
	cutilSafeCall(cudaMemcpyToSymbol(invDataValueDim, &invDim, sizeof(float), transferFunctionID*sizeof(float)));
}

__device__ float4 transferFunctionGetValue(const float val, const unsigned int textureID) { // val = [0,1]
	return tex1D(transferFunctionGetTex(textureID), val);
}

__device__ float transferFunctionValueNormalization(const float val, const unsigned int transferFuncID) { // normalizes to [0,1]
	return (val-minimumDataValue[transferFuncID])*invDataValueDim[transferFuncID];
}

#endif