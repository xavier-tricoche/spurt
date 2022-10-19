#ifndef CUDA_WRAPPER_H
#define CUDA_WRAPPER_H

#include "CudaConstants.cuh"

#ifdef COMPILE_WITH_CUDA

#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <cutil_inline.h>
#include <cutil_gl_inline.h>
#include <cuda_gl_interop.h>

#pragma comment(lib, "cudart.lib")

#if defined(_WIN64)  // note: VS syntax highlighting may grey out below even though it is defined and will be used
	#if defined(_DEBUG)
		#pragma comment(lib, "cutil64D.lib")
		#pragma comment(lib, "cudpp64d.lib")
	#else
		#pragma comment(lib, "cutil64.lib")
		#pragma comment(lib, "cudpp64.lib")
	#endif
#elif defined(_WIN32)
	#if defined(_DEBUG)
		#pragma comment(lib, "cutil32D.lib")
		#pragma comment(lib, "cudpp32d.lib")
	#else
		#pragma comment(lib, "cutil32.lib")
		#pragma comment(lib, "cudpp32.lib")
	#endif
#endif

namespace CudaWrapper {
	#include <GL/glew.h>

	static inline void init() { cudaGLSetGLDevice(cutGetMaxGflopsDeviceId()); }

	static inline void removeBufferObject(GLuint *bufferObject) {
		if (bufferObject != NULL && *bufferObject != 0) {
			cutilSafeCall(cudaGLUnregisterBufferObject(*bufferObject));
			glDeleteBuffersARB(1, bufferObject);
		}
	}

	static inline void addPixelBufferObject(GLuint *bufferObject, const unsigned int bufferSizeInBytes) {
		glGenBuffersARB(1, bufferObject);
		glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, *bufferObject);
		glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, bufferSizeInBytes, 0, GL_STREAM_DRAW_ARB);
		glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

		cutilSafeCall(cudaGLRegisterBufferObject(*bufferObject));
	}

	static inline void addVertexBufferObject(GLuint *bufferObject, const unsigned int bufferSizeInBytes) {
		glGenBuffersARB(1, bufferObject);
		glBindBuffer(GL_ARRAY_BUFFER, *bufferObject);
		glBufferData(GL_ARRAY_BUFFER, bufferSizeInBytes, 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		cutilSafeCall(cudaGLRegisterBufferObject(*bufferObject));
	}

	template <class T> static inline void mapBufferObjectToArray(const GLuint bufferObject, T **mappedArray)
		{ cutilSafeCall(cudaGLMapBufferObject((void**)mappedArray, bufferObject)); }

	template <class T> static inline void mapBufferObjectToArray(const GLuint bufferObject, T **mappedArray, const unsigned int numberOfArrayEntries, const int resetValue) {
		mapBufferObjectToArray(bufferObject, mappedArray);
		cutilSafeCall(cudaMemset(*mappedArray, resetValue, numberOfArrayEntries));
	}

	static inline void unmapBufferObject(const GLuint bufferObject) { cutilSafeCall(cudaGLUnmapBufferObject(bufferObject)); }

	template <class T> static inline void cudaNew(T **val, const size_t numberOfElements) { cutilSafeCall(cudaMallocHost((void**)val, numberOfElements*sizeof(T))); }
	template <class T> static inline void cudaDelete(T *val) { cutilSafeCall(cudaFreeHost((void*)val));  val = NULL; }

	static inline void convert(const Vector2f input, float2 &output) { output = make_float2(input.x, input.y); }
	static inline void convert(const Vector2i input, int2 &output) { output = make_int2(input.x, input.y); }
	static inline void convert(const Vector2ui input, uint2 &output) { output = make_uint2(input.x, input.y); }

	static inline void convert(const Vector4f input, float4 &output) { output = make_float4(input.x, input.y, input.z, input.w); }
	static inline void convert(const Vector4i input, int4 &output) { output = make_int4(input.x, input.y, input.z, input.w); }
	static inline void convert(const Vector4ui input, uint4 &output) { output = make_uint4(input.x, input.y, input.z, input.w); }

	static inline void convert(const Vector3f input, float4 &output) { output = make_float4(input.x, input.y, input.z, 0); }
	static inline void convert(const Vector3i input, int4 &output) { output = make_int4(input.x, input.y, input.z, 0); }
	static inline void convert(const Vector3ui input, uint4 &output) { output = make_uint4(input.x, input.y, input.z, 0); }

	static inline void convert(const Vector3f input, float3 &output) { output = make_float3(input.x, input.y, input.z); }
	static inline void convert(const Vector3i input, int3 &output) { output = make_int3(input.x, input.y, input.z); }
	static inline void convert(const Vector3ui input, uint3 &output) { output = make_uint3(input.x, input.y, input.z); }


	template <class T> static inline void convert(const T *input, const size_t size, T **output) {
		cudaNew(output, size);
		for (unsigned int i=0; i<size; i++) (*output)[i] = input[i];
	}
	
	template <class T> static inline void convert(const Vector2<T> *input, const size_t size, T **output) {
		cudaNew(output, size*2);
		for (unsigned int i=0; i<size; i++) {
			(*output)[ i*2 ] = input[i].x;
			(*output)[i*2+1] = input[i].y;
		}
	}

	template <class T> static inline void convert(const Vector4<T> *input, const size_t size, T **output) {
		cudaNew(output, size*4);
		for (unsigned int i=0; i<size; i++) {
			(*output)[ i*4 ] = input[i].x;
			(*output)[i*4+1] = input[i].y;
			(*output)[i*4+2] = input[i].z;
			(*output)[i*4+3] = input[i].w;
		}
	}

	template <class T> static inline void convert(const Vector3<T> *input, const size_t size, T **output, bool alignTo4) {
		unsigned int dim = alignTo4 ? 4 : 3;
		cudaNew(output, size*dim);

		for (unsigned int i=0; i<size; i++) {
			(*output)[ i*dim ] = input[i].x;
			(*output)[i*dim+1] = input[i].y;
			(*output)[i*dim+2] = input[i].z;
			if (alignTo4) (*output)[i*dim+3] = 0;
		}
	}

	template <class T, class S> static inline void convert(const vector<T> &input, S **output) { convert(&input[0], input.size(), output); }
	template <class T> static inline void convert(const vector< Vector3<T> > &input, T **output, bool alignTo4) { convert(&input[0], input.size(), output, alignTo4); }

	static inline void convert(const vector< Vector3f > &input, float3 **output) {
		cudaNew(output, input.size()*3);
		for (unsigned int i=0; i<input.size(); i++) { convert(input[i], (*output)[i]); }
	}
}

#endif

#endif