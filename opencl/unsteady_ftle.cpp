/*
    FTLE computation

    This program loads a file describing a flow field and performs an integration 
	for the computation of FTLE.
*/

// OpenGL Graphics Includes
#include <GL/glew.h>
#if defined(__APPLE__) || defined(__MACOSX)
    #include <OpenGL/OpenGL.h>
    #include <GLUT/glut.h>
#else
    #include <GL/freeglut.h>
    #ifdef UNIX
       #include <GL/glx.h>
    #endif
#endif

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <simd/vector.h>
//#include <vector_functions.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <teem/nrrd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

// Utilities, OpenCL and system includes
// #include <oclUtils.h>
// #include <shrQATest.h>

#if defined (__APPLE__) || defined(MACOSX)
   #define GL_SHARING_EXTENSION "cl_APPLE_gl_sharing"
#else
   #define GL_SHARING_EXTENSION "cl_khr_gl_sharing"
#endif

using namespace std;

// Constants, defines, typedefs and global declarations
//*****************************************************************************

int *pArgc = NULL;
char **pArgv = NULL;

typedef unsigned int uint;
typedef unsigned char uchar;

uint width, height, depth;
char directory[2048];
uint steps;
float* flow_stamps;
char** flow_files;
float4 minB;
float4 maxB;
float4 spacing;
float4* h_flowMap;

// OpenCL vars
cl_platform_id cpPlatform;
cl_device_id* cdDevices;
cl_uint uiDeviceUsed;
cl_uint uiDevCount;
cl_context cxGPUContext;
cl_command_queue cqCommandQueue;
cl_program cpProgram;
cl_kernel ckKernel;
cl_int ciErrNum;
char* cPathAndName = NULL;          // var for full paths to data, src, etc.
char* cSourceCL;                    // Buffer to hold source for compilation 
const char* cExecutableName = NULL;
cl_bool g_bImageSupport;
cl_mem d_flowMap;
cl_mem d_volumeFlow0;
cl_mem d_volumeFlow1;
cl_sampler volumeSampler0;
cl_sampler volumeSampler1;

// Helpers
void Cleanup(int iExitCode);
void (*pCleanup)(int) = &Cleanup;

// General utility functions
//*****************************************************************************
Nrrd *readNrrd(const char* filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename, NULL)) {
        printf("Error reading %s\n",filename);
                return NULL;
    }
    return nin;
}

void writeNrrd(void* data, const string& filename, int data_type, const vector< size_t >& dims, const vector< double >& spacing)
{
        Nrrd *nout = nrrdNew();

        if (nrrdWrap_nva(nout, data, data_type, dims.size(), &dims[0])) {
                cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
                exit(-1);
        }
        if (spacing.size() == dims.size()) {
                nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&spacing[0]);
        }
        if (nrrdSave(filename.c_str(), nout, NULL)) {
                cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
                exit(-1);
        }
}

void writeRawFile(float* data,const char *filename, int width, int height, int depth)
{
        vector< double > spacing;
        spacing.push_back(1.0);
        spacing.push_back(1.0);
        spacing.push_back(1.0);

        vector< size_t > dims;
        dims.push_back(width);
        dims.push_back(height);
        dims.push_back(depth);

        string file_string(filename);
        writeNrrd(data, file_string, nrrdTypeFloat, dims, spacing);

        printf("Write '%s'\n", filename);
}

void writeRawFile(float* data,const char *filename, int vec, int width, int height, int depth)
{
        vector< double > spacing;
        spacing.push_back(1.0);
        spacing.push_back(1.0);
        spacing.push_back(1.0);
        spacing.push_back(1.0);

        vector< size_t > dims;
        dims.push_back(vec);
        dims.push_back(width);
        dims.push_back(height);
        dims.push_back(depth);

        string file_string(filename);
        writeNrrd(data, file_string, nrrdTypeFloat, dims, spacing);

        printf("Write '%s'\n", filename);
}

// Functions for volume access
//*****************************************************************************
int cint(double x)
{
	double fractpart, intpart;
	fractpart = abs(modf (x , &intpart));

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

bool IsBoundary(const uint& x, const uint& y, const uint& z)
{
	if ((x < 1) || (y < 1) || (z < 1) || (x > width - 2) || (y > height - 2) || (z > depth - 2))
		return true;
	else return false;
}

bool IsValid(const int& x, const int& y, const int& z)
{
	if ((x < 0) || (y < 0) || (z < 0) || ((uint)x > width - 1) || ((uint)y > height - 1) || ((uint)z > depth - 1))
		return false;
	else return true;
}

bool InDomain(const float4& spt)
{
	if ((spt.x < minB.x) || (spt.x > maxB.x)
	 || (spt.y < minB.y) || (spt.y > maxB.y)
	 || (spt.z < minB.z) || (spt.z > maxB.z))
		return false;
	return true;
}

int Coord2Addr(const uint& x, const uint& y, const uint& z)
{
	return x + width * y + width * height * z;
}

int Coord2Addr(const int3& coord)
{
	return Coord2Addr(coord.x, coord.y, coord.z);
}

int3 Addr2Coord(const uint& i)
{
	uint x = i % width;
	uint y = ((i - x) / width) % height;
	uint z = (((i - x) / width) - y) / height;

	return make_int3(x, y, z);
}

float4 Grid2Space(const float4& gpt)
{
	float4 spt;

	spt.x = minB.x + (maxB.x - minB.x) * gpt.x / (width - 1.0);
	spt.y = minB.y + (maxB.y - minB.y) * gpt.y / (height - 1.0);
	spt.z = minB.z + (maxB.z - minB.z) * gpt.z / (depth - 1.0);

	return spt;
}

int3 Space2Grid(const float4& spt)
{
	int3 index;
	index.x = cint((width - 1.0) * (spt.x - minB.x) / (maxB.x - minB.x));
	index.y = cint((height - 1.0) * (spt.y - minB.y) / (maxB.y - minB.y));
	index.z = cint((depth - 1.0) * (spt.z - minB.z) / (maxB.z - minB.z));
	return index;
}

// Intitialize OpenCL
//*****************************************************************************
void createCLContext(int argc, const char** argv) {
	// For cards information: less /proc/driver/nvidia/gpus/0/information
	
    //Get the NVIDIA platform
    ciErrNum = oclGetPlatformID(&cpPlatform);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	printf("Platform id: %d\n", cpPlatform);

    // Get the number of GPU devices available to the platform
    ciErrNum = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 0, NULL, &uiDevCount);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	printf("Found %d OpenCL GPU devices.\n", uiDevCount);

    // Create the device list
    cdDevices = new cl_device_id [uiDevCount];
    ciErrNum = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, uiDevCount, cdDevices, NULL);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);

	// Print all devices info
	for (int i = 0; i < uiDevCount; i++)
	{
		//oclPrintDevInfo(LOGBOTH, cdDevices[i]);
	}
	
    // Get device requested on command line, if any
    uiDeviceUsed = 0;
    unsigned int uiEndDev = uiDevCount - 1;
    if(shrGetCmdLineArgumentu(argc, argv, "device", &uiDeviceUsed ))
    {
      uiDeviceUsed = CLAMP(uiDeviceUsed, 0, uiEndDev);
      uiEndDev = uiDeviceUsed; 
    } 

	// Check if the requested device (or any of the devices if none requested) supports context sharing with OpenGL
	// No GL interop
	cl_context_properties props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cpPlatform, 0};
	cxGPUContext = clCreateContext(props, 1, &cdDevices[uiDeviceUsed], NULL, NULL, &ciErrNum);


    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
}

// Load file for the time steps
//*****************************************************************************
void loadStepsFile(char* filename)
{
	printf("\n\n\nReading dataset information from file %s\n", filename);
	
	ifstream fp_in;  // declarations of streams fp_in and fp_out
	fp_in.open(filename, ios::in);    // open the streams
	
	fp_in >> width >> height >> depth;
	cout << "Dimensions: " << width << " " << height << " " << depth << "\n";
	fp_in >> directory;
	cout << "Directory: " << directory << "\n";
	fp_in >> steps;
	cout << "Number of time steps: " << steps << "\n";
	fp_in >> minB.x >> minB.y >> minB.z;
	cout << "BBox min: " << minB.x << " " << minB.y << " " << minB.z << "\n";
	fp_in >> maxB.x >> maxB.y >> maxB.z;
	cout << "BBox max: " << maxB.x << " " << maxB.y << " " << maxB.z << "\n";
	fp_in >> spacing.x >> spacing.y >> spacing.z;
	cout << "Spacing: " << spacing.x << " " << spacing.y << " " << spacing.z << "\n";
	minB.w = 0.0;
	maxB.w = 0.0;
	
	flow_stamps = (float*) malloc(steps * sizeof(float));
	flow_files = (char**) malloc(steps * sizeof(char*));
	
	for (uint i = 0; i < steps; i++)
	{
		printf("%d ", i);
		flow_files[i] = (char*) malloc(1024 * sizeof(char));
		fp_in >> flow_files[i];
		fp_in >> flow_stamps[i];
	}
	
	fp_in.close();   // close the streams
}

// Create textures and flow map
//*****************************************************************************
void initCLGPUMemory()
{
	ciErrNum = CL_SUCCESS;
	
	// create GPU flow map
	d_flowMap = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE, width * height * depth * sizeof(float4), NULL, &ciErrNum);
	ciErrNum |= clSetKernelArg(ckKernel, 0, sizeof(cl_mem), (void *) &d_flowMap);
	oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	
	// create CPU flow map
	h_flowMap = (float4*) malloc(width * height * depth * sizeof(float4));
	for (uint i = 0; i < width * height * depth; i++)
	{
		int3 g = Addr2Coord(i);
		float4 gpt = (float4)(g.x, g.y, g.z, 0.0);
		float4 spt = Grid2Space(gpt);
		h_flowMap[i] = (float4)(spt.x, spt.y, spt.z, 0.0);
	}
	
	// copy CPU to GPU
	ciErrNum = clEnqueueWriteBuffer(cqCommandQueue, d_flowMap, CL_FALSE, 0, width * height * depth * sizeof(float4), h_flowMap, 0, NULL, NULL);
	oclCheckError(ciErrNum, CL_SUCCESS);
	ciErrNum = clFinish(cqCommandQueue);
	oclCheckError(ciErrNum, CL_SUCCESS);
}

// Copy flow map from GPU to CPU
//*****************************************************************************
void copyGPUFlowMap()
{
	// copy GPU to CPU
	ciErrNum = clEnqueueReadBuffer(cqCommandQueue, d_flowMap, CL_FALSE, 0, width * height * depth * sizeof(float4), h_flowMap, 0, NULL, NULL);
	oclCheckError(ciErrNum, CL_SUCCESS);
	ciErrNum = clFinish(cqCommandQueue);
	oclCheckError(ciErrNum, CL_SUCCESS);
}

// Load time steps in textures
//*****************************************************************************
void loadFlowTextures(char* filename0, char* filename1)
{
	// free original textures
	if(volumeSampler0)clReleaseSampler(volumeSampler0);
	if(volumeSampler1)clReleaseSampler(volumeSampler1); 
	if(d_volumeFlow0)clReleaseMemObject(d_volumeFlow0);   
	if(d_volumeFlow1)clReleaseMemObject(d_volumeFlow1);   
	
	// read nrrd files
	Nrrd* nrrdflow_0 = readNrrd((string(directory) + string("/") + string(filename0)).c_str());
	Nrrd* nrrdflow_1 = readNrrd((string(directory) + string("/") + string(filename1)).c_str());
	
	// texture format
	cl_image_format volume_format;
	volume_format.image_channel_order = CL_RGBA;
	volume_format.image_channel_data_type = CL_FLOAT;
	
	// create textures
	float4* h_volumeFlow = (float4*)malloc(width * height * depth * sizeof(float4));
	
	for (uint i = 0; i < width * height * depth; i++)
	{
		h_volumeFlow[i].x = ((float*)nrrdflow_0->data)[3 * i + 0];
		h_volumeFlow[i].y = ((float*)nrrdflow_0->data)[3 * i + 1];
		h_volumeFlow[i].z = ((float*)nrrdflow_0->data)[3 * i + 2];
		h_volumeFlow[i].w = 0.0;
	}
	d_volumeFlow0 = clCreateImage3D(cxGPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, &volume_format, 
                                    width, height, depth,
                                    (width * sizeof(float4)), (width * height * sizeof(float4)),
                                    h_volumeFlow, &ciErrNum);
	oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	
	for (uint i = 0; i < width * height * depth; i++)
	{
		h_volumeFlow[i].x = ((float*)nrrdflow_1->data)[3 * i + 0];
		h_volumeFlow[i].y = ((float*)nrrdflow_1->data)[3 * i + 1];
		h_volumeFlow[i].z = ((float*)nrrdflow_1->data)[3 * i + 2];
		h_volumeFlow[i].w = 0.0;
	}
	d_volumeFlow1 = clCreateImage3D(cxGPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, &volume_format, 
                                    width, height, depth,
                                    (width * sizeof(float4)), (width * height * sizeof(float4)),
                                    h_volumeFlow, &ciErrNum);
	oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	
	free(h_volumeFlow);
	
	// create the sample objects
	volumeSampler0 = clCreateSampler(cxGPUContext, false, CL_ADDRESS_REPEAT, CL_FILTER_LINEAR, &ciErrNum);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	volumeSampler1 = clCreateSampler(cxGPUContext, false, CL_ADDRESS_REPEAT, CL_FILTER_LINEAR, &ciErrNum);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	
	// set parameters
	ciErrNum |= clSetKernelArg(ckKernel, 8, sizeof(cl_mem), (void *) &d_volumeFlow0);
	ciErrNum |= clSetKernelArg(ckKernel, 11, sizeof(cl_mem), (void *) &d_volumeFlow1);
	ciErrNum |= clSetKernelArg(ckKernel, 9, sizeof(cl_sampler), &volumeSampler0);
	ciErrNum |= clSetKernelArg(ckKernel, 12, sizeof(cl_sampler), &volumeSampler1);
	oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
	
	// free memory
	nrrdNuke(nrrdflow_0);
	nrrdNuke(nrrdflow_1);
}

// Invoke the kernel
//*****************************************************************************
float4 advcminB, advcmaxB, advspacing;
void advect(float st, float et, uint start_step)
{
	ciErrNum = CL_SUCCESS;
	
	// set parameters
	ciErrNum |= clSetKernelArg(ckKernel, 0, sizeof(cl_mem), (void *) &d_flowMap);
	ciErrNum |= clSetKernelArg(ckKernel, 1, sizeof(uint), &width);
	ciErrNum |= clSetKernelArg(ckKernel, 2, sizeof(uint), &height);
	ciErrNum |= clSetKernelArg(ckKernel, 3, sizeof(uint), &depth);
	
	advcminB = float4(minB.x, minB.y, minB.z, 0.0);
	advcmaxB = float4(maxB.x, maxB.y, maxB.z, 0.0);
	advspacing = float4(spacing.x, spacing.y, spacing.z, 0.0);
	ciErrNum |= clSetKernelArg(ckKernel, 4, sizeof(float4), &advcminB);
	ciErrNum |= clSetKernelArg(ckKernel, 5, sizeof(float4), &advcmaxB);
	ciErrNum |= clSetKernelArg(ckKernel, 6, sizeof(float4), &advspacing);
	
	ciErrNum |= clSetKernelArg(ckKernel, 7, sizeof(float), &flow_stamps[start_step]);
	ciErrNum |= clSetKernelArg(ckKernel, 10, sizeof(float), &flow_stamps[start_step + 1]);
	
	ciErrNum |= clSetKernelArg(ckKernel, 13, sizeof(float), &st);
	ciErrNum |= clSetKernelArg(ckKernel, 14, sizeof(float), &et);
	
	uint boundary = 0;
	ciErrNum |= clSetKernelArg(ckKernel, 15, sizeof(uint), &boundary);
	
	ciErrNum |= clSetKernelArg(ckKernel, 8, sizeof(cl_mem), (void *) &d_volumeFlow0);
	ciErrNum |= clSetKernelArg(ckKernel, 11, sizeof(cl_mem), (void *) &d_volumeFlow1);
	ciErrNum |= clSetKernelArg(ckKernel, 9, sizeof(cl_sampler), &volumeSampler0);
	ciErrNum |= clSetKernelArg(ckKernel, 12, sizeof(cl_sampler), &volumeSampler1);
	
	// execute kernel
	for (int p = 0; p < depth; p++)
	{
		ciErrNum |= clSetKernelArg(ckKernel, 16, sizeof(uint), &p);
		
		size_t localSize[] = { width};
		size_t gridSize[] = { width * height };
		ciErrNum |= clEnqueueNDRangeKernel(cqCommandQueue, ckKernel, 1, NULL, gridSize, localSize, 0, 0, 0);
		oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
		clFinish(cqCommandQueue);
	}
	
	
}

void invokeAdvection(float start_time, float duration)
{
	int start_step = -1;
	bool inrange = false;
	for (uint i = 0; i < steps - 1; i++)
	{
		if ((start_time >= flow_stamps[i]) && (start_time <= flow_stamps[i+1]))
		{
			start_step = i;
			inrange = true;
			break;
		}
	}
	
	if (!inrange)
	{
		printf("Error: start time is out of range.\n");
		return;
	}
	
	float from, to;
	if (duration > 0.0)
	{
		from = to = start_time;
		while ((start_step < steps - 1) && (to < start_time + duration))
		{
			loadFlowTextures(flow_files[start_step], flow_files[start_step + 1]);
			
			from = max(start_time, flow_stamps[start_step]);
			to = min(start_time + duration, flow_stamps[start_step + 1]);
			
			printf("From %f to %f\n", from, to);
			advect(from, to, start_step);
			start_step++;
		}
	}
	else
	{
		from = to = start_time;
		while ((start_step >= 0) && (to > start_time + duration))
		{
			loadFlowTextures(flow_files[start_step], flow_files[start_step + 1]);
			
			from = min(start_time, flow_stamps[start_step + 1]);
			to = max(start_time + duration, flow_stamps[start_step]);
			
			printf("From %f to %f\n", from, to);			
			advect(from, to, start_step);
			start_step--;
		}
	}
}

// Compute the FTLE field based on the flow map gradient
//*****************************************************************************
bool eigenValues(double A[3][3], double* evalr)
{
	double data[] = {A[0][0],A[0][1],A[0][2],A[1][0],A[1][1],A[1][2],A[2][0],A[2][1],A[2][2]};
	gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3);
	gsl_vector *eval = gsl_vector_alloc (3);
	gsl_matrix *evec = gsl_matrix_alloc (3, 3);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
	gsl_eigen_symmv (&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

	evalr[0] = gsl_vector_get (eval, 0);
	evalr[1] = gsl_vector_get (eval, 1);
	evalr[2] = gsl_vector_get (eval, 2);

	gsl_vector_free (eval);
	gsl_matrix_free (evec);

    return true;
}

void computeFTLE(int argc, const char** argv, double duration)
{
	// allocate memory for the FTLE
    float* ftle = (float*) malloc(width * height * depth * sizeof(float));
		
	// compute ftle
	float4 sp = spacing;
	#pragma omp parallel for schedule(dynamic)
	for (uint x = 0; x < width; x++)
	{
		for (uint y = 0; y < height; y++)
		{
			for (uint z = 0; z < depth; z++)
			{
				double m1[9];
				double ev;
				float4 tp1, tp2;
	
				uint addr = x + width * (y + height * z);
				if ((x == width - 1) || (y == height - 1) || (z == depth - 1))
				{
					ftle[addr] = 0.0f;
					continue;
				}

				float sumis = 0.0;
				float4 here = h_flowMap[addr];
				sumis += h_flowMap[addr].w;
				
				///////////////////////////////////////////////////////////////////////////////
				//////// Find the x finite difference part
				tp1 = h_flowMap[addr + 1];
				sumis += h_flowMap[addr + 1].w;
				tp2 = here;

				m1[0] = (tp1.x - tp2.x) / sp.x;
				m1[3] = (tp1.y - tp2.y) / sp.x;
				m1[6] = (tp1.z - tp2.z) / sp.x;

				/////////////////////////////////////////////////////////////////////////////////
				////////// Find the y finite difference part
				tp1 = h_flowMap[addr + width];
				sumis += h_flowMap[addr + width].w;
				tp2 = here;

				m1[1] = (tp1.x - tp2.x) / sp.y;
				m1[4] = (tp1.y - tp2.y) / sp.y;
				m1[7] = (tp1.z - tp2.z) / sp.y;

				/////////////////////////////////////////////////////////////////////////////////
				////////// Find the z finite difference part
				tp1 = h_flowMap[addr + width * height];
				sumis += h_flowMap[addr + width * height].w;
				tp2 = here;

				m1[2] = (tp1.x - tp2.x) / sp.z;
				m1[5] = (tp1.y - tp2.y) / sp.z;
				m1[8] = (tp1.z - tp2.z) / sp.z;

				/////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////
				////////// Initialize S and D
				double wr[3];
				double A[3][3];
				A[0][0] = m1[0] * m1[0] + m1[3] * m1[3] + m1[6] * m1[6];
				A[0][1] = m1[0] * m1[1] + m1[3] * m1[4] + m1[6] * m1[7];
				A[0][2] = m1[0] * m1[2] + m1[3] * m1[5] + m1[6] * m1[8];
				A[1][1] = m1[1] * m1[1] + m1[4] * m1[4] + m1[7] * m1[7];
				A[1][2] = m1[1] * m1[2] + m1[4] * m1[5] + m1[7] * m1[8];
				A[2][2] = m1[2] * m1[2] + m1[5] * m1[5] + m1[8] * m1[8];
				A[1][0] = A[0][1];
				A[2][0] = A[0][2];
				A[2][1] = A[1][2];

				eigenValues(A, wr);
				ev = wr[2];
				/////////////////////////////////////////////////////////////////////////////////

				// set the value in the texture
				ev = 0.1 * log(sqrt(ev)) / duration;

				//if (ev < 0.0)
				//	ev = 0.0;

				//if (sumis > 0.0)
				//	ev = 0.0;

				ftle[addr] = ev;
			}
		}
		//printf("Slice %d complete!\n", x);
	}

	char *filename;
    if (shrGetCmdLineArgumentstr(argc, (const char**)argv, "oftle", &filename)) {
        printf("FTLE file: %s\n", filename);
		writeRawFile(ftle, filename, width, height, depth);
    }

	free(ftle);
	
	printf("Done!");

}

// Main program
//*****************************************************************************
int main(int argc, char** argv) 
{
	cout << "Example: oclFTLE -device=2 -file=/scratch/sbarakat/Datasets/jet4/samer_jet4.timesteps -start=0.110 -length=0.02 -oftle=ftle.nrrd\n";
	
	pArgc = &argc;
	pArgv = argv;

	shrQAStart(argc, argv);

    // start logs
	cExecutableName = argv[0];
    shrSetLogFileName ("log.txt");
    shrLog("%s Starting...\n\n", argv[0]); 

    // Create OpenCL context, get device info, select device, select options for image/texture and CL-GL interop
    createCLContext(argc, (const char**)argv);

    // Print device info
    clGetDeviceInfo(cdDevices[uiDeviceUsed], CL_DEVICE_IMAGE_SUPPORT, sizeof(g_bImageSupport), &g_bImageSupport, NULL);
    shrLog("%s...\n\n", g_bImageSupport ? "Using Image (Texture)" : "No Image (Texture) Support");      
    shrLog("Detailed Device info:\n\n");
    oclPrintDevInfo(LOGBOTH, cdDevices[uiDeviceUsed]);

    // create a command-queue
    cqCommandQueue = clCreateCommandQueue(cxGPUContext, cdDevices[uiDeviceUsed], 0, &ciErrNum);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);

    // Program Setup
    size_t program_length;
    cPathAndName = shrFindFilePath("ftleComputation.cl", argv[0]);
    oclCheckErrorEX(cPathAndName != NULL, shrTRUE, pCleanup);
    cSourceCL = oclLoadProgSource(cPathAndName, "", &program_length);
    oclCheckErrorEX(cSourceCL != NULL, shrTRUE, pCleanup);

    // create the program
    cpProgram = clCreateProgramWithSource(cxGPUContext, 1,
					  (const char **)&cSourceCL, &program_length, &ciErrNum);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);
    
    // build the program
    std::string buildOpts = "-cl-fast-relaxed-math";
    buildOpts += g_bImageSupport ? " -DIMAGE_SUPPORT" : "";
    ciErrNum = clBuildProgram(cpProgram, 0, NULL, buildOpts.c_str(), NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        // write out standard error, Build Log and PTX, then cleanup and return error
        shrLogEx(LOGBOTH | ERRORMSG, ciErrNum, STDERROR);
        oclLogBuildInfo(cpProgram, oclGetFirstDev(cxGPUContext));
        oclLogPtx(cpProgram, oclGetFirstDev(cxGPUContext), "oclVolumeRender.ptx");
        Cleanup(EXIT_FAILURE); 
    }

    // create the kernel
    ckKernel = clCreateKernel(cpProgram, "d_advectParticles", &ciErrNum);
    oclCheckErrorEX(ciErrNum, CL_SUCCESS, pCleanup);

    // parse arguments
    char *filename;
    if (shrGetCmdLineArgumentstr(argc, (const char**)argv, "file", &filename)) {
        printf("Data file: %s\n", filename);
    }
    float start_time;
    if (shrGetCmdLineArgumentf(argc, (const char**)argv, "start", &start_time)) {
        printf("Start time: %f\n", start_time);
    }
	float length_time;
    if (shrGetCmdLineArgumentf(argc, (const char**)argv, "length", &length_time)) {
        printf("Time length: %f\n", length_time);
    }

    // init timer 1 for fps measurement 
    shrDeltaT(1);  
    
	// load data file 
	loadStepsFile(filename);

	// create textures and flow map space
	initCLGPUMemory();
	
    // do advection
	invokeAdvection(start_time, length_time);

	// compute FTLE
	copyGPUFlowMap();
	computeFTLE(argc, (const char**)argv, length_time);

	// Normally unused return path
    Cleanup(EXIT_SUCCESS);
}

// Function to clean up and exit
//*****************************************************************************
void Cleanup(int iExitCode)
{
	// free memory
	free(flow_stamps);
	free(flow_files);
	
    // cleanup allocated objects
    shrLog("\nStarting Cleanup...\n\n");
    if(cPathAndName)free(cPathAndName);
    if(cSourceCL)free(cSourceCL);
	if(ckKernel)clReleaseKernel(ckKernel);  
    if(cpProgram)clReleaseProgram(cpProgram);
    if(volumeSampler0)clReleaseSampler(volumeSampler0);
	if(volumeSampler1)clReleaseSampler(volumeSampler1);
	if(d_flowMap)clReleaseMemObject(d_flowMap);   
	if(d_volumeFlow0)clReleaseMemObject(d_volumeFlow0);   
	if(d_volumeFlow1)clReleaseMemObject(d_volumeFlow1);   
    if(cqCommandQueue)clReleaseCommandQueue(cqCommandQueue);
    if(cxGPUContext)clReleaseContext(cxGPUContext);
	
    // finalize logs and leave
	shrLogEx(LOGBOTH | CLOSELOG, 0, "%s Exiting...\nPress <Enter> to Quit\n", cExecutableName);
	#ifdef WIN32
		getchar();
	#endif

    exit (iExitCode);
}


