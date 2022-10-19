#include <iostream>
#include <image/nrrd_wrapper.hpp>
#include <teem/hest.h>
#include <math/fixed_vector.hpp>
#include <exception>

char* name_in, *name_out;
double length, eps;
size_t res[2];
double spc[2];

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1,  1,  &name_in,       NULL,   "input file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    hestOptAdd(&hopt, "l",      "length",           airTypeDouble,  1,  1,  &length,        NULL,   "integration length");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  1,  1,  &eps,           NULL,   "step size");
    hestOptAdd(&hopt, "s",      "size",             airTypeSize_t,  2,  2,  res,            NULL,   "image resolution");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "LIC computation over 2D NRRD vector field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T>
struct texture {
    texture(int width, int height) : __w(width), __h(height) {
        __tex = (T*)calloc(__w*__h*sizeof(T));
    }
    
    const T& operator()(int i, int j) const {
        if (i>=0 & i<__w && j>=0 && j<__h) {
            return __tex[i+j*__w];
        }
        throw std::runtime_error("out of bounds");
    }
    
    T& operator()(int i, int j) {
        if (i>=0 & i<__w && j>=0 && j<__h) {
            return __tex[i+j*__w];
        }
        throw std::runtime_error("out of bounds");
    }
    
    ~texture() {
        if (__tex) {
            delete[] __tex;
        }
    }
    
    int __w, __h;
    T* __tex;
};

// sample noise texture along prescribed segment using Bresenham algorithm
float bresenham(const fvec2& x, const fvec2& y,
                __read_only image2d_t noise, sampler_t noise_s,
                const int2 dim, const float4 bounds, const float2 step)
{
    int2 start = world_to_index(x, dim, bounds, step);
    if (!valid_pos(start)) {
        return 0;
    }
    int2 end = world_to_index(y, dim, bounds, step);
    if (!valid_pos(end)) {
        return 0;
    }
    
    int dx = abs(end.x - start.x);
    int dy = abs(end.y - start.y);
    if (dx == 0 && dy == 0) {
        return 0;
    }
    int sx = 1;
    int sy = 1;
    if (end.x < start.x) {
        sx = -1;
    }
    if (end.y < start.y) {
        sy = -1;
    }
    int err = dx-dy;
    
    float sum = 0;
    
    for (int i=0 ; true ; ++i) {
        float4 n = read_imagef(noise, noise_s, float2(x0+0.5f, y0+0.5f));
        if (i) {
            sum += n.x;    // skip noise value at segment's start to avoid repeats
        }
        if (x0 == x1 && y0 == y1) {
            break;
        }
        int e2 = 2*err;
        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y0 += sy;
        }
    }
}

void integrate(float2 start, __global float2* lic_image
               __read_only image2d_t flow, sampler_t flow_spl,
               __read_only image2d_t noise, sampler_t noise_spl
               const float h, const int nsteps,
               const int2 dim, const float4 bounds, const float2 step)
{

    float2 k1, k2, x=start;
    float sum=0;
    bool ok = true;
    for (int i=0 ; i<nsteps && ok ; ++i) {
        // RK2
        k1 = h*rhs(x, flow, flow_spl);
        k2 = h*rhs(x+0.5f*k1, flow, flow_spl);
        if (!k1.z || !k2.z) {
            ok = false;
        } else {
            sum += bresenham(x, x+k2, noise, noise_spl);
            x += k2;
        }
    }
    *lic_image += sum;
}

__kernel void
filter(__global float* lic_image, uint width, uint height,
       const int2 dim, const float4 bounds, const float2 step, float h, int n,
       __read_only image2d_t rhs, sampler_t rhs_s,
       __read_only image2d_t noise, sampler_t noise_s)
{

    uint i = get_global_id(0);
    float2 seed = int_to_world(i, dim, bounds, step);
    int2 idx = int_to_index(i, dim);
    float4 n = read_imagef(noise, noise_s, float2(idx.x+0.5f, idx.y+0.5f));
    *lic_image = n.x;
    integrate(seed, lic_image, rhs, rhs_s, noise, noise_s, h, nsteps, dim, bounds, step);
    integrate(seed, lic_image, rhs, rhs_s, noise, noise_s, -h, nsteps, dim, bounds, step);
}
