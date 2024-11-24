#include <iostream>

#include <teem/hest.h>
#include <teem/unrrdu.h>

#include <math/types.hpp>
#include <misc/progress.hpp>
#include <image/nrrd_field.hpp>
#include "ftle.hpp"
#include <tensor/double_point.hpp>
#include <tensor/eigenvector_field.hpp>

#include <image/nrrd_wrapper.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
double length, sampling, dx, rel, min_length, max_length, d_length;
int dir, scaleLen, eigen;
bool eval_weight;
size_t nsamples[3];
char* tl_out;
int nb_seeds;

typedef Eigen::Matrix<double, 7, 1> tensor_type;
typedef spurt::nrrd_field< tensor_type, 3, double > tensor_field_type;

#define __EXPORT_FTLE__
// #define __EXPORT_SEED__
#define __EXPORT_EVEC__
// #define __EXPORT_ERROR__
// #define __EXPORT_ENDPOINT__
#define __EXPORT_LENGTH__

#define __SLICE__

using namespace spurt;

void initialize(int argc, const char* argv[], hestOpt* hopt)
{
    hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  0,  1,  &name_in,       "none", "input file name, if none is provided, a synthetic double point load dataset will be procedurally generated");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    // hestOptAdd(&hopt, "min", "min length",       airTypeDouble,  1,  1,  &min_length,    NULL,   "min integration length for flow map computation");
    // hestOptAdd(&hopt, "max", "max length",       airTypeDouble,  1,  1,  &max_length,    NULL,   "max integration length for flow map computation");
    // hestOptAdd(&hopt, "step",    "length step",      airTypeDouble,  1,  1,  &d_length,      NULL,   "integration length interval for flow map computation");
    hestOptAdd(&hopt, "l",      "length",           airTypeDouble,  1,  1,  &length,        NULL,   "integration length for flow map computation");
    hestOptAdd(&hopt, "h",      "step size",        airTypeDouble,  1,  1,  &dx,            NULL,   "integration step size");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis",  &scaleLen);
    hestOptAdd(&hopt, "r",      "dpl rel distance", airTypeDouble,  0,  1,  &rel,           "0.5",  "relative distance between single point loads in procedural double point load model. Ignored if an input file is selected");
    hestOptAdd(&hopt, "e",      "eigenvector",      airTypeInt,     0,  1,  &eigen,         "0",    "eigenvector field along which integration takes place");
    hestOptAdd(&hopt, "w",      "eval weight",      airTypeBool,    0,  1,  &eval_weight,   "0",    "use associated eigenvalue as norm for the integration");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE in eigenvector field of 3D symmetric second-order tensor field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct field_wrapper {
    field_wrapper(const DoublePointLoad& field)
        : dpl(new EigenvectorField<DoublePointLoad>(field)),
          nrrd(0), procedural(true), bbox(field.bounds()) {}
          
    field_wrapper(const tensor_field_type& field)
        : dpl(0), nrrd(new EigenvectorField<tensor_field_type>(field)),
          procedural(false), bbox(field.bounds()) {}
          
    vec3 interpolate(const vec3& x) const {
        vec4 ev;
        if (procedural) {
            ev = dpl->eigen(x, eigen);
        } else {
            ev = nrrd->eigen(x, eigen);
        }
        vec3 e(ev[0], ev[1], ev[2]);
        if (eval_weight) {
            e *= ev[3];
        }
        return e;
    }
    
    const bbox3& bounds() const {
        return bbox;
    }
    
    const EigenvectorField<DoublePointLoad>*    dpl;
    const EigenvectorField<tensor_field_type>*  nrrd;
    const bool                                  procedural;
    bbox3                                       bbox;
};

int main(int argc, const char* argv[])
{
    hestOpt* hopt;
    initialize(argc, argv, hopt);
    
    if (eigen < 0 || eigen > 2) {
        std::cerr << "ERROR: invalid eigenvector field selected" << std::endl;
        hestUsage(stderr, hopt, argv[0], 0);
        return -1;
    }
    
    bool procedural = !strcmp(name_in, "none");
    field_wrapper*   efield            = 0;
    DoublePointLoad* dpl               = 0;
    tensor_field_type*  nrrd_tensor    = 0;
    if (procedural) {
        dpl = new DoublePointLoad(50, 50, 20, rel);
        dpl->set_check_inside(false);
        efield = new field_wrapper(*dpl);
        std::cerr << "generating a synthetic double point load tensor field" << std::endl;
    } else {
        nrrd_tensor = new tensor_field_type(name_in);
        efield = new field_wrapper(*nrrd_tensor);
        std::cerr << "processing nrrd file: " << name_in << std::endl;
    }
    
    ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    rgrid3 sampling_grid(res, efield->bounds());
    spurt::raster3d<vec3> flowmaps[2]
        = { spurt::raster3d<vec3>(sampling_grid),
            spurt::raster3d<vec3>(sampling_grid)
          };
          
    spurt::timer timer;
    int npoints = sampling_grid.size();
    
    const vec3 zero = 0;
    for (int i = 0 ; i < npoints ; ++i) {
        ivec3 c = flowmaps[0].grid()[i];
        flowmaps[0](c) = flowmaps[1](c) = zero;
    }
    
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
#ifdef __EXPORT_EVEC__
    float* evec = (float*)calloc(3 * npoints, sizeof(float));
    std::cerr << "eigenvectors will be exported" << std::endl;
#endif
    
#ifdef __EXPORT_SEED__
    float* start = (float*)calloc(3 * npoints, sizeof(float));
    std::cerr << "seed point coordinates will be exported" << std::endl;
#endif
    
#ifdef __EXPORT_FTLE__
    float* ftle = (float*)calloc(npoints, sizeof(float));
    std::cerr << "FTLE values will be exported" << std::endl;
#endif
    
#ifdef __EXPORT_ENDPOINT__
    float* endpt_f = (float*)calloc(3 * npoints, sizeof(float));
    float* endpt_b = (float*)calloc(3 * npoints, sizeof(float));
    lexicographical_order lexorder;
    std::cerr << "top reached positions will be exported" << std::endl;
#endif
    
#ifdef __EXPORT_ERROR__
    int* err = (int*)calloc(3 * npoints, sizeof(int));
    std::cerr << "error code will be exported" << std::endl;
#endif
    
#ifdef __EXPORT_LENGTH__
    int* len = (int*)calloc(npoints, sizeof(int));
    std::cerr << "integration length will be exported" << std::endl;
#endif
    
#ifdef __SLICE__
    int nb_layers = sampling_grid.resolution()[2];
    int nb_pts_per_layer = sampling_grid.resolution()[0] * sampling_grid.resolution()[1];
#endif
    
    int nprocessed = 0;
    int counter = 0;
    
    int nsteps_necessary = (int)(length/dx);
    int nmax = 2*nsteps_necessary;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        
#ifdef __SLICE__
        for (int l = 0 ; l < nb_layers ; ++l) {
            for (int m = 0 ; m < nb_pts_per_layer ; ++m) {
                // progress tracking
                ++nprocessed;
                
                int n = m + l * nb_pts_per_layer;
                int pct = 100 * nprocessed / npoints;
#else
        for (int n = 0 ; n < npoints ; ++n) {  {
                int pct = 100 * n / npoints;
#endif
                
                int thread = 0;
                int __c = n;
#if _OPENMP
                thread = omp_get_thread_num();
                #pragma omp atomic
                ++counter;
                
                __c = counter;
#endif
                
                ivec3 c = sampling_grid[n];
                vec3 seed = sampling_grid(c);
                
                if (!thread) {
                    std::cerr << '\r' << 100*__c / npoints << "% completed in "
                              << timer.elapsed()
                              << "s. "
                              << "integrator at " << seed << "                 \r"
                              << std::flush;
                }
                
                
#ifdef __EXPORT_SEED__
                for (int k = 0 ; k < 3 ; ++k) {
                    start[3*n+k] = seed[k];
                }
#endif
                
#ifdef __EXPORT_EVEC__
                {
                    vec3 ev;
                    try {
                        ev = efield->interpolate(seed);
                    } catch (...) {
                        ev = vec3(0, 0, 0);
                    }
                    for (int k = 0 ; k < 3 ; ++k) {
                        evec[3*n+k] = fabs(ev[k]);
                    }
                }
#endif
                for (int dir = 0 ; dir < 2 ; ++dir) {
                    int error = -3;
                    double h = (dir ? dx : -dx);
                    try {
                        double length_io = length;
                        flowmaps[dir](c) = ftle::eigen_flow_map(*efield, seed, h, length_io, error, nmax);
                        const vec3& z = flowmaps[dir](c);
                        if (any(isinvalid(z)) {
                            flowmaps[dir](c) = zero;
                        }
#ifdef __EXPORT_LENGTH__
                        else if (!dir) {
                            len[n] = length_io;
                        }
#endif
                    } catch (...) {
                        flowmaps[dir](c) = zero;
                    }
                    
#ifdef __EXPORT_ERROR__
                    if (!dir) {
                        switch (error) {
                            case 0:
                                err[3*n+1] = 255;
                                break;
                            case 1:
                                err[3*n+2] = 255;
                                break;
                            case 100:
                                err[3*n] = 255;
                            default:
                                break;
                        }
                    }
#endif
                }
#ifdef __EXPORT_ENDPOINT__
                vec3 end_f, end_b;
                if (lexorder(flowmaps[0](c), flowmaps[1](c))) {
                    end_f = flowmaps[1](c);
                    end_b = flowmaps[0](c);
                } else {
                    end_f = flowmaps[0](c);
                    end_b = flowmaps[1](c);
                }
                for (int k = 0 ; k < 3 ; ++k) {
                    endpt_f[3*n+k] = end_f[k];
                    endpt_b[3*n+k] = end_b[k];
                }
#endif
            }
        }
    }
    
    std::cout << "\ntotal computation time for eigen flow map was " << timer.elapsed() << '\n';
    
    spurt::nrrd_utils::nrrd_params<int, 3>    p3i;
    spurt::nrrd_utils::nrrd_params<int, 4>    p4i;
    spurt::nrrd_utils::nrrd_params<float, 3>  p3f;
    spurt::nrrd_utils::nrrd_params<float, 4>  p4f;
    const vec3& s = sampling_grid.spacing();
    p4i.sizes()[0] = p4f.sizes()[0] = airNaN();
    for (int i = 0 ; i < 3 ; ++i) {
        p3i.sizes()[i] =
        p3f.sizes()[i] = 
        p4i.sizes()[i+1] = 
        p4f.sizes()[i+1] = res[i];
        p3i.spacings()[i] =
        p3f.spacings()[i] = 
        p4i.spacings()[i+1] = 
        p4f.spacings()[i+1] = s[i];
        p3i.mins()[i] =
        p3f.mins()[i] = 
        p4i.mins()[i+1] = 
        p4f.mins()[i+1] = efield->bounds().min()[i];
    }
    
#ifdef __EXPORT_ENDPOINT__
    spurt::nrrd_utils::writeNrrd(endpt_f, "endpoints_f.nrrd", p4f);
    spurt::nrrd_utils::writeNrrd(endpt_b, "endpoints_b.nrrd", p4f);
#endif
    
#ifdef __EXPORT_EVEC__
    spurt::nrrd_utils::writeNrrdFromParams(evec, "evec0.nrrd", p4f);
#endif
    
#ifdef __EXPORT_SEED__
    spurt::nrrd_utils::writeNrrdFromParams(start, "seeds.nrrd", p4f);
#endif
    
#ifdef __EXPORT_ERROR__
    spurt::nrrd_utils::writeNrrdFromParams(err, "error.nrrd", p4i);
#endif
    
#ifdef __EXPORT_LENGTH__
    spurt::nrrd_utils::writeNrrdFromParams(len, "length.nrrd", p3i);
#endif
    
#ifdef __EXPORT_FTLE__
    timer.restart();
    counter = 0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n = 0 ; n < npoints ; ++n) {
            int thread = 0;
            int __c = n;
#if _OPENMP
            thread = omp_get_thread_num();
            #pragma omp atomic
            ++counter;
            
            __c = counter;
#endif
            if (!thread) {
                int pct = 100 * __c / npoints;
                std::cerr << '\r' << pct << "% completed in "
                          << timer.elapsed()
                          << "s.                 \r"
                          << std::flush;
            }
            
            vec2 ans;
            try {
                ans = ftle::eigenftle(n, *efield, flowmaps, length);
            } catch (...) {
                continue;
            }
            if (std::isinf(ans[0]) || std::isinf(ans[1]) ||
                    std::isnan(ans[0]) || std::isnan(ans[1])) {
                continue;
            }
            ftle[n] = std::max(ans[0], ans[1]);
        }
    }
    std::cout << "\ntotal computation time for ftle was " << timer.elapsed() << '\n';
    
    spurt::nrrd_utils::writeNrrdFromParams(ftle, name_out, p3f);
#endif
    
    if (procedural) {
        delete dpl;
    } else {
        delete nrrd_tensor;
    }
    
    return 0;
}






















































































































































