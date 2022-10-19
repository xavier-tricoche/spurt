#include <iostream>
#include <sstream>

#include <teem/hest_helper.hpp>
#include <teem/unrrdu.h>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
#include "image/nrrd_field.hpp"
#include "ftle.hpp"
#include "tensor/double_point.hpp"
#include "tensor/eigenvector_field.hpp"

#include <image/nrrd_wrapper.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
double length, sampling, dx, lmin, lmax, lstep;
int dir, scaleLen;
size_t nsamples[3];

#define __EXPORT_FTLE__
#define __EXPORT_ERROR__
#define __EXPORT_ENDPOINT__
#define __EXPORT_LENGTH__

// #define __SLICE__

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
    hestOptAdd(&hopt, "lmin",   "min length",       airTypeDouble,  1,  1,  &lmin,          NULL,   "lower bound on integration length for flow map computation");
    hestOptAdd(&hopt, "lmax",   "max length",       airTypeDouble,  1,  1,  &lmax,          NULL,   "upper bound on integration length for flow map computation");
    hestOptAdd(&hopt, "lstep",  "length step",      airTypeDouble,  1,  1,  &lstep,         NULL,   "step size integration length for flow map computation");
    hestOptAdd(&hopt, "h",      "step size",        airTypeDouble,  1,  1,  &dx,            NULL,   "integration step size");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis",  &scaleLen);
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Compute FTLE in major eigenvector field of 3D symmetric second-order tensor field",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct field_wrapper {
    field_wrapper(const DoublePointLoad& field)
        : dpl(new EigenvectorField<DoublePointLoad>(field)),
          nrrd(0), procedural(true), bbox(field.bounds()) {}
          
    field_wrapper(const nrrd_field<7>& field)
        : dpl(0), nrrd(new EigenvectorField<nrrd_field<7> >(field)),
          procedural(false), bbox(field.bounds()) {}
          
    spurt::vec3 interpolate(const spurt::vec3& x) const {
        if (procedural) {
            return dpl->major(x);
        } else {
            return nrrd->major(x);
        }
    }
    
    const spurt::bbox3& bounds() const {
        return bbox;
    }
    
    const EigenvectorField<DoublePointLoad>*     dpl;
    const EigenvectorField<nrrd_field<7> >*      nrrd;
    const bool                                  procedural;
    spurt::bbox3                                 bbox;
};

int main(int argc, const char* argv[])
{
    hestOpt* hopt;
    initialize(argc, argv, hopt);
    
    file_name name(name_out);
    
    bool procedural = !strcmp(name_in, "none");
    field_wrapper*   efield         = 0;
    DoublePointLoad* dpl            = 0;
    nrrd_field<7>*   nrrd_tensor    = 0;
    if (procedural) {
        dpl = new DoublePointLoad(50, 50, 20);
        dpl->set_check_inside(false);
        efield = new field_wrapper(*dpl);
        std::cerr << "generating a synthetic double point load tensor field" << std::endl;
    } else {
        nrrd_tensor = new nrrd_field<7>(name_in);
        efield = new field_wrapper(*nrrd_tensor);
        std::cerr << "processing nrrd file: " << name_in << std::endl;
    }
    
    spurt::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    RasterGrid<3> sampling_grid(res, efield->bounds());
    RasterData<spurt::vec3, 3> flowmaps[2]
        = { RasterData<spurt::vec3, 3>(sampling_grid),
            RasterData<spurt::vec3, 3>(sampling_grid)
          };
          
    // initialize flow map to seed locations
    for (int i = 0 ; i < sampling_grid.size() ; ++i) {
        spurt::ivec3 idx = sampling_grid.coord(i);
        spurt::vec3 seed = sampling_grid(idx);
        flowmaps[0](idx) = flowmaps[1](idx) = seed;
    }
    
    spurt::timer timer;
    int npoints = sampling_grid.size();
    
    const spurt::vec3 zero(0.);
    
    std::cout << "nb points = " << npoints << '\n';
    
    double lastl = 0.;
    int nlengths = 0;
    std::vector<bool> blocked_fwd(sampling_grid.size(), false);
    std::vector<bool> blocked_bwd(sampling_grid.size(), false);
    std::vector<double> last_step_fwd(sampling_grid.size(), dx);
    std::vector<double> last_step_bwd(sampling_grid.size(), -dx);
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    int nblocked = 0;
    
    for (length = lmin ; length <= lmax ; length += lstep, ++nlengths) {
    
        std::cerr << nblocked << " blocked trajectories so far" << std::endl;
        
#ifdef __EXPORT_FTLE__
        float* ftle = (float*)calloc(npoints, sizeof(float));
        std::cerr << "FTLE values will be exported" << std::endl;
#endif
        
#ifdef __EXPORT_ENDPOINT__
        float* endpt_f = (float*)calloc(3 * npoints, sizeof(float));
        float* endpt_b = (float*)calloc(3 * npoints, sizeof(float));
        spurt::lexicographical_order lexorder;
        std::cerr << "top reached positions will be exported" << std::endl;
#endif
        
#ifdef __EXPORT_ERROR__
        int* err_f = (int*)calloc(3 * npoints, sizeof(int));
        int* err_b = (int*)calloc(3 * npoints, sizeof(int));
        std::cerr << "error code will be exported" << std::endl;
#endif
        
#ifdef __EXPORT_LENGTH__
        float* len_f = (float*)calloc(npoints, sizeof(int));
        float* len_b = (float*)calloc(npoints, sizeof(int));
        std::cerr << "integration length will be exported" << std::endl;
#endif
        
#ifdef __SLICE__
        int nb_layers = sampling_grid.resolution()[2];
        int nb_pts_per_layer = sampling_grid.resolution()[0] * sampling_grid.resolution()[1];
#endif
        
        double requested_length = length - lastl;
        
        std::cout << "advancing flow map from length " << lastl << " to length " << length << std::endl;
        
        int nprocessed = 0;
        int lastpct = -1;
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
                    
#if _OPENMP
                    const int thread = omp_get_thread_num();
                    if (!thread && pct > lastpct) {
#else
                    if (pct > lastpct) {
#endif
                        lastpct = pct;
                        std::cerr << '\r' << pct << "% completed in "
                                  << timer.elapsed() << "s.                         " << std::flush;
                    }
                    spurt::ivec3 c = sampling_grid.coord(n);
                    
                    for (int dir = 0 ; dir < 2 ; ++dir) {
                        if ((!dir && blocked_fwd[n]) || (dir && blocked_bwd[n])) {
                            continue;
                        }
                        
                        const spurt::vec3& seed = flowmaps[dir](c);
                        int error = -3;
                        std::vector<double>& last_step = (dir ? last_step_bwd : last_step_fwd);
                        std::vector<bool>& blocked = (dir ? blocked_bwd : blocked_fwd);
                        RasterData<spurt::vec3, 3>& fmap = flowmaps[dir];
                        try {
                            double length_io = requested_length;
                            spurt::vec3 z = ftle::eigen_flow_map(*efield, seed, last_step[n], length_io, error);
                            bool valok = true;
                            if (std::isinf(spurt::norm(z)) || std::isnan(spurt::norm(z))) {
                                blocked[n] = true;
                                valok = false;
                                ++nblocked;
                            } else {
                                fmap(c) = z;
                            }
#ifdef __EXPORT_LENGTH__
                            float* len = (dir ? len_b : len_f);
                            len[n] = length_io;
#endif
                        } catch (...) {
                            blocked[n] = true;
                            ++nblocked;
                        }
                        
                        if (error != 0) {
                            blocked[n] = true;
                            ++nblocked;
                        }
                        
#ifdef __EXPORT_ERROR__
                        int* err = (dir ? err_b : err_f);
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
#endif
                    }
#ifdef __EXPORT_ENDPOINT__
                    spurt::vec3 end_f, end_b;
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
        std::cout << "\n\ntotal computation time for eigen flow map was " << timer.elapsed() << '\n';
        
        std::vector<size_t> size_4d(4), size_3d(3);
        std::vector<double> step_4d(4), step_3d(3);
        const spurt::vec3& s = sampling_grid.step();
        step_4d[0] = airNaN();
        for (int i = 0 ; i < 3 ; ++i) {
            size_3d[i] = size_4d[i+1] = res[i];
            step_3d[i] = step_4d[i+1] = s[i];
        }
        
        
#ifdef __EXPORT_ENDPOINT__
        size_4d[0] = 3;
        spurt::writeNrrd(endpt_f, name.make("endpoints-fwd", nlengths, length), nrrdTypeFloat, size_4d, step_4d);
        delete[] endpt_f;
        spurt::writeNrrd(endpt_b, name.make("endpoints-bwd", nlengths, length), nrrdTypeFloat, size_4d, step_4d);
        delete[] endpt_b;
#endif
        
#ifdef __EXPORT_ERROR__
        size_4d[0] = 3;
        spurt::writeNrrd(err_f, name.make("error-fwd", nlengths, length), nrrdTypeInt, size_4d, step_4d);
        delete[] err_f;
        spurt::writeNrrd(err_b, name.make("error-bwd", nlengths, length), nrrdTypeInt, size_4d, step_4d);
        delete[] err_b;
#endif
        
#ifdef __EXPORT_LENGTH__
        spurt::writeNrrd(len_f, name.make("length-fwd", nlengths, length), nrrdTypeFloat, size_3d, step_3d);
        delete[] len_f;
        spurt::writeNrrd(len_b, name.make("length-bwd", nlengths, length), nrrdTypeFloat, size_3d, step_3d);
        delete[] len_b;
#endif
        
#ifdef __EXPORT_FTLE__
        timer.restart();
#pragma openmp parallel for
        for (int n = 0 ; n < npoints ; ++n) {
            spurt::vec2 ans;
            try {
                ans = ftle::eigenftle(n, *efield, flowmaps, length);
            } catch (...) {
                blocked_fwd[n] = blocked_bwd[n] = true;
                continue;
            }
            if (std::isinf(ans[0]) || std::isinf(ans[1]) ||
                    std::isnan(ans[0]) || std::isnan(ans[1])) {
                blocked_fwd[n] = blocked_bwd[n] = true;
                continue;
            }
            ftle[n] = std::max(ans[0], ans[1]);
        }
        std::cout << "total computation time for ftle was " << timer.elapsed() << '\n';
        spurt::writeNrrd(ftle, name.make("ftle", nlengths, length), nrrdTypeFloat, size_3d, step_3d);
        delete[] ftle;
#endif
        
        lastl = length;
    }
    
    if (procedural) {
        delete dpl;
    } else {
        delete nrrd_tensor;
    }
    
    return 0;
}
