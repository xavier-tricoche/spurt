#include <iostream>

#include <teem/hest_helper.hpp>
#include <teem/unrrdu.h>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
#include "nrrd_field.hpp"
#include "ftle.hpp"
#include "double_point.hpp"
#include "eigenvector_field.hpp"

#include <image/nrrd_wrapper.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
double length, sampling, dx, rel, min_length, max_length, d_length;
int dir, scaleLen, eigen;
size_t nsamples[3];

#define __EXPORT_FTLE__
// #define __EXPORT_ERROR__
// #define __EXPORT_ENDPOINT__
// #define __EXPORT_LENGTH__

#define __SLICE__

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
    hestOptAdd(&hopt, "min",    "min length",       airTypeDouble,  1,  1,  &min_length,    NULL,   "min integration length for flow map computation");
    hestOptAdd(&hopt, "max",    "max length",       airTypeDouble,  1,  1,  &max_length,    NULL,   "max integration length for flow map computation");
    hestOptAdd(&hopt, "step",   "length step",      airTypeDouble,  1,  1,  &d_length,      NULL,   "integration length interval for flow map computation");
    // hestOptAdd(&hopt, "l",       "length",           airTypeDouble,  1,  1,  &length,        NULL,   "integration length for flow map computation");
    hestOptAdd(&hopt, "h",      "step size",        airTypeDouble,  1,  1,  &dx,            NULL,   "integration step size");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis",  &scaleLen);
    hestOptAdd(&hopt, "r",      "dpl rel distance", airTypeDouble,  0,  1,  &rel,           "0.5",  "relative distance between single point loads in procedural double point load model. Ignored if an input file is selected");
    hestOptAdd(&hopt, "e",      "eigenvector",      airTypeInt,     0,  1,  &eigen,         "0",    "eigenvector field along which integration takes place");
    
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
          
    nvis::vec3 interpolate(const nvis::vec3& x) const {
        if (procedural)
            switch (eigen) {
                case 0:
                    return dpl->major(x);
                case 1:
                    return dpl->medium(x);
                case 2:
                    return dpl->minor(x);
                default:
                    assert(false);
            }
        else
            switch (eigen) {
                case 0:
                    return nrrd->major(x);
                case 1:
                    return nrrd->medium(x);
                case 2:
                    return nrrd->minor(x);
                default:
                    assert(false);
            }
    }
    
    const nvis::bbox3& bounds() const {
        return bbox;
    }
    
    const EigenvectorField<DoublePointLoad>*     dpl;
    const EigenvectorField<nrrd_field<7> >*      nrrd;
    const bool                                  procedural;
    nvis::bbox3                                 bbox;
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
    field_wrapper*   efield         = 0;
    DoublePointLoad* dpl            = 0;
    nrrd_field<7>*   nrrd_tensor    = 0;
    if (procedural) {
        dpl = new DoublePointLoad(50, 50, 20, rel);
        dpl->set_check_inside(false);
        efield = new field_wrapper(*dpl);
        std::cerr << "generating a synthetic double point load tensor field" << std::endl;
    } else {
        nrrd_tensor = new nrrd_field<7>(name_in);
        efield = new field_wrapper(*nrrd_tensor);
        std::cerr << "processing nrrd file: " << name_in << std::endl;
    }
    
    nvis::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    RasterGrid<3> sampling_grid(res, efield->bounds());
    RasterData<nvis::vec3, 3> flowmaps[2]
        = { RasterData<nvis::vec3, 3>(sampling_grid),
            RasterData<nvis::vec3, 3>(sampling_grid)
          };
          
    nvis::timer timer;
    int npoints = sampling_grid.size();
    
    const nvis::vec3 zero(0.);
    for (int i = 0 ; i < npoints ; ++i) {
        nvis::ivec3 c = flowmaps[0].grid().coord(i);
        flowmaps[0](c) = flowmaps[1](c) = zero;
    }
    
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
#ifdef __SLICE__
    int nb_layers = sampling_grid.resolution()[2];
    int nb_pts_per_layer = sampling_grid.resolution()[0] * sampling_grid.resolution()[1];
#endif
    
    int nprocessed = 0;
    
    std::vector<std::vector<bool> > last_dir(2);
    last_dir[0].resize(npoints, true);
    last_dir[1].resize(npoints, false);
    
    std::vector<std::vector<nvis::vec3> > last_pos(2);
    last_pos[0].resize(npoints);
    last_pos[1].resize(npoints);
    for (int i = 0 ; i < npoints ; ++i) {
        nvis::ivec3 c = sampling_grid.coord(i);
        last_pos[0][i] = last_pos[1][i] = sampling_grid(c);
    }
    std::vector<std::vector<bool> > valid(2);
    valid[0].resize(npoints, true);
    valid[1].resize(npoints, true);
    
    file_name name(name_out);
    
    double current_length, next_length;
    current_length = 0;
    next_length = min_length;
    
    for (unsigned int frame_id = 0; current_length <= max_length ; ++frame_id, next_length += d_length) {
    
        timer.restart();
        
#ifdef __EXPORT_FTLE__
        float* ftle = (float*)calloc(npoints, sizeof(float));
        std::cerr << "FTLE values will be exported" << std::endl;
#endif
        
#ifdef __EXPORT_ENDPOINT__
        float* endpt_f = (float*)calloc(3 * npoints, sizeof(float));
        float* endpt_b = (float*)calloc(3 * npoints, sizeof(float));
        nvis::lexicographical_order lexorder;
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
                                  << timer.elapsed()
                                  << "s.                 \r"
                                  << std::flush;
                    }
                    
                    nvis::ivec3 c = sampling_grid.coord(n);
                    
                    for (int dir = 0 ; dir < 2 ; ++dir) {
                        if (!valid[dir][n]) {
                            continue;
                        }
                        int error = -3;
                        double h = (last_dir[dir][n] ? dx : -dx);
                        try {
                            bool fwd;
                            double length_io = next_length - current_length;
                            flowmaps[dir](c) = ftle::eigen_flow_map(*efield, last_pos[dir][n], h,
                                                                    length_io, error, fwd);
                            last_dir[dir][n] = fwd;
                            const nvis::vec3& z = flowmaps[dir](c);
                            if (std::isinf(nvis::norm(z)) || std::isnan(nvis::norm(z))) {
                                flowmaps[dir](c) = zero;
                                valid[dir][n] = false;
                            }
                            
#ifdef __EXPORT_LENGTH__
                            else if (!dir) {
                                len[n] = length_io;
                            }
#endif
                            if (length_io < next_length - current_length) {
                                valid[dir][n] = false;
                            }
                            
                            last_pos[dir][n] = z;
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
                    nvis::vec3 end_f, end_b;
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
        
        std::vector<size_t> size(4), size3(3);
        std::vector<double> step(4), step3(3);
        const nvis::vec3& s = sampling_grid.step();
        step[0] = airNaN();
        for (int i = 0 ; i < 3 ; ++i) {
            size3[i] = size[i+1] = res[i];
            step3[i] = step[i+1] = s[i];
        }
        
#ifdef __EXPORT_ENDPOINT__
        size[0] = 3;
        spurt::nrrd_utils::writeNrrd(endpt_f, name.make("endpoints_f.nrrd", frame_id, next_length), nrrdTypeFloat, size, step);
        spurt::nrrd_utils::writeNrrd(endpt_b, name.make("endpoints_b.nrrd", frame_id, next_length), nrrdTypeFloat, size, step);
#endif
        
#ifdef __EXPORT_ERROR__
        size[0] = 3;
        spurt::nrrd_utils::writeNrrd(err, name.make("error.nrrd", frame_id, next_length), nrrdTypeInt, size, step);
#endif
        
#ifdef __EXPORT_LENGTH__
        {
            std::vector< size_t > _size(3);
            std::vector<double> _spac(3);
            for (int k = 0 ; k < 3 ; ++k) {
                _size[k] = size[k+1];
                _spac[k] = step[k+1];
            }
            spurt::nrrd_utils::writeNrrd(len, name.make("length.nrrd", frame_id, next_length), nrrdTypeInt, _size, _spac);
        }
#endif
        
#ifdef __EXPORT_FTLE__
        timer.restart();
#pragma openmp parallel for
        for (int n = 0 ; n < npoints ; ++n) {
            nvis::vec2 ans;
            try {
                ans = ftle::eigenftle(n, *efield, flowmaps, next_length);
            } catch (...) {
                continue;
            }
            if (std::isinf(ans[0]) || std::isinf(ans[1]) ||
                    std::isnan(ans[0]) || std::isnan(ans[1])) {
                // std::cerr << "skipping invalid value #" << n << std::endl;
                continue;
            }
            // std::cerr << "valid FTLE values at vertex #" << n << std::endl;
            ftle[n] = std::max(ans[0], ans[1]);
        }
        std::cout << "total computation time for ftle was " << timer.elapsed() << '\n';
        
        spurt::nrrd_utils::writeNrrd(ftle, name.make("ftle", frame_id, next_length), nrrdTypeFloat, size3, step3);
#endif
        
        current_length = next_length;
    }
    
    if (procedural) {
        delete dpl;
    } else {
        delete nrrd_tensor;
    }
    
    return 0;
}




























































































































































