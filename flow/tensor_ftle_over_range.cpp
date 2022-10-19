#include <iostream>

#include <teem/hest_helper.hpp>
#include <teem/unrrdu.h>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
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
size_t nsamples[3];

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
        : dpl(0), nrrd_7(new EigenvectorField<nrrd_field<7> >(field, true)), _7d(true),
          procedural(false), bbox(field.bounds()) {}
          
    field_wrapper(const nrrd_field<9>& field)
        : dpl(0), nrrd_9(new EigenvectorField<nrrd_field<9> >(field, false)), _7d(false)
        procedural(false), bbox(field.bounds()) {}
        
    spurt::vec3 interpolate(const spurt::vec3& x) const {
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
        else if (_7d)
            switch (eigen) {
                case 0:
                    return nrrd_7->major(x);
                case 1:
                    return nrrd_7->medium(x);
                case 2:
                    return nrrd_7->minor(x);
                default:
                    assert(false);
            }
        else
            switch (eigen) {
                case 0:
                    return nrrd_9->major(x);
                case 1:
                    return nrrd_9->medium(x);
                case 2:
                    return nrrd_9->minor(x);
                default:
                    assert(false);
            }
    }
    
    const spurt::bbox3& bounds() const {
        return bbox;
    }
    
    const EigenvectorField<DoublePointLoad>*     dpl;
    const EigenvectorField<nrrd_field<7> >*      nrrd_7;
    const EigenvectorField<nrrd_field<9> >*      nrrd_9;
    const bool                                  procedural, _7d;
    spurt::bbox3                                 bbox;
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
    
    file_name name(name_out);
    
    bool procedural = !strcmp(name_in, "none");
    field_wrapper*   efield         = 0;
    DoublePointLoad* dpl            = 0;
    nrrd_field<7>*   nrrd_tensor_7D = 0;
    nrrd_field<9>*   nrrd_tensor_9D = 0;
    if (procedural) {
        dpl = new DoublePointLoad(50, 50, 20, rel);
        dpl->set_check_inside(false);
        efield = new field_wrapper(*dpl);
        std::cerr << "generating a synthetic double point load tensor field" << std::endl;
    } else {
        Nrrd* nin = nrrdNew();
        if (nrrdLoad(nin, name_in, NULL)) {
            std::cerr << biffGetDone(NRRD) << std::endl;
            exit(-1);
        }
        size_t sz = nin->axis[0].size;
        nrrdNuke(nin);
        
        if (sz == 7) {
            nrrd_tensor_7D = new nrrd_field<7>(name_in);
            efield = new field_wrapper(*nrrd_tensor_7D, true);
        } else {
            nrrd_tensor_9D = new nrrd_field<9>(name_in);
            efield = new field_wrapper(*nrrd_tensor_9D, false);
        }
        std::cerr << "processing nrrd file: " << name_in << std::endl;
    }
    
    spurt::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    RasterGrid<3> sampling_grid(res, efield->bounds());
    RasterData<spurt::vec3, 3> flowmaps[2]
        = { RasterData<spurt::vec3, 3>(sampling_grid),
            RasterData<spurt::vec3, 3>(sampling_grid)
          };
          
    spurt::timer timer;
    int npoints = sampling_grid.size();
    
    const spurt::vec3 zero(0.);
    for (int i = 0 ; i < npoints ; ++i) {
        spurt::ivec3 c = flowmaps[0].grid().coord(i);
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
    std::vector<bool> valid[2];
    valid[0].resize(npoints, true);
    valid[1].resize(npoints, true);
    
    // data structure to keep track of the intermediate positions reached
    std::vector<std::pair<spurt::vec3, spurt::vec3> > end_pt[2];
    end_pt[0].resize(npoints);
    end_pt[1].resize(npoints);
    for (int i = 0 ; i < npoints ; ++i) {
        spurt::ivec3 c = sampling_grid.coord(i);
        spurt::vec3 x = sampling_grid(c);
        spurt::vec3 e;
        try {
            e = efield->interpolate(x);
        } catch(...) {
            valid[0][i] = valid[1][i] = false;
        }
        end_pt[0][i] = std::make_pair(x, e);
        end_pt[1][i] = std::make_pair(x, -1.*e);
    }
    
    double cur_length, next_length = min_length;
    unsigned int counter = 0;
    for (cur_length = 0 ; cur_length < max_length ; next_length += d_length, ++counter) {
        double length = std::min(max_length, next_length) - cur_length;
        if (length <= 0) {
            break;
        }
        
        float* ftle = (float*)calloc(npoints, sizeof(float));
        std::cerr << "FTLE values will be exported" << std::endl;
        
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
                    
                    spurt::ivec3 c = sampling_grid.coord(n);
                    
                    for (int dir = 0 ; dir < 2 ; ++dir) {
                        if (!valid[dir][n]) {
                            continue;
                        }
                        int error = -3;
                        try {
                            double length_io = length;
                            flowmaps[dir](c) = ftle::eigen_flow_map(*efield, end_pt[dir][n].first, dx,
                                                                    length_io, error, end_pt[dir][n].second);
                            const spurt::vec3& z = flowmaps[dir](c);
                            if (length_io < length) {
                                valid[dir][n] = false;
                            }
                            if (std::isinf(spurt::norm(z)) || std::isnan(spurt::norm(z))) {
                                flowmaps[dir](c) = zero;
                                valid[dir][n] = false;
                            }
#ifdef __EXPORT_LENGTH__
                            else if (!dir) {
                                len[n] = length_io;
                            }
#endif
                        } catch (...) {
                            flowmaps[dir](c) = zero;
                        }
                        
                        end_pt[dir][n].first = flowmaps[dir](c);
                    }
                }
            }
        }
        
        std::cout << "\ntotal computation time for eigen flow map was " << timer.elapsed() << '\n';
        
        std::vector<size_t> size(4), size3(3);
        std::vector<double> step(4), step3(3);
        const spurt::vec3& s = sampling_grid.step();
        step[0] = airNaN();
        for (int i = 0 ; i < 3 ; ++i) {
            size3[i] = size[i+1] = res[i];
            step3[i] = step[i+1] = s[i];
        }
        
#ifdef __EXPORT_LENGTH__
        {
            std::vector< size_t > _size(3);
            std::vector<double> _spac(3);
            for (int k = 0 ; k < 3 ; ++k) {
                _size[k] = size[k+1];
                _spac[k] = step[k+1];
            }
            spurt::writeNrrd(len, name.make("length.nrrd", counter, cur_length), nrrdTypeInt, _size, _spac);
        }
#endif
        
        timer.restart();
#pragma openmp parallel for
        for (int n = 0 ; n < npoints ; ++n) {
            spurt::vec2 ans;
            try {
                ans = ftle::eigenftle(n, *efield, flowmaps, length);
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
        
        spurt::writeNrrd(ftle, name.make("ftle", counter, cur_length), nrrdTypeFloat, size3, step3);
        
        cur_length = next_length;
    }
    
    if (procedural) {
        delete dpl;
    } else {
        delete nrrd_tensor;
    }
    
    return 0;
}
