#include <iostream>

#include <teem/hest_helper.hpp>
#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/wall_timer.hpp>
#include <vis/streamline.hpp>

#include "ftle.hpp"
#include <data/raster.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
double length, eps;
size_t nsamples[3];
double abc[3];

void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    hestOptAdd(&hopt, "l",      "length",           airTypeDouble,  1,  1,  &length,        NULL,   "integration length");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  1,  1,  &eps,           NULL,   "integration precision");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis");
    hestOptAdd(&hopt, "abc",    "A B C",            airTypeDouble,  3,  3,  abc,            "1.73205080756888 1.41421356237310 1",  "ABC constants");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Compute flow map and FTLE field in ABC vector field",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct ABC_field {
    ABC_field(double a, double b, double c) : _a(a), _b(b), _c(c),
        _bounds(spurt::vec3(0,0,0), spurt::vec3(2*M_PI, 2*M_PI, 2*M_PI)) {}
        
    bool operator()(double, const spurt::vec3& x, spurt::vec3& f) const {
        f[0] = _a*sin(x[2]) + _c*cos(x[1]);
        f[1] = _b*sin(x[0]) + _a*cos(x[2]);
        f[2] = _c*sin(x[1]) + _b*cos(x[0]);
        return true;
    }
    
    const spurt::bbox3& bounds() const {
        return _bounds;
    }
    
    double _a, _b, _c;
    spurt::bbox3 _bounds;
};


struct integration_monitor {

    enum state { OK = 0,
                 MAX_NB_STEPS,
                 MAX_TIME,
                 STEP_SIZE_UNDERFLOW
               };
               
    integration_monitor() {}
    
    void reset() {}
    
    bool operator()(const spurt::dopri5<spurt::fixed_vector<double, 5ul> >::step&) {
        return true;
    }
};

int main(int argc, const char* argv[])
{
    using namespace spurt;
    
    initialize(argc, argv);
    
    ABC_field wrapper(abc[0], abc[1], abc[2]);
    
    spurt::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    spurt::RasterGrid<3> sampling_grid(res, wrapper.bounds());
    spurt::RasterData<spurt::vec3, 3> flowmaps[2]
        = { spurt::RasterData<spurt::vec3, 3>(sampling_grid),
            spurt::RasterData<spurt::vec3, 3>(sampling_grid)
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
    
    float* __ftle = (float*)calloc(2 * npoints, sizeof(float));
    float* __fmap = (float*)calloc(6 * npoints, sizeof(float));
    std::cerr << "flow map and FTLE values will be exported" << std::endl;
    
    int counter = 0;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        
        for (int n = 0 ; n < npoints ; ++n) {
        
            #pragma omp atomic
            ++counter;
            
            float pct = 100 * (float)counter / (float)npoints;
            std::ostringstream os;
            os << "\r" << counter << " (" << pct << "%)              \r" << std::flush;
            std::cout << os.str();
            
            spurt::ivec3 c = sampling_grid.coord(n);
            spurt::vec3 seed = sampling_grid(c);
            spurt::streamline sl(seed);
            sl.record = false;
            sl.reltol = sl.abstol = eps;
            sl.stepsz = 0;
            integration_monitor nostop;
            try {
                std::ostringstream os;
                spurt::streamline::state state = sl.advance(wrapper, length, nostop);
                flowmaps[0](c) = sl(sl.t_max());
                state = sl.advance(wrapper, -length, nostop);
                flowmaps[1](c) = sl(sl.t_min());
                
                for (int k = 0 ; k < 3 ; ++k) {
                    __fmap[6*n+k] = flowmaps[0](c)[k];
                    __fmap[6*n+3+k] = flowmaps[1](c)[k];
                }
                
            } catch (...) {
                std::ostringstream os;
                os << "exception caught\n";
                std::cerr << os.str();
                continue;
            }
        }
    }
    
    std::cout << "\ntotal computation time for flow map was " << timer.elapsed() << '\n';
    
    std::vector<size_t> size(4);
    std::vector<double> step(4);
    const spurt::vec3& s = sampling_grid.step();
    step[0] = airNaN();
    size[0] = 2;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i+1] = res[i];
        step[i+1] = s[i];
    }
    
    timer.restart();
#pragma openmp parallel for
    for (int n = 0 ; n < npoints ; ++n) {
        try {
            __ftle[2*n] = ftle::ftle(n, flowmaps[0], length);
        } catch (...) {
        }
        try {
            __ftle[2*n+1] = ftle::ftle(n, flowmaps[1], length);
        } catch (...) {
        }
        if (spurt::invalid(__ftle[2*n]) || spurt::invalid(__ftle[2*n+1])) {
            continue;
        }
    }
    std::cout << "total computation time for ftle was " << timer.elapsed() << '\n';
    
    std::ostringstream os;
    os << name_out << "-ftle-T=" << length << ".nrrd";
    spurt::writeNrrd(__ftle, os.str(), nrrdTypeFloat, size, step);
    
    size[0] = 6;
    os.clear();
    os.str("");
    os << name_out << "-flowmap-T=" << length << ".nrrd";
    spurt::writeNrrd(__fmap, os.str(), nrrdTypeFloat, size, step);
    
    return 0;
}


























