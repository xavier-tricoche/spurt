#include <iostream>

#include <teem/hest.h>
#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/wall_timer.hpp>
#include <vis/streamline.hpp>

#include "ftle.hpp"
#include "data/raster.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in;
char* name_out;
double length, eps;
size_t nsamples[3];

void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1,  1,  &name_in,       NULL,   "input file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    hestOptAdd(&hopt, "l",      "length",           airTypeDouble,  1,  1,  &length,        NULL,   "integration length");
    // hestOptAdd(&hopt, "l",       "length",           airTypeDouble,  1,  1,  &length,        NULL,   "integration length for flow map computation");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  1,  1,  &eps,           NULL,   "integration precision");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE in NRRD vector field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef xavier::nrrd_data_traits<Nrrd*>  field_type;

struct field_wrapper {
    field_wrapper(const field_type& field) : __field(field) {}
    
    bool operator()(double, const nvis::vec3& x, nvis::vec3& f) const {
        return __field.get_value(x, f);
    }
    
    const nvis::bbox3& bounds() const {
        return __field.bounds();
    }
    
    const field_type&   __field;
};


struct integration_monitor {

    enum state { OK = 0,
                 MAX_NB_STEPS,
                 MAX_TIME,
                 STEP_SIZE_UNDERFLOW,
                 OUT_OF_BOUNDS,
               };
               
    typedef nvis::streamline::int_step  step_type;
    
    integration_monitor(size_t max_steps, double min_step_size, double max_time,
                        const nvis::bbox3& bounds)
        : __counter(0), __length(0), __max_nb_steps(max_steps), __max_time(max_time),
          __min_step_size(min_step_size), __bounds(bounds), __state(OK) {}
          
    void reset() {
        __counter = 0;
        __length = 0;
        __state = OK;
        __last_pos = nvis::vec3(0);
    }
    
    bool operator()(const step_type& step) {
    
        std::cerr << "check\n";
        
        ++__counter;
        double dx = nvis::norm(nvis::subv<0, 3>(step.y1()) - nvis::subv<0, 3>(step.y0()));
        __length += dx;
        if (!__bounds.inside(nvis::subv<0, 3>(step.y1()))) {
            // binary search for exit point
            double tmin = step.t0();
            double tmax = step.t1();
            while (tmax - tmin > 1.0e-6) {
                double t = 0.5 * (tmin + tmax);
                if (!__bounds.inside(nvis::subv<0, 3>(step.y(t)))) {
                    tmax = t;
                } else {
                    tmin = t;
                }
            }
            __state = OUT_OF_BOUNDS;
            __last_pos = nvis::subv<0, 3>(step.y(0.5 * (tmin + tmax)));
            return true;
        } else if (__counter >= __max_nb_steps) {
            __state = MAX_NB_STEPS;
            __last_pos = nvis::subv<0, 3>(step.y1());
            std::cerr << "too many steps\n";
            return true;
        } else if (step.t1() > __max_time) {
            if (step.t0() > __max_time) {
                std::cerr << "invalid step: required length already exceeded!";
                throw std::runtime_error("invalid step");
            }
            __state = MAX_TIME;
            __last_pos = nvis::subv<0, 3>(step.y(__max_time));
            return true;
        } else if (dx < __min_step_size) {
            __state = STEP_SIZE_UNDERFLOW;
            __last_pos = nvis::subv<0, 3>(step.y1());
            std::cerr << "step size underflow\n";
            return true;
        } else {
            __last_pos = nvis::subv<0, 3>(step.y1());
            return false;
        }
    }
    
    size_t      __counter, __max_nb_steps;
    double      __length, __max_time, __min_step_size;
    state       __state;
    nvis::bbox3 __bounds;
    nvis::vec3  __last_pos;
};

std::ostream& operator<<(std::ostream& os, const integration_monitor& monitor)
{
    switch (monitor.__state) {
        case integration_monitor::OK:
            os << "successful << (" << monitor.__counter << ", " << monitor.__length << ")";
            break;
        case integration_monitor::MAX_NB_STEPS:
            os << "max nb steps reached";
            break;
        case integration_monitor::MAX_TIME:
            os << "max time reached";
            break;
        case integration_monitor::STEP_SIZE_UNDERFLOW:
            os << "step size underflow";
            break;
        case integration_monitor::OUT_OF_BOUNDS:
            os << "stepped out of bounds";
            break;
        default:
            os << "invalid state";
    }
    return os;
}

int main(int argc, const char* argv[])
{
    using namespace xavier;
    
    initialize(argc, argv);
    
    Nrrd* nin = xavier::nrrd_utils::readNrrd(name_in);
    field_type vf(nin);
    field_wrapper wrapper(vf);
    
    nvis::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    xavier::raster_grid<3> sampling_grid(res, wrapper.bounds());
    xavier::image<nvis::vec3, 3> flowmaps[2]
        = { xavier::image<nvis::vec3, 3>(sampling_grid),
            xavier::image<nvis::vec3, 3>(sampling_grid)
          };
          
    nvis::timer timer;
    int npoints = sampling_grid.size();
    
    const nvis::vec3 zero(0.);
    for (int i = 0 ; i < npoints ; ++i) {
        nvis::ivec3 c = flowmaps[0].grid()(i);
        flowmaps[0](c) = flowmaps[1](c) = zero;
    }
    
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    float* __ftle = (float*)calloc(2 * npoints, sizeof(float));
    float* __fmap = (float*)calloc(6 * npoints, sizeof(float));
    float* __length = (float*)calloc(2 * npoints, sizeof(float));
    int* __states = (int*)calloc(4 * npoints, sizeof(int));
    std::cerr << "FTLE values will be exported" << std::endl;
    
    std::fill(__states, __states + 4*npoints, -1);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        
        for (int n = 0 ; n < npoints ; ++n) {
            std::cerr << n << '\n';
            
            int pct = 100 * n / npoints;
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
            
            nvis::ivec3 c = sampling_grid(n);
            nvis::vec3 seed = sampling_grid(c);
            nvis::streamline sl(seed);
            sl.record = false;
            sl.reltol = sl.abstol = eps;
            sl.stepsz = 0;
            integration_monitor stopc(2000, 1.0e-6, length, wrapper.bounds());
            try {
                std::ostringstream os;
                nvis::streamline::state state = sl.advance(wrapper, length, stopc);
                // flowmaps[0](c) = stopc.__last_pos;
                flowmaps[0](c) = sl(sl.t_max());
                // if (!wrapper.bounds().inside(flowmaps[0](c))) {
                //  os << "Warning, " << flowmaps[0](c) << " is out of bounds\n";
                // }
                __length[2*n] = sl.a_max();
                __states[4*n] = state;
                __states[4*n+2] = stopc.__state;
                double l = stopc.__length;
                stopc.reset();
                state = sl.advance(wrapper, -length, stopc);
                // flowmaps[1](c) = stopc.__last_pos;
                flowmaps[1](c) = sl(sl.t_min());
                // if (!wrapper.bounds().inside(flowmaps[1](c))) {
                //  os << "Warning, " << flowmaps[1](c) << " is out of bounds\n";
                // }
                __length[2*n+1] = sl.a_min();
                __states[4*n+1] = state;
                __states[4*n+3] = stopc.__state;
                l = std::max(l, stopc.__length);
                // os << "max length was " << l << std::endl;
                // std::cerr << os.str() << std::flush;
                
                for (int k = 0 ; k < 3 ; ++k) {
                    __fmap[6*n+k] = flowmaps[0](c)[k];
                    __fmap[6*n+3+k] = flowmaps[1](c)[k];
                }
                
            } catch (...) {
                std::ostringstream os;
                os << "exception caught: stop state is " << stopc << '\n';
                std::cerr << os.str();
                continue;
            }
        }
    }
    
    std::cout << "\ntotal computation time for flow map was " << timer.elapsed() << '\n';
    
    std::vector<size_t> size(4);
    std::vector<double> step(4);
    const nvis::vec3& s = sampling_grid.spacing();
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
        if (xavier::nrrd_utils::invalid(__ftle[2*n]) || xavier::nrrd_utils::invalid(__ftle[2*n+1])) {
            // std::cerr << "skipping invalid value #" << n << std::endl;
            continue;
        }
    }
    std::cout << "total computation time for ftle was " << timer.elapsed() << '\n';
    
    std::ostringstream os;
    os << name_out << "-ftle-T=" << length << ".nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(__ftle, os.str(), /*nrrdTypeFloat,*/ size, step);
    
    os.clear();
    os.str("");
    os << name_out << "-length-T=" << length << ".nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(__length, os.str(), /*nrrdTypeFloat,*/ size, step);
    
    os.clear();
    os.str("");
    os << name_out << "-states-T=" << length << ".nrrd";
    size[0] = 4;
    xavier::nrrd_utils::writeNrrdFromContainers(__states, os.str(), /*nrrdTypeInt,*/ size, step);
    
    size[0] = 6;
    os.clear();
    os.str("");
    os << name_out << "-flowmap-T=" << length << ".nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(__fmap, os.str(), /*nrrdTypeFloat,*/ size, step);
    
    return 0;
}


























