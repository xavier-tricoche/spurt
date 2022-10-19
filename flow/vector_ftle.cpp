#include <iostream>

#include <teem/hest.h>
#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <misc/time_helper.hpp>
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

typedef spurt::nrrd_data_traits<Nrrd*>  field_type;

struct field_wrapper {
    field_wrapper(const field_type& field) : m_field(field) {}
    
    bool operator()(double, const spurt::vec3& x, spurt::vec3& f) const {
        return m_field.get_value(x, f);
    }
    
    const spurt::bbox3& bounds() const {
        return m_field.bounds();
    }
    
    const field_type&   m_field;
};


struct integration_monitor {

    enum state { OK = 0,
                 MAX_NB_STEPS,
                 MAX_TIME,
                 STEP_SIZE_UNDERFLOW,
                 OUT_OF_BOUNDS,
               };
               
    typedef spurt::streamline::int_step  step_type;
    
    integration_monitor(size_t max_steps, double min_step_size, double max_time,
                        const spurt::bbox3& bounds)
        : m_counter(0), m_length(0), m_max_nb_steps(max_steps), m_max_time(max_time),
          m_min_step_size(min_step_size), m_bounds(bounds), m_state(OK) {}
          
    void reset() {
        m_counter = 0;
        m_length = 0;
        m_state = OK;
        m_last_pos = spurt::vec3(0);
    }
    
    bool operator()(const step_type& step) {
    
        std::cerr << "check\n";
        
        ++m_counter;
        double dx = spurt::norm(spurt::subv<0, 3>(step.y1()) - spurt::subv<0, 3>(step.y0()));
        m_length += dx;
        if (!m_bounds.inside(spurt::subv<0, 3>(step.y1()))) {
            // binary search for exit point
            double tmin = step.t0();
            double tmax = step.t1();
            while (tmax - tmin > 1.0e-6) {
                double t = 0.5 * (tmin + tmax);
                if (!m_bounds.inside(spurt::subv<0, 3>(step.y(t)))) {
                    tmax = t;
                } else {
                    tmin = t;
                }
            }
            m_state = OUT_OF_BOUNDS;
            m_last_pos = spurt::subv<0, 3>(step.y(0.5 * (tmin + tmax)));
            return true;
        } else if (m_counter >= m_max_nb_steps) {
            m_state = MAX_NB_STEPS;
            m_last_pos = spurt::subv<0, 3>(step.y1());
            std::cerr << "too many steps\n";
            return true;
        } else if (step.t1() > m_max_time) {
            if (step.t0() > m_max_time) {
                std::cerr << "invalid step: required length already exceeded!";
                throw std::runtime_error("invalid step");
            }
            m_state = MAX_TIME;
            m_last_pos = spurt::subv<0, 3>(step.y(m_max_time));
            return true;
        } else if (dx < m_min_step_size) {
            m_state = STEP_SIZE_UNDERFLOW;
            m_last_pos = spurt::subv<0, 3>(step.y1());
            std::cerr << "step size underflow\n";
            return true;
        } else {
            m_last_pos = spurt::subv<0, 3>(step.y1());
            return false;
        }
    }
    
    size_t       m_counter, m_max_nb_steps;
    double       m_length, m_max_time, m_min_step_size;
    state        m_state;
    spurt::bbox3 m_bounds;
    spurt::vec3  m_last_pos;
};

std::ostream& operator<<(std::ostream& os, const integration_monitor& monitor)
{
    switch (monitor.m_state) {
        case integration_monitor::OK:
            os << "successful << (" << monitor.m_counter << ", " << monitor.m_length << ")";
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
    using namespace spurt;
    
    initialize(argc, argv);
    
    Nrrd* nin = spurt::readNrrd(name_in);
    field_type vf(nin);
    field_wrapper wrapper(vf);
    
    spurt::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    spurt::raster_grid<3> sampling_grid(res, wrapper.bounds());
    spurt::image<spurt::vec3, 3> flowmaps[2]
        = { spurt::image<spurt::vec3, 3>(sampling_grid),
            spurt::image<spurt::vec3, 3>(sampling_grid)
          };
          
    spurt::timer timer;
    int npoints = sampling_grid.size();
    
    const spurt::vec3 zero(0.);
    for (int i = 0 ; i < npoints ; ++i) {
        spurt::ivec3 c = flowmaps[0].grid()(i);
        flowmaps[0](c) = flowmaps[1](c) = zero;
    }
    
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    float* m_ftle = (float*)calloc(2 * npoints, sizeof(float));
    float* m_fmap = (float*)calloc(6 * npoints, sizeof(float));
    float* m_length = (float*)calloc(2 * npoints, sizeof(float));
    int* m_states = (int*)calloc(4 * npoints, sizeof(int));
    std::cerr << "FTLE values will be exported" << std::endl;
    
    std::fill(m_states, m_states + 4*npoints, -1);
    
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
            
            spurt::ivec3 c = sampling_grid(n);
            spurt::vec3 seed = sampling_grid(c);
            spurt::streamline sl(seed);
            sl.record = false;
            sl.reltol = sl.abstol = eps;
            sl.stepsz = 0;
            integration_monitor stopc(2000, 1.0e-6, length, wrapper.bounds());
            try {
                std::ostringstream os;
                spurt::streamline::state state = sl.advance(wrapper, length, stopc);
                // flowmaps[0](c) = stopc.m_last_pos;
                flowmaps[0](c) = sl(sl.t_max());
                // if (!wrapper.bounds().inside(flowmaps[0](c))) {
                //  os << "Warning, " << flowmaps[0](c) << " is out of bounds\n";
                // }
                m_length[2*n] = sl.a_max();
                m_states[4*n] = state;
                m_states[4*n+2] = stopc.m_state;
                double l = stopc.m_length;
                stopc.reset();
                state = sl.advance(wrapper, -length, stopc);
                // flowmaps[1](c) = stopc.m_last_pos;
                flowmaps[1](c) = sl(sl.t_min());
                // if (!wrapper.bounds().inside(flowmaps[1](c))) {
                //  os << "Warning, " << flowmaps[1](c) << " is out of bounds\n";
                // }
                m_length[2*n+1] = sl.a_min();
                m_states[4*n+1] = state;
                m_states[4*n+3] = stopc.m_state;
                l = std::max(l, stopc.m_length);
                // os << "max length was " << l << std::endl;
                // std::cerr << os.str() << std::flush;
                
                for (int k = 0 ; k < 3 ; ++k) {
                    m_fmap[6*n+k] = flowmaps[0](c)[k];
                    m_fmap[6*n+3+k] = flowmaps[1](c)[k];
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
    const spurt::vec3& s = sampling_grid.spacing();
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
            m_ftle[2*n] = ftle::ftle(n, flowmaps[0], length);
        } catch (...) {
        }
        try {
            m_ftle[2*n+1] = ftle::ftle(n, flowmaps[1], length);
        } catch (...) {
        }
        if (spurt::invalid(m_ftle[2*n]) || spurt::invalid(m_ftle[2*n+1])) {
            // std::cerr << "skipping invalid value #" << n << std::endl;
            continue;
        }
    }
    std::cout << "total computation time for ftle was " << timer.elapsed() << '\n';
    
    std::ostringstream os;
    os << name_out << "-ftle-T=" << length << ".nrrd";
    spurt::writeNrrdFromContainers(m_ftle, os.str(), /*nrrdTypeFloat,*/ size, step);
    
    os.clear();
    os.str("");
    os << name_out << "-length-T=" << length << ".nrrd";
    spurt::writeNrrdFromContainers(m_length, os.str(), /*nrrdTypeFloat,*/ size, step);
    
    os.clear();
    os.str("");
    os << name_out << "-states-T=" << length << ".nrrd";
    size[0] = 4;
    spurt::writeNrrdFromContainers(m_states, os.str(), /*nrrdTypeInt,*/ size, step);
    
    size[0] = 6;
    os.clear();
    os.str("");
    os << name_out << "-flowmap-T=" << length << ".nrrd";
    spurt::writeNrrdFromContainers(m_fmap, os.str(), /*nrrdTypeFloat,*/ size, step);
    
    return 0;
}


























