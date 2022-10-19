#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <fstream>
#include <new>

#include <teem/hest.h>
#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <misc/time_helper.hpp>

#include "data/raster.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

char* name_in, *pref;
char* name_out;
double length, eps, t0, T;
int dim;
int mem;
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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1,  1,  &name_in,       NULL,   "list of input files");
    hestOptAdd(&hopt, "p",      "path",             airTypeString,  1,  1,  &pref,          NULL,   "path to input files");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    hestOptAdd(&hopt, "t0",     "init time",        airTypeDouble,  1,  1,  &t0,            NULL,   "seed time");
    hestOptAdd(&hopt, "T",      "integration time", airTypeDouble,  1,  1,  &T,             NULL,   "integration time");
    hestOptAdd(&hopt, "mem",    "memory",           airTypeInt,     0,  1,  &mem,           "1024", "available memory (in MB)");
    hestOptAdd(&hopt, "e",      "epsilon",          airTypeDouble,  1,  1,  &eps,           NULL,   "integration precision");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute flowmap in NRRD vector field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef spurt::nrrd_data_traits<Nrrd*>  field_type;

std::map<float, std::string> file_names;
std::vector<float> _times;

typedef spurt::vec3 point_type;
typedef spurt::vec3 deriv_type;
typedef double scalar_type;

struct Observer {
    Observer(point_type& seed, double& t, double& d) : last_p(seed), last_t(t), distance(d) {}
    void operator()(const point_type& p, double t) {
        distance += spurt::norm(last_p-p);
        last_p = p;
        last_t = t;
    }
    
    point_type& last_p;
    double& last_t;
    double& distance;
};
typedef Observer observer_t;

struct RHS {
    RHS(const field_type& field) : __field(field) {}

    bool operator()(const spurt::vec3& x, spurt::vec3& f, double t) const {
        return __field.get_value(x, f);
    }

    const spurt::bbox3& bounds() const {
        return __field.bounds();
    }

    const field_type&   __field;
};
typedef RHS steady_rhs_type;

inline size_t KB(const size_t s)
{
    return s >> 10;
}

inline size_t MB(const size_t s)
{
    return s >> 20;
}

inline size_t GB(const size_t s)
{
    return s >> 30;
}

inline size_t TB(const size_t s)
{
    return s >> 40;
}

inline size_t PB(const size_t s)
{
    return s >> 50;
}

size_t how_many_fit(size_t one_size)
{
    size_t mb = 1 << 20;
    return mem*mb/one_size;
}

template<typename Tag, typename Obj, class Alloc>
class Repository {
public:
    typedef Tag                                 tag_type;
    typedef Obj                                 object_type;
    typedef Alloc                               alloc_function;
    typedef std::map<tag_type, object_type*>    container_type;
    typedef std::queue<tag_type>                queue_type;

    Repository(size_t max_size=10)
        : __size(max_size) {}

    object_type* get(const tag_type& tag) {
        typename container_type::iterator it = __container.find(tag);
        if (it!=__container.end()) {
            return it->second;
        }

        alloc_function alloc;

#pragma openmp atomic
        object_type* obj = alloc(tag);
#pragma openmp atomic
        __container[tag] = obj;
#pragma openmp atomic
        __queue.push(tag);
        it = __container.find(tag);
        return it->second;
    }

    void load(size_t start, size_t end) {

        alloc_function alloc;

        for (size_t i=start ; i<=end ; ++i) {
#pragma openmp atomic
            object_type* obj = alloc(_times[i]);
#pragma openmp atomic
            __container[_times[i]] = obj;
#pragma openmp atomic
            __queue.push(_times[i]);
        }
    }

    void clear() {
        while (!__queue.empty()) {
            __queue.pop();
        }
        for (typename container_type::iterator i=__container.begin() ;
                i!=__container.end() ; ++i) {
            delete i->second;
            i->second = 0;
        }
        __container.clear();
    }

private:
    container_type      __container;
    queue_type          __queue;
    size_t              __size;
};


struct load_file {
    field_type* operator()(float t) {

        spurt::timer dt;
        std::map<float, std::string>::iterator it = file_names.find(t);
        if (it==file_names.end()) {
            throw std::runtime_error("invalid time step");
        }
        Nrrd* nin = spurt::readNrrd(it->second);

        std::ostringstream os;
        os << it->second << " imported in " << dt.elapsed() << " s." << std::endl;
        std::cerr << os.str();

        return new field_type(nin);
    }
};

typedef Repository<float, field_type, load_file>    repository_type;

class time_dependent_field {
    bool tbounds(float& before, float& after, float t) const {
        if (t<__tmin || t>__tmax) {
            return false;
        }
        size_t rank = std::distance(_times.begin(), std::lower_bound(_times.begin() + __start, _times.begin() + __end+1, t));
        before = _times[rank-1];
        after = _times[rank];
        return true;
    }

public:
    time_dependent_field(const spurt::bbox3& bounds) : __bounds(bounds) {
        __tmin = std::numeric_limits<float>::max();
        __tmax = std::numeric_limits<float>::min();
    }

    void import(size_t start, size_t end) {
        __repo.load(start, end);
        __tmin = _times[start];
        __tmax = _times[end];
        __start = start;
        __end = end;
    }

    void clear() {
        __repo.clear();
    }

    bool operator()(double t, const spurt::vec3& x, spurt::vec3& f) const {
        float prev, next;
        if (!tbounds(prev, next, t)) {
            return false;
        }

        field_type* field0 = __repo.get(prev);
        field_type* field1 = __repo.get(next);

        spurt::vec3 f0, f1;
        if (!field0->get_value(x, f0)) {
            return false;
        }
        field1->get_value(x, f1);
        float u = (t-prev)/(next-prev);
        f = (1.-u)*f0 + u*f1;
        std::cout << "f(" << x << ")=" << f << '\n';
        return true;
    }

    const spurt::bbox3& bounds() const {
        return __bounds;
    }

private:
    mutable repository_type     __repo;
    spurt::bbox3                __bounds;
    float                       __tmin, __tmax;
    size_t                      __start, __end;
};

int main(int argc, const char* argv[])
{
    using namespace spurt;

    initialize(argc, argv);

    // read file names
    std::fstream in(name_in, std::ios::in);
    int nb_files;
    in >> nb_files;
    float _time;
    std::string name;
    for (int i=0 ; i<nb_files ; ++i) {
        in >> _time >> name;
        std::string path = pref;
        path += '/' + name;
        file_names[_time] = path;
        _times.push_back(_time);
        std::cerr << path << " -> " << _time << std::endl;
    }
    in.close();

    Nrrd* nin = spurt::readNrrd(file_names.begin()->second);
    field_type vf(nin);
    steady_rhs_type wrapper(vf);

    size_t field_size = vf.size();
    size_t nb_alloc = how_many_fit(field_size);
    std::cerr << nb_alloc << " time steps fit in\n";

    time_dependent_field field(wrapper.bounds());

    spurt::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    spurt::raster_grid<3> sampling_grid(res, wrapper.bounds());
    int npoints = sampling_grid.size();
    float* flowmap = (float*)calloc(3*npoints, sizeof(float));

    spurt::timer timer;

    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';

    // initialize coordinates
#pragma openmp parallel
    for (int n=0 ; n<npoints ; ++n) {
        spurt::vec3 x = sampling_grid(sampling_grid.coordinates(n));
        flowmap[3*n  ] = x[0];
        flowmap[3*n+1] = x[1];
        flowmap[3*n+2] = x[2];
    }
    std::vector<bool> stopped(npoints, false);

#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif

    float tmin = std::min(t0, t0+T);
    float tmax = std::max(t0, t0+T);
    size_t low = std::distance(_times.begin(), std::lower_bound(_times.begin(), _times.end(), tmin));
    size_t high = std::distance(_times.begin(), std::upper_bound(_times.begin(), _times.end(), tmin));
    if (low > 0) {
        --low;
    }
    if (high < _times.size()-1) {
        ++high;
    }

    std::cerr << "low = " << low << ", high = " << high << std::endl;

    int incr = (T > 0) ? +1 : -1;
    int cur = (T > 0) ? low : high;
    std::vector<int> steps;
    steps.push_back(cur);
    for (int count=1 ; (T>0 && _times[cur]<tmax) || (T<0 && _times[cur]>tmin) ; cur += incr, ++count) {
        if (count % nb_alloc == 0) {
            steps.push_back(cur);
        }
    }
    steps.push_back(cur);
    std::cerr << "integration milestones are:\n";
    for (int i=0 ; i<steps.size() ; ++i) {
        std::cerr << _times[steps[i]] << std::endl;
    }

    double cur_time = t0;
    for (int k=0 ; k<steps.size()-1 ; ++k) {
        field.clear();
        double dT;
        if (T<0) {
            size_t min = std::max(steps[k+1]-1, 0);
            size_t max = std::min((size_t)steps[k]+1, _times.size()-1);
            field.import(min, max);
            dT = std::max((double)(_times[steps[k]] - _times[steps[k+1]]), t0+T - cur_time);
        } else {
            size_t min = std::max(steps[k]-1, 0);
            size_t max = std::min((size_t)steps[k+1]+1, _times.size()-1);
            field.import(steps[k], steps[k+1]);
            dT = std::min((double)(_times[steps[k+1]] - _times[steps[k]]), t0+T - cur_time);
        }
        std::cerr << "current time = " << cur_time << ", dT = " << dT << std::endl;

        size_t nbcomputed;
        spurt::progress_display progress(true);

        progress.start(npoints);

     #pragma omp parallel
        {
        #pragma omp for schedule(dynamic,1)
        for (size_t n = 0 ; n < npoints ; ++n) {
            #if _OPENMP
            const int thread=omp_get_thread_num();
            #else
            const int thread=0;
            #endif

            if (!thread) progress.update(n);

                if (!(n%10)) {
                    float pct = 100.*(float)n/(float)npoints;
                    std::cerr << n << " trajectories (" << pct << "%) computed in "
                              << timer.elapsed()
                              << "s.                 \n"
                              << std::flush;
                }

                if (stopped[n]) {
                    continue;
                }

                spurt::vec3 seed(flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]);
                std::cout << "seed=" << seed << '\n';
                sl.record = false;
                sl.reltol = sl.abstol = eps;
                sl.stepsz = 0;
                try {
                    std::ostringstream os;
                    spurt::streamline::state state = sl.advance(field, cur_time+dT, none);
                    std::cout << "end state = " << state << '\n';
                    spurt::vec3 p = sl(sl.t_max());
                    flowmap[3*n  ] = p[0];
                    flowmap[3*n+1] = p[1];
                    flowmap[3*n+2] = p[2];
                    std::cout << spurt::norm(p-seed) << "\n";
                } 
                catch (std::runtime_error& e) {
                    std::ostringstream os;
                    os << "unable to integrate from " << seed << '\n';
                    os << "error message was:\n" << e.what() << std::endl;
                    std::cerr << os.str();
                    stopped[n] = true;
                    continue;
                }
            }
        }

        cur_time += dT;
    }

    std::cout << "\ntotal computation time for flow map was " << timer.elapsed() << '\n';

    std::vector<size_t> size(4);
    std::vector<double> step(4);
    const spurt::vec3& s = sampling_grid.spacing();
    step[0] = 1 /*airNaN()*/;
    size[0] = 3;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i+1] = res[i];
        step[i+1] = s[i];
    }

    std::ostringstream os;
    os << name_out << "-flowmap-t0=" << t0 << "-T=" << T << ".nrrd";
    spurt::writeNrrdFromContainers(flowmap, os.str(), size, step);

    return 0;
}
