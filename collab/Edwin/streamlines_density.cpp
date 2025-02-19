#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include <vis/streamline.hpp>
#include <boost/shared_ptr.hpp>
#include <data/field_wrapper.hpp>
#include "Garcia_vis_helper.hpp"
#include "math/inverse_transform.hpp"

inline void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

int     dataset;
char*    out;
int     nblines;
float   length;
float   center[3];
float   radius[3];
float   eps;
float   h;
float   lmax;
float   min_step;
float   col_gamma;
int     discretization;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;

    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "id",     "dataset ID",       airTypeInt,     1, 1, &dataset,             NULL,       "dataset ID: 0: textured, 1: untextured, 2: MC_r00b09, 3: MC_r06b09, 4: MC_r10b09");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &out,                 NULL,       "output file");
    hestOptAdd(&hopt, "n",      "nb seeds",         airTypeInt,     0, 1, &nblines,             "100",      "number of lines");
    hestOptAdd(&hopt, "h",      "step size",        airTypeFloat,   0, 1, &h,                   "1.0e-3",   "integration precision");
    hestOptAdd(&hopt, "max",    "max step",         airTypeFloat,   0, 1, &lmax,                "0.5",      "max integration step length");
    hestOptAdd(&hopt, "min",    "min step",         airTypeFloat,   0, 1, &min_step,            "1.0e-8",   "step size underflow threshold");
    hestOptAdd(&hopt, "d",      "discretization",   airTypeInt,     0, 1, &discretization,      "20",       "number of lines");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &col_gamma,           "1.",       "gamma factor for color scale");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute streamline density per grain for bottom face seeding",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T>
inline T sign(const T& t)
{
    return (t < 0 ? -1 : 1);
}

struct check {

    check(size_t max_count, const nvis::bbox3& box)
        : _counter(0), _max_count(max_count), _bounds(box) {}

    void reset() {
        _counter = 0;
    }

    bool operator()(const nvis::streamline::int_step& step) {
        // std::cerr << _counter << nvis::subv<0, 3>(step.y1()) << '\n';

        if (++_counter >= _max_count) {
            return true;
        } else if (!_bounds.inside(nvis::subv<0, 3>(step.y1()))) {
            return true;
        }

        // std::cerr << "continue\n";
        return false;
    }

    size_t _counter, _max_count;
    nvis::bbox3 _bounds;
};

struct i2x {
    i2x(const Nrrd* nrrd) {
        for (int i = 0 ; i < 3 ; ++i) {
            step[i] = nrrd->axis[nrrd->dim-3+i].spacing;
            size[i] = nrrd->axis[nrrd->dim-3+i].size;
        }
        bounds = spurt::nrrd_utils::get_bounds<3>(nrrd);
    }

    nvis::vec3 operator()(int id) const {
        int i = id % size[0];
        id /= size[0];
        int j = id % size[1];
        int k = id / size[1];
        return bounds.min() + nvis::vec3(i, j, k)*step;
    }

    nvis::vec3  step;
    nvis::ivec3 size;
    nvis::bbox3 bounds;
};

template<typename T>
struct right_hand_side {
    typedef T   field_type;

    right_hand_side(const field_type& field)
        : _field(field) {}

    bool operator()(const nvis::vec3& x, nvis::vec3& f) const {
        return _field.get_value(x, f);
    }

    const field_type& _field;
};

template<typename RHS>
struct euler {

    typedef RHS     rhs_type;

    enum state {
        OK = 0,
        LEFT_DOMAIN,
        STOPPED,
    };

    euler(const rhs_type& rhs, double h, double lmax, double eps = 1.0e-6) : _rhs(rhs), _h(h), _lmax(lmax), _eps(eps) {}

    state advance(const nvis::vec3& in, std::list<nvis::vec3>& out, double length, double& actual_length) const {
        double h = (length < 0 ? -_h : _h);
        actual_length = 0;
        nvis::vec3 f, x = in;
        while (actual_length < fabs(length)) {
            // std::cerr << "at " << x << ", h = " << h << '\n';
            if (_rhs(x, f)) {
                out.push_back(x);
                if (nvis::norm(f) < _eps) {
                    return STOPPED;
                } else if (fabs(h)*nvis::norm(f) > _lmax) {
                    h = sign(length) * 0.5 * _lmax / nvis::norm(f);
                    x += h * f;
                } else {
                    x += h * f;
                    if (fabs(h)*nvis::norm(f) < 10.*_lmax) {
                        h *= 10;
                    }
                }
                actual_length += fabs(h);
            } else {
                return LEFT_DOMAIN;
            }
        }
        return OK;
    }

    const rhs_type& _rhs;
    double _h, _lmax, _eps;
};

typedef spurt::nrrd_data_traits<Nrrd*>      field_type;
typedef right_hand_side<field_type>     rhs_type;
typedef euler<rhs_type>                 euler_type;
typedef std::pair<nvis::vec3, double>   curve_point;
typedef std::list<curve_point>          curve_type;
typedef nvis::vec3                      value_type;
typedef spurt::nrrd_data_traits<Nrrd*>      nrrd_data_traits;


struct point_location {
    point_location(const Nrrd* nrrd) {
        nvis::bbox3 b = spurt::nrrd_utils::get_bounds<3>(nrrd);
        min = b.min();
        step = spurt::nrrd_utils::step<3>(nrrd);
        for (int i = 0 ; i < 3 ; ++i) {
            size[i] = nrrd->axis[nrrd->dim-3+i].size;
        }
    }

    int index(const nvis::vec3& x) const {
        nvis::vec3 y = (x - min) / step;
        nvis::ivec3 n;
        for (int i = 0 ; i < 3 ; ++i) {
            n[i] = round(y[i]);
        }
        return n[0] + size[0]*(n[1] + size[1]*n[2]);
    }

    nvis::vec3  min, step;
    nvis::ivec3 size;
};

int main(int argc, char* argv[])
{
    initialize(argc, argv);

    using namespace Garcia_vis_helper;

    set_paths();

    dataset_info* __info;
    switch (dataset) {
        case 0:
            __info = &textured_info;
            break;
        case 1:
            __info = &untextured_info;
            break;
        case 2:
            __info = &mc_info_00;
            break;
        case 3:
            __info = &mc_info_06;
            break;
        case 4:
            __info = &mc_info_10_09;
            break;
        case 5:
            __info = &mc_info_10_05;
            break;
        default:
            std::cerr << "unknown dataset\n";
            return 1;
    }
    const dataset_info& info = *__info;

    std::cerr << "base dir = " << info.base_dir << '\n';

    std::string name = info.base_dir + "Efield.nrrd";
    Nrrd* nin_vec = spurt::nrrd_utils::readNrrd(name);
    if (nin_vec->dim != 4) {
        throw;
    }
    field_type      vf(nin_vec);
    rhs_type        rhs(vf);
    euler_type      integrator(rhs, h, lmax);

    std::string base(out);
    std::ostringstream os;

    name = info.microstruct;
    Nrrd* nin_tag = spurt::nrrd_utils::readNrrd(name);
    point_location pl(nin_tag);
    int* tags = (int*)nin_tag->data;
    int nb_voxels = pl.size[0] * pl.size[1] * pl.size[2];
    std::map<int, int> grain_size;
    for (int i = 0 ; i < nb_voxels ; ++i) {
        int id = tags[i];
        if (id < 0) {
            continue;
        } else if (grain_size.find(id) == grain_size.end()) {
            grain_size[id] = 1;
        } else {
            ++grain_size[id];
        }
    }

    nvis::bbox3 bounds = spurt::nrrd_utils::get_bounds<3>(nin_vec);
    nvis::vec3 min = bounds.min();
    nvis::vec3 diameter = bounds.size();

    std::map<int, int> hit_counter;

    srand48(time(0));
    double intg_time = 0;
    float length = std::numeric_limits<double>::max();

    int nx, ny, nz;
    nx = ny = nz = round(pow(nblines, 1./3.));
    nblines = nx * ny * nz;
    double dx = diameter[0] / (double)nx;
    double dy = diameter[1] / (double)ny;
    double dz = diameter[2] / (double)nz;

    for (int i = 0 ; i < nblines ; ++i) {

        int _i = i % nx;
        int _j = (i / nx) % ny;
        int _k = i / (nx*ny);
        double _x = min[0] + ((double)_i + 0.5 + drand48()) * dx;
        double _y = min[1] + ((double)_j + 0.5 + drand48()) * dy;
        double _z = min[2] + ((double)_k + 0.5 + drand48()) * dz;

        // nvis::vec3 seed = min + nvis::vec3(drand48(), drand48(), 0.001) * diameter;
        nvis::vec3 seed(_x, _y, _z);
        std::cerr << "seeding streamline #" << i << "/" << nblines << " at " << seed << '\n';

        std::list<nvis::vec3>   line, aux;
        double                  actual_l;

        integrator.advance(seed, line, length, actual_l);
        integrator.advance(seed, aux, -length, actual_l);
        if (aux.size()) {
            aux.pop_front();
        }
        line.insert(line.begin(), aux.rbegin(), aux.rend());

        std::cerr << line.size() << " points in line\n";

        int last_id = -5;
        for (std::list<nvis::vec3>::const_iterator it = line.begin() ; it != line.end() ; ++it) {
            nvis::vec3 x = *it;
            int id = pl.index(x);
            if (id < 0 || id == last_id) {
                continue;
            }
            last_id = id;
            int tag = tags[id];
            if (hit_counter.find(tag) == hit_counter.end()) {
                hit_counter[tag] = 1;
            } else {
                ++hit_counter[tag];
            }
        }
    }

    // normalize results
    std::map<int, float> freq;
    for (std::map<int, int>::const_iterator it = hit_counter.begin() ;
            it != hit_counter.end() ; ++it) {
        int id = it->first;
        float fq = (float)it->second / (float)grain_size[id] / (float)nblines;
        freq[id] = fq;
    }

    float* __freq_in_place = (float*)calloc(nb_voxels, sizeof(float));
    for (int i = 0 ; i < nb_voxels ; ++i) {
        int id = tags[i];
        if (id < 0) {
            continue;
        }
        float f = freq[id];
        __freq_in_place[i] = f;
    }
    size_t __size[] = {(size_t)pl.size[0], (size_t)pl.size[1], (size_t)pl.size[2]};
    Nrrd* __nout = nrrdNew();
    os << base << "-n=" << nblines << "-inplace.nrrd";
    nrrdWrap_nva(__nout, __freq_in_place, nrrdTypeFloat, 3, __size);
    nrrdSave(os.str().c_str(), __nout, NULL);
    if (__nout) {
        nrrdNuke(__nout);
    }

    float* __freq = (float*)calloc(2 * freq.size(), sizeof(float));
    int counter = 0;
    for (std::map<int, float>::const_iterator it = freq.begin() ; it != freq.end() ; ++it) {
        __freq[counter++] = it->first;
        __freq[counter++] = it->second;
    }
    Nrrd* nout = nrrdNew();
    size_t size[] = {2, freq.size()};
    nrrdWrap_nva(nout, __freq, nrrdTypeFloat, 2, size);
    os.clear();
    os.str("");
    os << base << "-n=" << nblines << ".nrrd";
    nrrdSave(os.str().c_str(), nout, NULL);

    if (nin_vec) {
        nrrdNuke(nin_vec);
    }
    if (nin_tag) {
        nrrdNuke(nin_tag);
    }
    if (nout) {
        nrrdNuke(nout);
    }

    return 0;
}
