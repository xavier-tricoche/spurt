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
double length, h, scale;
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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1,  1,  &name_in,       NULL,       "input file name");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  0,  1,  &name_out,      "none",     "output name");
    hestOptAdd(&hopt, "l",      "length",           airTypeDouble,  1,  1,  &length,        NULL,       "integration length");
    hestOptAdd(&hopt, "scl",    "scale",            airTypeDouble,  0,  1,  &scale,         "1.0e+6",   "scaling factor for vector field");
    hestOptAdd(&hopt, "h",      "step",             airTypeDouble,  1,  1,  &h,             NULL,       "step size");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,       "number of samples per axis");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE in NRRD vector field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T>
struct scaled_field {
    typedef T   field_type;
    
    scaled_field(const field_type& field, double scale = 1)
        : _field(field), _s(scale) {}
        
    bool operator()(const nvis::vec3& x, nvis::vec3& f) const {
        bool valid = _field.get_value(x, f);
        f *= _s;
        return valid;
    }
    
    const field_type& _field;
    double _s;
};

typedef xavier::nrrd_data_traits<Nrrd*>  field_type;
typedef scaled_field<field_type>    rhs_type;

template<typename RHS>
struct euler {

    typedef RHS     rhs_type;
    
    enum state {
        OK = 0,
        LEFT_DOMAIN,
    };
    
    euler(const rhs_type& rhs, double h) : _rhs(rhs), _h(h) {}
    
    state fmap(const nvis::vec3& in, nvis::vec3& out, double length, double& actual_length) const {
        double h = (length < 0 ? -_h : _h);
        actual_length = 0;
        out = in;
        nvis::vec3 f;
        int n = floor(fabs(length) / _h);
        for (int i = 0 ; i < n ; ++i) {
            if (_rhs(out, f)) {
                out += h * f;
                actual_length += fabs(h);
            } else {
                return LEFT_DOMAIN;
            }
        }
        double dt = length - n * h;
        if (_rhs(out, f)) {
            out += dt * f;
            actual_length += dt;
        } else {
            return LEFT_DOMAIN;
        }
        return OK;
    }
    
    const rhs_type& _rhs;
    double _h;
};

int main(int argc, const char* argv[])
{
    using namespace xavier;
    
    initialize(argc, argv);
    
    Nrrd* nin = xavier::nrrd_utils::readNrrd(name_in);
    field_type  vf(nin);
    rhs_type    rhs(vf, scale);
    euler<rhs_type> intg(rhs, h);
    
    nvis::fixed_vector<size_t, 3> res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    xavier::rgrid3d sampling_grid(res, vf.bounds());
    xavier::image3d<nvis::vec3> flowmaps[2] = {
        xavier::image3d<nvis::vec3>(sampling_grid),
        xavier::image3d<nvis::vec3>(sampling_grid)
    };
    
    nvis::timer timer;
    int npoints = sampling_grid.size();
    
    const nvis::vec3 zero(0.);
    for (int i = 0 ; i < npoints ; ++i) {
        nvis::ivec3 c = flowmaps[0].grid().coordinates(i);
        flowmaps[0](c) = flowmaps[1](c) = zero;
    }
    
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    float* __ftle = (float*)calloc(2 * npoints, sizeof(float));
    std::cerr << "FTLE values will be exported" << std::endl;
    
    float* __fmap = (float*)calloc(6 * npoints, sizeof(float));
    float* __leng = (float*)calloc(2 * npoints, sizeof(float));
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        
        for (int n = 0 ; n < npoints ; ++n) {
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
            
            nvis::ivec3 c = sampling_grid.coordinates(n);
            nvis::vec3 seed = sampling_grid(c);
            double l_fwd, l_bwd;
            intg.fmap(seed, flowmaps[0](c), length, l_fwd);
            intg.fmap(seed, flowmaps[1](c), -length, l_bwd);
            // std::ostringstream os;
            // os << '\n' << flowmaps[0](c) << " <---> " << flowmaps[1](c) << '\n';
            // std::cerr << os.str() << std::flush;
            
            for (int i = 0 ; i < 3 ; ++i) {
                __fmap[6*n+i  ] = flowmaps[0](c)[i];
                __fmap[6*n+3+i] = flowmaps[1](c)[i];
                // std::cerr << __fmap[6*n+i] << ", " << __fmap[6*n+3+i] << '\n';
            }
            __leng[2*n  ] = l_fwd;
            __leng[2*n+1] = l_bwd;
        }
    }
    
    std::cout << "\ntotal computation time for flow map was " << timer.elapsed() << '\n';
    
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
    }
    std::cout << "total computation time for ftle was " << timer.elapsed() << '\n';
    
    std::string out;
    if (!strcmp(name_out, "none")) {
        std::string in(name_in);
        out = in.substr(0, in.size()-5);
        std::cerr << "out = " << out << '\n';
    } else {
        out = std::string(name_out);
    }
    
    const nvis::vec3& s = sampling_grid.spacing();
    size_t size[] = {2, res[0], res[1], res[2]};
    double spc[] = {airNaN(), s[0], s[1], s[2]};
    std::ostringstream os;
    os << out << "-ftle-T=" << length << ".nrrd";
    Nrrd* nout = nrrdNew();
    if (nrrdWrap_nva(nout, __ftle, nrrdTypeFloat, 4, size)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    
    os.clear();
    os.str("");
    os << out << "-actual_T.nrrd";
    nout = nrrdNew();
    if (nrrdWrap_nva(nout, __leng, nrrdTypeFloat, 4, size)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    
    size[0] = 6;
    os.clear();
    os.str("");
    os << out << "-fmap-T=" << length << ".nrrd";
    nout = nrrdNew();
    if (nrrdWrap_nva(nout, __fmap, nrrdTypeFloat, 4, size)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    if (nrrdSave(os.str().c_str(), nout, NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    nrrdNuke(nout);
    
    return 0;
}















































