// stdlib
#include <exception>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <sstream>
#include <string>
// nvis
#include <math/bounding_box.hpp>
#include <math/fixed_vector.hpp>
#include <util/timer.hpp>
#include <vis/integral_curve.hpp>
// openmp
#ifdef _OPENMP
#include <omp.h>
#endif
// spurt
#include "denis.hpp"
#include <image/nrrd_wrapper.hpp>
#include <poincare/metric.hpp>

typedef spurt::denis::discrete_holmes   map_type;

// Equation parameters
double      _omega, _T, _gamma, _rho, _g_star, _a;
size_t      _N;

// Computation parameters and variables
nvis::ivec2               _res;
nvis::bbox2               _bounds;
int                       _subres;
std::vector<nvis::vec2>   _raster;
nvis::vec2                _step;
bool                      _forward;

// Algorithm parameters
bool _verbose;

const double TINY = 1.0e-8;

struct uniform_weight {
    double operator()(const nvis::vec2&) const {
        return 1.;
    }
};

template<typename Sampler, typename Value = typename Sampler::value_type, 
         typename Weight = uniform_weight>
inline Value jittered_sampling(const nvis::bbox2& region, int res, 
                               const Sampler& f, const Weight& weight) {
    if (res == 1) return f(0.5*(region.min() + region.max()));
    const nvis::vec2 step = region.size() / nvis::vec2(res, res);
    Value v(0);
    double tot_w = 0;
    for (int n=0 ; n<res*res ; ++n) {
        int i = n%res;
        int j = n/res;
        nvis::vec2 x = region.min() + 
            nvis::vec2(i + drand48(), j + drand48())*step;
        double w = weight(x);
        v += w*f(x);
        tot_w += w;
    }
    v /= tot_w;
    return v;
}

struct tangent_probe {
    typedef nvis::vec2  value_type;
    
    tangent_probe(const map_type& phi) : _phi(phi) {}
    
    nvis::vec2 operator()(const nvis::vec2& x) const {
        nvis::vec2 f = _phi(x, true, false) - x;
        nvis::vec2 b = x - _phi(x, false, false);
        double lf = std::max(nvis::norm(f), TINY);
        double lb = std::max(nvis::norm(b), TINY);
        if (!_forward) {
            double u = lf/(lb + lf);
            return u*b + (1-u)*f;
        }
        else {
            return f;
        }
    }
    
    const map_type& _phi;
};

struct tangent_weight {    
    double operator()(const nvis::vec2& v) const {
        return 1./std::max(nvis::norm(v), TINY);
    }
};

void probe(const nvis::ivec2& p, int subres, const map_type& phi) {    
    nvis::vec2 corner = _bounds.min() + nvis::vec2(p)*_step;
    nvis::bbox2 cell(corner, corner + _step);
    
    tangent_weight weight;
    tangent_probe probe(phi);
    
    int n = p[0] + p[1]*_res[0];
    _raster[n] = jittered_sampling<tangent_probe, nvis::vec2, tangent_weight>
                    (cell, subres, probe, weight);
}

std::string me;
void usage(const std::string& what = "")
{
    if (what != "") {
        std::cerr << "ERROR: " << what << '\n';
    }
    std::cerr
        << "DESCRIPTION: Compute vector field approximation of discrete Holmes\n"
        << "             model of 1D particle column undergoing tapping.\n"
        << '\n'
        << "USAGE:" << me << " [parameters] [options]\n"
        << '\n'
        << "PARAMETERS:\n"
        << " -o | --output             Output file\n" 
        << '\n'
        << "OPTIONS:\n"
        << " -h | --help               Print this information\n"
        << " -r | --rho <float>        Coefficient of restitution\n"
        << " -g | --gamma <float>      Gamma constant\n"
        << " -a | --amplitude <float>  Normalized tap amplitude\n"
        << " -w | --omega <float>      Omega constant\n"
        << " -n | --number <int>       Number of particles in column\n"
        << " -T | --period <float>     Overall tap period\n"
        << " -b | --bounds <float> x2  Velocity bounds\n"
        << " -s | --size <int> x2      Sampling resolution\n"
        << " -j | --jitter <int>       Jittering resolution\n"
        << " -f | --forward            Forward only\n"
        << " -v | --verbose            Turn on verbose mode\n"
        << '\n';

    exit(!what.empty());
}

int main(int argc, char* argv[])
{    
    double _vmin  = -1.;
    double _vmax  =  1.;
    
    _N       = 20;
    _T       = 10.;
    _rho     = 0.8;
    _a       = -1;    // invalid value
    _gamma   = -1;    // invalid value
    _omega   = -1;    // invalid value
    _verbose = false;
    _res[0]  = 1400;
    _res[1]  = 850;
    _subres  = 3;
    _forward = false;
    
    std::string output;
    
    me = argv[0];
    for (int i=0 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-a" || arg == "--amplitude") {
            if (i == argc-1) usage("Missing amplitude value");
            _a = atof(argv[++i]);
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) usage("Missing output file name");
            output = argv[++i];
        }
        else if (arg == "-r" || arg == "--rho") {
            if (i == argc-1) usage("Missing rho value");
            _rho = atof(argv[++i]);
        }
        else if (arg == "-g" || arg == "--gamma") {
            if (i == argc-1) usage("Missing gamma value");
            _gamma = atof(argv[++i]);
        }
        else if (arg == "-w" || arg == "--omega") {
            if (i == argc-1) usage("Missing omega value");
            _omega = atof(argv[++i]);
        }
        else if (arg == "-N" || arg == "-n" || arg == "--particles") {
            if (i == argc-1) usage("Missing number of particles");
            _N = atoi(argv[++i]);
        }
        else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) usage("Missing period value");
            _T = atof(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-2) usage("Missing velocity bounds");
            _vmin = atof(argv[++i]);
            _vmax = atof(argv[++i]);
        }
        else if (arg == "-s" || arg == "--size") {
            if (i >= argc-2) usage("Missing sampling resolution");
            _res[0] = atoi(argv[++i]);
            _res[1] = atoi(argv[++i]);
        }
        else if (arg == "-j" || arg == "--jitter") {
            if (i == argc-1) usage("Missing jittering resolution");
            _subres = atoi(argv[++i]);
        }
        else if (arg == "-f" || arg == "--forward") {
            _forward = true;
        }
        else if (arg == "-v" || arg == "--verbose") {
            _verbose = true;
        }
    }
    
    if (output.empty()) usage("No output filename provided");
    
    _g_star = map_type::g/(double)_N;
    
    // \gamma := 2*\omega^2*a*(1 + \rho)/g
    // a = g*\gamma/(2*\omega^2*(1 + \rho))
    if (_a < 0 && _gamma > 0 && _omega > 0) {
        // compute amplitude as a function of gamma and omega
        _a = _gamma*_g_star/(2*_omega*_omega*(1+_rho)); 
    }
    else if (_gamma < 0 && _a > 0 && _omega > 0) {
        _gamma = 2*_omega*_omega*_a*(1+_rho)/_g_star;
    }
    else if (_omega < 0 && _a > 0 && _gamma > 0) {
    }
    else {
        if (_a < 0) _a = 0.1;
        if (_gamma < 0) _a = 0.1;
    }
    
    // linear section of the phase portrait:
    // \pi <= \theta + v <= \omega*T
    
    map_type phi(_T, _a, _rho, _gamma, _N);
    _omega = phi.omega();
    
    if (_verbose) {
        std::cout << "parameters:\n"
                  << "omega = " << _omega << ", pi/omega = " << M_PI/_omega << '\n'
                  << "gamma = " << _gamma << '\n'
                  << "amp   = " << _a << '\n'
                  << "T     = " << _T << '\n'
                  << "N     = " << _N << '\n';
    }
    
    _bounds.min() = nvis::vec2(0,         _vmin);
    _bounds.max() = nvis::vec2(_omega*_T, _vmax);
    
    if (_verbose) std::cout << "computation bounds are " << _bounds << '\n';
    
    srand48(time(0));
    
    _raster.resize(_res[0]*_res[1]);
    _step = _bounds.size() / (nvis::vec2(_res) - nvis::vec2(1,1));
    
    nvis::timer _timer;
    
#   pragma openmp parallel for 
    for (int n=0 ; n<_raster.size() ; ++n) {
        int i = n%_res[0];
        int j = n/_res[0];
        nvis::ivec2 p(i, j);
        probe(p, _subres, phi);
        // progress display
        int thread_id = 0;
#       if _OPENMP        
        thread_id = omp_get_thread_num();
#       endif
        if (!thread_id) {
            double dt = _timer.elapsed();
            std::cout << "\r" << std::setprecision(3) << std::setw(4) 
            << n << " / " << _raster.size() 
            << " (" << (float)n*100./(float)_raster.size() << "%) in " 
            << dt << " s. (" << n/dt << " Hz)          \r";
        }
    }
    std::cout << "Total computation time: " << _timer.elapsed() << '\n';
    
    double* data = (double*)calloc(_raster.size()*2, sizeof(double));
    for (int i=0 ; i<_raster.size() ; ++i) {
        data[2*i  ] = _raster[i][0];
        data[2*i+1] = _raster[i][1];
    }
    
    double _nan = airNaN();
    if (_nan == 0) _nan = sqrt(-1);
    
    spurt::nrrd_params<double, 3> param;
    param.mins()[0]     = _nan;
    param.mins()[1]     = _bounds.min()[0];
    param.mins()[2]     = _bounds.min()[1];
    param.spacings()[0] = _nan;
    param.spacings()[1] = _step[0];
    param.spacings()[2] = _step[1];
    param.sizes()[0]    = 2;
    param.sizes()[1]    = _res[0];
    param.sizes()[2]    = _res[1];
    
    spurt::writeNrrdFromParams(data, output, param);
    std::cout << output << " has been successfully exported\n";
    
    return 0;
}
