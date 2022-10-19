// stdlib
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
// openmp
#ifdef _OPENMP
#include <omp.h>
#endif
// nvis
#include <math/fixed_vector.hpp>
#include <util/timer.hpp>
// teem
#include <teem/nrrd.h>
#include "denis.hpp"

using namespace xavier::denis;

// constant of gravity
const double _G  = 9.80665;

// Equation parameters
double  _omega, _T, _gamma, _rho, _g_star, _a, _period;
size_t  _N;

// Computation parameters
double _vmin, _vmax;
size_t _niter, _npass, _width, _height;

inline double _modulo(double x, double y=_period)
{
    if (!x) return 0;
    else return x - y*floor(x/y);
}

inline double W(double s)
{
    return (s <= M_PI ? cos(s) : 0 );
}

inline nvis::vec2 phi(const nvis::vec2& x)
{
    const double& theta = x[0];
    const double& v     = x[1];
    double theta_next   = _modulo(theta+v);
    return nvis::vec2(theta_next, _rho*v + _gamma*W(theta_next));
}

std::string me;
void usage(const std::string& what = "",
           bool doExit = true)
{
    if (what != "") {
        std::cerr << "An error occurred: " << what << '\n';
    }
    std::cout
            << "Description: Visualization of discrete dynamical system described\n"
            << "             in Equation 34 of Chaos submission.\n"
            << "Usage  : " << me << " [parameters] [options]\n"
            << "Parameteds:\n"
            << " -o | --output <string>     Ouput image name\n"
            << "Options:\n"
            << " -h | --help                Print this information\n"
            << " -i | --iter <int>          Number of iterations of the map (default: " << _niter << ")\n"
            << " -p | --pass <int>          Number of averaging passes (default: " << _npass << ")\n"
            << " -s | --size <int> x2       Resolution of output image (default: " << _width << ", " << _height << ")\n"
            << " -b | --bounds <float> x2   Bouds for v (default: " << _vmin << ", " << _vmax << ")\n"
            << " Mode parameters:\n"
            << " -a | --amplitude <float>   Tap amplitude (default: " << _a << ")\n"
            << " -g | --gamma <float>       Gamma coefficient (default: " << _gamma << ")\n"
            << " -r | --rho <float>         Coefficient of restitution (default: " << _rho << ")\n"
            << " -n | --number <int>        Number of particles in column (default: " << _N << ")\n"
            << " -T | --period <float>      Time interval between taps (s) (default: " << _T << ")\n"
            << std::endl;
            
    if (doExit) {
        if (!what.empty()) {
            exit(1);
        } else {
            exit(0);
        }
    }
}

typedef nvis::fvec3 color_type;

inline nvis::vec2 i2w(const nvis::ivec2& i)
{
    static const double scale_x = _period/(double)_width;
    static const double scale_y = (_vmax-_vmin)/(double)_height;
    
    return nvis::vec2((double)(i[0] + 0.5)*scale_x,
                      _vmin + (double)(i[1] + 0.5)*scale_y);
}

inline nvis::ivec2 w2i(const nvis::vec2& w)
{
    static const double scale_x = (double)_width/_period;
    static const double scale_y = (double)_height/(_vmax-_vmin);
    
    return nvis::ivec2(floor(w[0]*scale_x),
                       floor((w[1]-_vmin)*scale_y));
}

inline nvis::ivec2 id2i(size_t id)
{
    return nvis::ivec2(id%_width, id/_width);
}

inline size_t i2id(const nvis::ivec2& i)
{
    return i[0] + _width*i[1];
}

inline nvis::vec2 iterate(const nvis::vec2& x,
                          bool fwd,
                          const discrete_holmes& phi)
{
    return phi(x, fwd, false);
}

inline color_type orbit_average(const color_type* input,
                                size_t ic , size_t niter,
                                const discrete_holmes& phi)
{
    color_type acc = input[ic];
    nvis::vec2 x = i2w(id2i(ic));
    int counter = 0;
    for (size_t n=0 ; n<niter ; ++n, ++counter) {
        x = iterate(x, false, phi);
        nvis::vec2 y(_modulo(x[0]), x[1]);
        if (x[1]<_vmin || x[1]>_vmax) break;
        acc += input[i2id(w2i(y))];
    }
    x = i2w(id2i(ic));
    for (size_t n=0 ; n<niter ; ++n, ++counter) {
        x = iterate(x, true, phi);
        nvis::vec2 y(_modulo(x[0]), x[1]);
        if (x[1]<_vmin || x[1]>_vmax) break;
        acc += input[i2id(w2i(y))];
    }
    return acc / (float)(counter+1);
}

void denis_oa(color_type* out, const color_type* in, const discrete_holmes& phi) {
    srand48(12345678);
    
    size_t npixels = _width*_height;
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif

    nvis::timer _timer;
    
#pragma openmp parallel for
    for (size_t i=0 ; i<npixels ; ++i) {
        out[i] = orbit_average(in, i, _niter, phi);
        
        int thread_id = 0;
#if _OPENMP        
        thread_id = omp_get_thread_num();
#endif
        if (!thread_id) {
            std::cout << "\r" << std::setprecision(3) << std::setw(4) 
            << i << " / " << npixels 
            << " (" << (float)i*100./(float)npixels << "%) in " 
            << _timer.elapsed() << " s.\r";
        }
    }
    
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    me = argv[0];
  
    _vmin   = -1;
    _vmax   = 1.;
    _N      = 20;
    _width  = 1400;
    _height = 512;
    _T      = 10;
    _a      = 0.1;
    _rho    = 0.8;
    _gamma  = 0.1;
    _niter  = 500;
    _npass  = 3;
    
    std::string _output = "undefined!";
    
    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            usage();
        } else if (arg == "-a" || arg == "--amplitude") {
            if (i == argc-1) {
                usage("missing amplitude value");
            }
            _a = atof(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                usage("missing output name");
            }
            _output = argv[++i];
        } else if (arg == "-i" || arg == "--iterations") {
            if (i == argc-1) {
                usage("missing number of iterations");
            }
            _niter = atoi(argv[++i]);      
        } else if (arg == "-p" || arg == "--pass") {
            if (i == argc-1) {
                usage("missing number of passes");
            }
            _npass = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--size") {
            if (i >= argc-2) {
                usage("missing image resolution");
            }
            _width = atoi(argv[++i]);
            _height = atoi(argv[++i]);
        } else if (arg == "-n" || arg == "--number") {
            if (i == argc-1) {
                usage("missing number of particles");
            }
            _N = atoi(argv[++i]);
        } else if (arg == "-r" || arg == "--rho") {
            if (i == argc-1) {
                usage("missing rho value");
            }
            _rho = atof(argv[++i]);
        } else if (arg == "-g" || arg == "--gamma") {
            if (i == argc-1) {
                usage("missing gamma value");
            }
            _gamma = atof(argv[++i]);
        } else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) {
                usage("missing period value");
            }
            _T = atof(argv[++i]);
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-2) {
                usage("missing v boundary value(s)");
            }
            _vmin = atof(argv[++i]);
            _vmax = atof(argv[++i]);
        } else {
            usage(std::string("unrecognized input directive: ") + arg);
        }       
    }
    
    if (_output == "undefined!") {
        usage("missing output name");
    } 
    
      
    
    discrete_holmes phi(_T, _a, _rho, _gamma, _N);
    _omega = phi.omega();
    
    _period = _omega*_T;
    
    double H = _vmax - _vmin;
    double ratio = H/(_omega*_T);
    double newheight = ratio*_width;
    double newwidth = _height/ratio;
    if (newheight <= _height) _height = newheight;
    else if (newwidth <= _width) _width = newwidth;
    else std::cerr << "unable to resize window to fit!\n";
    
    size_t npixels = _width*_height;
    
    color_type *buf1 = (color_type*)calloc(npixels, sizeof(color_type));  
    color_type *buf2 = (color_type*)calloc(npixels, sizeof(color_type));  
    for (size_t i=0 ; i<npixels ; ++i) {
        buf1[i][0] = drand48();
        buf1[i][1] = drand48();
        buf1[i][2] = drand48();
    }
    
    for (int i=0 ; i<_npass ; ++i) {
        denis_oa(buf2, buf1, phi);
        color_type* tmp = buf1;
        buf1 = buf2;
        buf2 = tmp;
    }
    delete[] buf2;
    
    Nrrd *nout = nrrdNew();
    size_t size[] = {3, _width, _height};
    if (nrrdWrap_nva(nout, (float*)buf1, nrrdTypeFloat, 3, size)) {
        std::cerr << biffGetDone(NRRD) << '\n';
    }

    double spc[] = {AIR_NAN, _period/(double)(_width-1), (_vmax-_vmin)/(double)(_height-1)};
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
    int center[] = {nrrdCenterUnknown, nrrdCenterNode, nrrdCenterNode};
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, center);
    if (nrrdSave(_output.c_str(), nout, NULL)) {
        std::cerr << biffGetDone(NRRD) << '\n';
    }
    
    delete[] buf1;
    
    return 0;
}
