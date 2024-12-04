// stdlib
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
// openmp
#ifdef _OPENMP
#include <omp.h>
#endif
// nvis
#include <math/fixed_vector.hpp>
#include <util/timer.hpp>
// teem
#include <teem/nrrd.h>

// constant of gravity
const double _G  = 9.80665;

// Equation parameters
double  _omega, _T, _gamma, _rho, _g_star, _a, _period;
size_t  _N;

// Computation parameters
double _vmin, _vmax;
size_t _niter, _nseeds, _width, _height;

inline double _modulo(double x, double y=_period)
{
    return x - y*floor(x/y);
}

inline double W(double s)
{
    return (s <= M_PI ? cos(s) : 0 );
}

inline nvis::vec2 phi(const nvis::vec2& x)
{
    const double& theta = x[0];
    const double& v     = x[1];
    double theta_next = _modulo(theta+v);
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
            << " -n | --nseeds <int>        Number of seeds (default: " << _nseeds << ")\n"
            << " -s | --size <int> x2       Resolution of output image (default: " << _width << ", " << _height << ")\n"
            << " -b | --bounds <float> x2   Bouds for v (default: " << _vmin << ", " << _vmax << ")\n"
            << " Mode parameters:\n"
            << " -a | --amplitude <float>   Tap amplitude (default: " << _a << ")\n"
            << " -g | --gamma <float>       Gamma coefficient (default: " << _gamma << ")\n"
            << " -r | --rho <float>         Coefficient of restitution (default: " << _rho << ")\n"
            << " -N | --number <int>        Number of particles in column (default: " << _N << ")\n"
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

inline nvis::vec2 iterate(const nvis::vec2& x)
{
    return phi(x);
}

inline void do_orbit(color_type* out,
                     size_t ic , size_t niter)
{
	color_type col(drand48(), drand48(), drand48());
	out[ic] = col;
    nvis::vec2 x = i2w(id2i(ic));
    for (size_t n=0 ; n<niter ; ++n) {
        x = iterate(x);
        out[i2id(w2i(x))] = col;
    }
}

void all_orbits(color_type* out) {
    srand48(12345678);
    
    size_t npixels = _width*_height;
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif

    nvis::timer _timer;
    
#pragma openmp parallel for
    for (size_t i=0 ; i<_nseeds ; ++i) {
        size_t ic = floor(drand48()*npixels);
        do_orbit(out, ic, _niter);
        
        int thread_id = 0;
        
#ifdef _OPENMP
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
    _width  = 512;
    _height = 512;
    _T      = 10.;
    _a      = 0.1;
    _rho    = 0.8;
    _gamma  = 0.1;
    _niter  = 500;
    _nseeds = 1000;
    
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
        } else if (arg == "-n" || arg == "--nseeds") {
            if (i == argc-1) {
                usage("missing number of seeds");
            }
            _nseeds = atoi(argv[++i]);
        } else if (arg == "-s" || arg == "--size") {
            if (i >= argc-2) {
                usage("missing image resolution");
            }
            _width = atoi(argv[++i]);
            _height = atoi(argv[++i]);
        } else if (arg == "-N" || arg == "--number") {
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
    
    // set value of dependent parameters
    // g_* = g/N
    _g_star = _G/(double)_N;
    // \gamma = 2\omega^2 a(1+\rho)/g_*
    _omega  = sqrt(_gamma*_g_star/(2.*_a*(1.+_rho)));
    _period = _omega*_T;
    
    
    std::cout << "period omega*T = " << _period << '\n';
    double Mmax = 2.*_a*_omega/(_g_star*_T)*(1.+_rho)/(1.-_rho);
    std::cout << "with the selected parameters, the upper bound on M is " << Mmax << '\n';
    double Mmax2 = _gamma/(_omega*_T*(1.-_rho));
    std::cout << "alternative computation of M yields " << Mmax2 << '\n';
    for (int m=1 ; m<=Mmax ; ++m) {
        double theta_m = acos(m*_omega*_T*(1.-_rho)/_gamma);
        double v_m = m*_omega*_T;
        std::cout << "for m = " << m << ", singularity location is " << nvis::vec2(theta_m, v_m) << '\n';
        std::cout << "\tcheck: phi" << nvis::vec2(theta_m, v_m) << " = " << phi(nvis::vec2(theta_m, v_m)) << '\n'; 
    } 
    
    color_type *image = (color_type*)calloc(_width*_height, sizeof(color_type));
    all_orbits(image);
    
    Nrrd *nout = nrrdNew();
    size_t size[] = {3, _width, _height};
    if (nrrdWrap_nva(nout, (float*)image, nrrdTypeFloat, 3, size)) {
        std::cerr << biffGetDone(NRRD) << '\n';
    }

    double spc[] = {AIR_NAN, _period/(double)(_width-1), (_vmax-_vmin)/(double)(_height-1)};
	nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
    int center[] = {nrrdCenterUnknown, nrrdCenterNode, nrrdCenterNode};
	nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, center);
    if (nrrdSave(_output.c_str(), nout, NULL)) {
        std::cerr << biffGetDone(NRRD) << '\n';
    }
    
    return 0;
}
