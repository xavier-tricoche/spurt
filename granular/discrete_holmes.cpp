// stdlib
#include <algorithm>
#include <exception>
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

// spurt
#include <format/filename.hpp>
#include <image/nrrd_wrapper.hpp>
#include <poincare/metric.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif


typedef nvis::fvec3 Color;
const Color black = Color(0, 0, 0);
const Color red   = Color(1, 0, 0);
const Color green = Color(0, 1, 0);
const Color blue  = Color(0, 0, 1);
const Color white = Color(1, 1, 1);

// Note: constants are capitalized
// constant of gravity
const double _G  = 9.80665;

// human readable Pi constant
const double _PI = M_PI;

// Equation parameters
double      _omega, _T, _gamma, _rho, _g_star, _a;
size_t      _N;

// Computation parameters
size_t _nseeds, _niter;
nvis::bbox2 _bounds;
Color _bg;

// Visualization parameters
double _kappa, _fraction;
nvis::ivec2 _res;
bool _verbose;
int _cmap_id;

spurt::map_metric _metric;

template<typename T1, typename T2, size_t N>
class raster {
public:
    typedef T1                                   coord_type;
    typedef typename nvis::fixed_vector<T1, 2>   pos_type;
    typedef nvis::ivec2                          index_type;
    typedef nvis::bounding_box<pos_type>         box_type;
    typedef T2                                   scalar_type;
    typedef typename nvis::fixed_vector<T2, N>   value_type;

private:
    int to_i(const index_type& id) const {
        return id[0] + id[1]*_res[0];
    }
    index_type to_i(int i) const {
        return index_type(i%_res[0], i/_res[0]);
    }
    int to_i(const pos_type& p) const {
        return to_i(index_type((p-_bounds.min())/_spacing));
    }
    bool check_inside_throw(int i) const {
        if (!inside(i)) {
            std::ostringstream os;
            os << "index " << i << " out of bounds";
            throw std::runtime_error(os.str());
        }
        else return true;
    }
public:
    raster(const index_type& res, const box_type& bounds)
        : _res(res), _bounds(bounds) {
        _data = (scalar_type*)calloc(_res[0]*_res[1]*N, sizeof(scalar_type));
        _spacing = _bounds.size() / pos_type(_res - index_type(1,1));
        _size = _res[0] * _res[1];
    }
    
    size_t size() const {
        return _size;
    }
    
    const box_type& bounds() const {
        return _bounds;
    }
    
    static const size_t channels() {
        return N;
    }
    
    const scalar_type* get_array() const {
        return _data;
    }
    
    const pos_type& spacing() const {
        return _spacing;
    }
    
    const index_type& resolution() const {
        return _res;
    }
    
    const value_type& operator()(int i) const {
        check_inside_throw(i);
        return ((const value_type*)_data)[i];
    }
    
    const value_type& operator()(const index_type& id) const {
        return (*this)(to_i(id));
    }
    
    const value_type& operator()(const pos_type& p) const {
        return (*this)(to_i(p));
    }
    
    value_type& operator()(int i) {
        check_inside_throw(i);
        return ((value_type*)_data)[i];
    }
    
    value_type& operator()(const index_type& id) {
        return (*this)(to_i(id));
    }
    
    value_type& operator()(const pos_type& p) {
        return (*this)(to_i(p));
    }
    
    int shift(const index_type& id) const {
        int i = to_i(id);
        check_inside_throw(i);
        return i*N;
    }
    
    int shift(const pos_type& p) const {
        int i = to_i(p);
        check_inside_throw(i);
        return i*N;
    }
    
    index_type index(int i) const {
        check_inside_throw(i);
        return to_i(i);
    }

    index_type index(const pos_type& p) const {
        int i = to_i(p);
        check_inside_throw(i);
        return to_i(i);
    }
    
    pos_type point(const index_type& id) const {
        int i = to_i(id);
        check_inside_throw(i);
        return _bounds.min() + pos_type(id)*_spacing;
    }
    
    pos_type point(int i) const {
        index_type id = to_i(i);
        check_inside_throw(i);
        return _bounds.min() + pos_type(id)*_spacing;
    }
    
    bool inside(int i) const {
        return (i>=0 && i<_size);
    }
    
    bool inside(const index_type& id) const {
        return inside(to_i(id));
    }
    
    bool inside(const pos_type& p) const {
        return _bounds.inside(p);
    }
    
private:
    size_t       _size;
    index_type   _res;
    scalar_type* _data;
    box_type     _bounds;
    pos_type     _spacing;
};

typedef raster<double, float, 3> color_raster;
typedef raster<double, float, 1> gray_raster;
typedef raster<double, int, 1>   int_raster;

template<typename Func>
struct decorator {
    decorator(const Func& f) : _f(f) {}

    template<typename Raster>
    Raster& operator()(Raster& raster) const {
        typedef typename Raster::index_type index_type;
        typedef typename Raster::pos_type   pos_type;
        typedef typename Raster::value_type value_type;
        for (size_t i=0 ; i<raster.size() ; ++i) {
            index_type id = raster.index(i);
            raster(id) = _f(id, raster.point(i));
        }
        
        return raster;
    }
    
    const Func&  _f;
};

struct color_2d {   
    template<typename Pos>
    Color operator()(const nvis::ivec2&, const Pos& x) const {
        nvis::vec2 loc = nvis::vec2(x) - bounds.min();
        loc /= bounds.size();
        const double& u = loc[0];
        const double& v = loc[1];
        return (1-u)*(1-v)*c0 + u*(1-v)*c1 + u*v*c2 + (1-u)*v*c3;
    }
    
    Color c0, c1, c2, c3;
    nvis::bbox2 bounds;
};

struct color_h {
    template<typename Pos>
    Color operator()(const nvis::ivec2&, const Pos& x) const {
        nvis::vec2 loc = nvis::vec2(x) - bounds.min();
        loc /= bounds.size();
        const double& u = loc[0];
        return (1-u)*c0 + u*c1;
    }
    
    Color c0, c1;
    nvis::bbox2 bounds;   
};

struct color_v {
    template<typename Pos>
    Color operator()(const nvis::ivec2&, const Pos& x) const {
        nvis::vec2 loc = nvis::vec2(x) - bounds.min();
        loc /= bounds.size();
        const double& v = loc[1];
        return (1-v)*c0 + v*c1;
    }
    
    Color c0, c1;
    nvis::bbox2 bounds;   
};

struct color_noise {
    template<typename Pos>
    Color operator()(const nvis::ivec2&, const Pos& x) const {
        return Color(drand48(), drand48(), drand48());
    }
};

struct color_split {
    
    template<typename Pos>
    Color operator()(const nvis::ivec2&, const Pos& x) const {
        double dv = x[1] - vsplit;
        if (dv > 0) {
            double v = (bounds.max()[1] - x[1])/(bounds.max()[1] - vsplit);
            return (1-v)*red + v*white;
        }
        else {
            double v = (x[1] - bounds.min()[1])/(vsplit - bounds.min()[1]);
            return (1-v)*blue + v*white;
        }
    }
    
    double vsplit;
    nvis::bbox2 bounds;
};

inline double W(double t)
{   
    double s = _metric.modulo(t, 0);
    if (s <= _PI) {
        return cos(s);
    } else {
        return 0.;
    }
}

inline nvis::vec2 phi(const nvis::vec2& x)
{   
    const double& theta = x[0];
    const double& v     = x[1];
    nvis::vec2 next(theta+v, _rho*v + _gamma*W(theta+v));
    return next;
}

void iterate(nvis::vec2 x0, size_t niter,
             color_raster::value_type& sum,
             std::vector<color_raster::index_type>& visited,
             const color_raster& input)
{
    typedef color_raster::value_type value_type;
    typedef color_raster::index_type index_type;
    
    nvis::vec2 x = x0;
    index_type id = input.index(_metric.modulo(x));
    sum = input(id);
    visited.push_back(id);
    
    size_t n;
    for (n=0 ; n<niter ; ++n) {
        x = phi(x);
        nvis::vec2 y = _metric.modulo(x);
        if (!input.inside(y)) break;
        id = input.index(y);
        sum += input(id);
        visited.push_back(id);
    }
    sum /= (double)n;
}

std::string me;
void usage(const std::string& what = "")
{
    if (what != "") {
        std::cerr << "Error: " << what << '\n';
    }
    std::cerr
        << "DESCRIPTION: Visualize the discrete dynamical system describing\n"
        << "the motion of the center of mass of a column of particles under\n"
        << "tapping boundary conditions.\n"
        << '\n'
        << "USAGE  : " << me << " [parameters] [options]\n"
        << '\n'
        << "PARAMETERS:\n"
        << " -o | --output <string>    Name of output file\n"
        << '\n'
        << "OPTIONS:\n"
        << " -h | --help               Print this information\n"
        << " -i | --iterations <int>   Number of iterations of the map\n"
        << " -r | --rho <float>        Coefficient of restitution\n"
        << " -g | --gamma <float>      Gamma constant\n"
        << " -a | --amplitude <float>  Normalized tap amplitude\n"
        << " -N | --particles <int>    Number of particles in column\n"
        << " -f | --fraction <float>   Fraction of raster used as seed\n"
        << " -T | --period <float>     Overall tap period\n"
        << " -b | --bounds <float> x2  Velocity bounds\n"
        << " -s | --size <int> x2      Image size\n"
        << " -c | --color <int>        Input color map ID\n"
        << "                           (0: 2d, 1: horiz., 2: vert., 3: noise)\n"
        << " -v | --verbose            Turn on verbose mode\n"
        << std::endl;

    exit(!what.empty());
}

int main(int argc, char* argv[])
{
    srand48(time(0));
    
    _bounds.min() = nvis::vec2(0, -1);
    _bounds.max() = nvis::vec2(0, 1);
    _N          = 20;
    _nseeds     = 200;
    _niter      = 200;
    _T          = 0.5;
    _a          = 0.1;
    _rho        = 0.8;
    _gamma      = 0.1;
    _verbose    = false;
    _res[0]     = 800;
    _res[1]     = 800;
    _cmap_id    = 0;
    _fraction   = -1;
    std::string output_name = "undefined";
    
    me = argv[0];
    for (int i=0 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) usage("Missing output name");
            output_name = spurt::filename::remove_extension(argv[++i]);
        }
        else if (arg == "-a" || arg == "--amplitude") {
            if (i == argc-1) usage("Missing amplitude value");
            _a = atof(argv[++i]);
        }
        else if (arg == "-r" || arg == "--rho") {
            if (i == argc-1) usage("Missing rho value");
            _rho = atof(argv[++i]);
        }
        else if (arg == "-g" || arg == "--gamma") {
            if (i == argc-1) usage("Missing gamma value");
            _gamma = atof(argv[++i]);
        }
        else if (arg == "-N" || arg == "--particles") {
            if (i == argc-1) usage("Missing number of particles");
            _N = atoi(argv[++i]);
        }
        else if (arg == "-f" || arg == "--fraction") {
            if (i == argc-1) usage("Missing number of seeds");
            _fraction = atof(argv[++i]);
        }
        else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) usage("Missing period value");
            _T = atof(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-2) usage("Missing bounds");
            _bounds.min()[1] = atof(argv[++i]);
            _bounds.max()[1] = atof(argv[++i]);
        }
        else if (arg == "-s" || arg == "--size") {
            if (i >= argc-2) usage("Missing image size");
            _res[0] = atoi(argv[++i]);
            _res[1] = atoi(argv[++i]);
        }
        else if (arg == "-i" || arg == "--iterations") {
            if (i == argc-1) usage("Missing number of iterations");
            _niter = atoi(argv[++i]);
        }
        else if (arg == "-c" || arg == "--color") {
            if (i == argc-1) usage("Missing color map ID");
            _cmap_id = atoi(argv[++i]);
        }
        else if (arg == "-v" || arg == "--verbose") {
            _verbose = true;
        }
    }
    
    if (output_name == "undefined") {
        usage("No output file name was provided");
    }
    
    // set value of dependent parameters
    // g_* = g/N
    _g_star = _G/(double)_N;
    // \gamma = 2\omega^2 a(1+\rho)/g_*
    _omega  = sqrt(_gamma*_g_star/(2.*_a*(1.+_rho)));
    
    // set metric for display
    _metric.bounds().min()[0] = 0;
    _metric.bounds().max()[0] = _omega*_T;
    _metric.periodic(0) = true;
    _metric.periodic(1) = false;
    
    _bounds.min()[0] = 0;
    _bounds.max()[0] = _omega*_T;
    
    std::cout << "computation bounds are " << _bounds << '\n';
    
    color_raster input(_res, _bounds);
    if (_cmap_id == 0) {
        color_2d cmap;
        cmap.c0 = black;
        cmap.c1 = blue;
        cmap.c2 = white;
        cmap.c3 = red;
        cmap.bounds = _bounds;
        decorator<color_2d> deco(cmap);
        input = deco(input);
    }  
    else if (_cmap_id == 1) {
        color_h cmap;
        cmap.c0 = green;
        cmap.c1 = red;
        cmap.bounds = _bounds;
        decorator<color_h> deco(cmap);
        input = deco(input);
    }  
    else if (_cmap_id == 2) {
        color_v cmap;
        cmap.c0 = green;
        cmap.c1 = red;
        cmap.bounds = _bounds;
        decorator<color_v> deco(cmap);
        input = deco(input);
    }  
    else if (_cmap_id == 3) {
        color_noise noise;
        decorator<color_noise> deco(noise);
        input = deco(input);
    }
    
    color_raster output(_res, _bounds);
    color_raster phi_raster(_res, _bounds);
    color_raster vec_raster(_res, _bounds);
    int_raster counter(_res, _bounds);
    
    size_t npixels = output.size();

    int nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif
    if (_verbose) std::cout << nthreads << " thread(s) available\n";
    
    typedef std::pair< color_raster::value_type, 
                       std::vector<color_raster::index_type> >
                       orbit_type;
    typedef std::vector<orbit_type> thread_work;
    std::vector<thread_work> work(nthreads);
    
    std::vector<int> seeds(npixels);
    for (int i=0 ; i<npixels ; ++i) seeds[i] = i;
    if (_fraction > 0 && _fraction < 1) {
        std::random_shuffle(seeds.begin(), seeds.end());
        seeds.resize((size_t)floor(npixels*_fraction));
    }
    else if (_fraction == 2) {
        nvis::vec2 x(0,0);
        if (_bounds.inside(x)) {
            seeds.clear();
            nvis::ivec2 id = input.index(x);
            for (int i=0 ; i<_res[0] ; ++i) {
                seeds.push_back(i+id[1]*_res[0]);
            }
        }
    }
    int nseeds = seeds.size();
    
    nvis::timer total_time, time_since_last;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<nseeds ; ++n) {
            int thread = 0;
#if _OPENMP
            thread = omp_get_thread_num();
#endif
            int i = seeds[n];
            nvis::vec2 x0 = output.point(i);
            work[thread].push_back(orbit_type());
            orbit_type& orbit = work[thread].back();            
            iterate(x0, _niter, orbit.first, orbit.second, input);
            
            double dt = total_time.elapsed();
            if (!thread && time_since_last.elapsed()>0.1) {
                std::cout << "\r" << n << "/" << nseeds 
                          << " (" << 100.*n/(double)nseeds << "%) "
                          << "in " << dt << " s. "
                          << "(" << (double)n/dt << " Hz)"
                          << "               \r";
                time_since_last.restart();
            }
        }
    }
    std::cout << '\n';
    
    // reduce
    if (_verbose) std::cout << "writing result to rasters..." << std::flush;
    for (int i=0 ; i<work.size() ; ++i) {
        for (int n=0 ; n<work[i].size() ; ++n) {
            const color_raster::value_type& sum = work[i][n].first;
            const std::vector<color_raster::index_type>& steps = work[i][n].second;
            double _rand = drand48();
            Color _col(_rand, _rand, 1);
            for (int j=0 ; j<steps.size() ; ++j) {
                nvis::vec2 x = _metric.modulo(steps[j]);
                if (!output.inside(x)) continue;
                counter(x)[0]++;
                if (_fraction != 2) 
                    output(x) = _col;
                else 
                    output(x) = sum;
            }
        }
    }
    if (_verbose) std::cout << " done\n";
    
    time_since_last.restart();
    total_time.restart();
    if (_verbose) 
        std::cout << "computing phi transform of phase portrait...\n";
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int n=0 ; n<npixels ; ++n) {
            int thread = 0;
#if _OPENMP
            thread = omp_get_thread_num();
#endif
            color_raster::index_type id = vec_raster.index(n);
            nvis::vec2 x = input.point(id);
            nvis::vec2 y = phi(x);
            nvis::vec2 v = y-x;
            nvis::vec2 z = _metric.modulo(y);
            if (phi_raster.inside(z)) phi_raster(z) = input(x); 
            vec_raster(id)[0] = v[0];
            vec_raster(id)[1] = 0;
            vec_raster(id)[2] = v[1];
            
            double dt = total_time.elapsed();
            if (!thread && time_since_last.elapsed()>0.1) {
                std::cout << "\r" << n << "/" << npixels 
                          << " (" << 100.*n/(double)npixels << "%) "
                          << "in " << dt << " s. "
                          << "(" << (double)n/dt << " Hz)"
                          << "               \r";
                time_since_last.restart();
            }
        }
    }
    std::cout << '\n';
    
    for (int n=0 ; n<npixels ; ++n) {
        const nvis::vec2 x = vec_raster.point(n);
        color_raster::index_type id = vec_raster.index(n);
        Color v = vec_raster(id);
        if (nvis::norm(v) == 0) {
            std::cout << "null vector at " << id << " = " << x << '\n';
            std::cout << "phi(" << x << ") = " << phi(x) << '\n';
            std::cout << "||x-y|| = " << nvis::norm(phi(x)-x) << '\n';
        }
    }
    
    spurt::nrrd_params<float, 3> params;
    params.sizes()[0] = 3;
    params.sizes()[1] = output.resolution()[0];
    params.sizes()[2] = output.resolution()[1];
    params.mins()[0] = 0;
    params.mins()[1] = output.bounds().min()[0];
    params.mins()[2] = output.bounds().min()[1];
    params.spacings()[0] = 1;
    params.spacings()[1] = output.spacing()[0];
    params.spacings()[2] = output.spacing()[1];
    
    spurt::writeNrrdFromParams(const_cast<float*>(output.get_array()), output_name + "-col.nrrd", params);
    std::cout << output_name + "-col.nrrd was exported\n";

    spurt::writeNrrdFromParams(const_cast<float*>(phi_raster.get_array()), output_name + "-phi.nrrd", params);
    std::cout << output_name + "-phi.nrrd was exported\n";

    spurt::writeNrrdFromParams(const_cast<float*>(input.get_array()), output_name + "-ref.nrrd", params);
    std::cout << output_name + "-ref.nrrd was exported\n";
    
    spurt::writeNrrdFromParams(const_cast<float*>(vec_raster.get_array()), output_name + "-vec.nrrd", params);
    std::cout << output_name + "-vec.nrrd was exported\n";
    
    spurt::nrrd_params<int, 2> paramsi;
    paramsi.sizes()[0] = output.resolution()[0];
    paramsi.sizes()[1] = output.resolution()[1];
    paramsi.mins()[0] = output.bounds().min()[0];
    paramsi.mins()[1] = output.bounds().min()[1];
    paramsi.spacings()[0] = output.spacing()[0];
    paramsi.spacings()[1] = output.spacing()[1];
    
    spurt::writeNrrdFromParams(const_cast<int*>(counter.get_array()), output_name + "-cnt.nrrd", paramsi);
    std::cout << output_name + "-cnt.nrrd was exported\n";
    
    return 0;
}
