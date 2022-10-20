#include <iostream>
#include <utility>
#include <algorithm>

#include <image/nrrd_wrapper.hpp>
#include <image/probe.hpp>
#include <math/fixed_vector.hpp>
// #include <vis/streamline.hpp>
#include <misc/progress.hpp>
#include <data/raster.hpp>
#include <format/filename.hpp>
#include <misc/option_parse.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

using namespace boost::numeric::odeint;

typedef spurt::raster_grid<2, double, size_t> grid_t;
typedef spurt::raster_data<vec3, 2, double, size_t> raster_t;
typedef vec3 state_t; // (x, y, inside)


std::string in_name, out_name;
double t0=0, t;
double eps=1.0e-6;
ivec2 res(256, 256);
bool verbose=false;
vec4 bounds_as_array(0);
bbox2 bounds;

struct nrrd_vector_field {
    typedef spurt::gage_interface::vector_wrapper wrapper_t;
    
    nrrd_vector_field(Nrrd* nin) 
        : m_wrapper(nin, spurt::gage_interface::BC_INTERP, false) {
            m_wrapper.use_world();
        }
            
    bool operator()(const vec3& x, vec3& f) const {
        return m_wrapper.value(x, f);
    }
    
    wrapper_t m_wrapper;
};

struct odeint_rhs {    
    odeint_rhs( const nrrd_vector_field& field ) 
        : m_field(field) {}
    
    odeint_rhs( const odeint_rhs& other ) 
        : m_field( other.m_field ) {}
    
    void operator()( const state_t& x, state_t& dxdt, double t) const {
        vec3 v;
        bool ok = m_field( vec3(x[0], x[1], t), v );
        if (!ok) {
            dxdt = state_t(0, 0, 1);
        }
        else {
            dxdt[0] = v[0];
            dxdt[1] = v[1];
            dxdt[2] = 0;
        }
    }
    
    const nrrd_vector_field& m_field;
};

struct left_domain_condition {
    bool operator()(const state_t& x) {
        return x[2] != 0;
    }
};


template<class RHS, class Condition>
std::pair<double, state_t>
find_condition(state_t& x0, RHS rhs, Condition cond, const double t_start,
               const double t_end, const double dt, const double prec=1.0e-6) {
    auto stepper = make_dense_output(eps, eps, runge_kutta_dopri5<state_t>());
    auto ode_range = make_adaptive_range(std::ref(stepper), rhs, x0, t_start,
                                         t_end, dt);
    auto found_iter = std::find_if(ode_range.first, ode_range.second, cond);
    if (found_iter == ode_range.second) {
        return std::make_pair(t_end+dt, x0);
    }
    
    double t0 = stepper.previous_time();
    double t1 = stepper.current_time();
    double t_m;
    state_t x_m;
    while (t1 - t0 > prec) {
        t_m = 0.5 * (t0 + t1);  // get the mid point time
        stepper.calc_state(t_m, x_m); // obtain the corresponding state
        if (cond(x_m))
            t1 = t_m;  // condition changer lies before midpoint
        else
            t0 = t_m;  // condition changer lies after midpoint
    }
    // we found the interval of size prec, take it's midpoint as final guess
    t_m = 0.5 * (t0 + t1);
    stepper.calc_state(t_m, x_m);
    return std::make_pair(t_m, x_m);
}

void initialize(int argc, char* argv[]) {
    namespace xcl = spurt::command_line;
    xcl::option_traits 
        required_group(true, false, "Required Options"), 
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");
        
    xcl::option_parser parser(argv[0],
        "Integrate flow map in time-dependent vector field stored as 3D NRRD");
        
    try {
        // parser.use_default_symbol();
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", in_name, "Input filename", required_group);
        parser.add_value("output", out_name, "Output filename", optional_group);
        parser.add_value("tinit", t0, t0, "Initial time", optional_group);
        parser.add_value("t", t, "Integration length", required_group);
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_tuple<4>("bounds", bounds_as_array, "Selected region", optional_group);
        parser.add_flag("verbose", verbose, "Verbose output", optional_group);
        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options enteredso far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

int main(int argc, char* argv[]) {
    initialize(argc, argv);
    
    double dt = (t > 0 ? eps : -eps); 
    
    int nthreads = 1;
#if _OPENMP
    nthreads = omp_get_max_threads();
#endif

    // bool fwd = (t > 0);
    
    spurt::ProgressDisplay progress;
    
    Nrrd* nin = spurt::nrrd_utils::readNrrd(in_name);
    
    std::vector<nrrd_vector_field*> nrrd_vfs(nthreads);
    for (int i=0; i<nthreads; ++i) {
        nrrd_vfs[i] = new nrrd_vector_field(nin);
    }
    
    if (bounds_as_array[2]>bounds_as_array[0] &&
        bounds_as_array[3]>bounds_as_array[1]) {
        bounds.min()=vec2(bounds_as_array[0], bounds_as_array[1]);
        bounds.max()=vec2(bounds_as_array[2], bounds_as_array[3]);    
    }
    else {
        bounds.min()[0] = nin->axis[1].min;
        bounds.min()[1] = nin->axis[2].min;
        bounds.max()[0] = nin->axis[1].min + (nin->axis[1].size-1)*nin->axis[1].spacing;
        bounds.max()[1] = nin->axis[2].min + (nin->axis[2].size-1)*nin->axis[2].spacing;
    }
    
    double minz = nin->axis[3].min;
    double maxz = minz + (nin->axis[3].size-1)*nin->axis[3].spacing;
    if (verbose) {
        std::cout << "bounds: " << bounds << '\n';
        std::cout << "z range: " << minz << " -> " << maxz << '\n';
    }
    
    grid_t sampling_grid(res, bounds);
    raster_t fmap(sampling_grid);
    std::vector<bool> stopped(sampling_grid.size(), false);
    
    size_t nfailed = 0;
    
    left_domain_condition left_c;
    
    progress.start(sampling_grid.size(), "flow map");
#pragma omp parallel
    {
    #pragma omp for schedule(dynamic,1)
    for (size_t i=0; i<sampling_grid.size(); ++i) {        
        #if _OPENMP
        const int thread=omp_get_thread_num();
        #else
        const int thread=0;
        #endif
        
        if (!thread) progress.update(i);

        odeint_rhs rhs(*nrrd_vfs[thread]);
        
        nvis::vec2 x0 = sampling_grid[i];
        state_t x(x0[0], x0[1], 0);
        state_t x_cond;
        double t_cond;
        std::tie(t_cond, x_cond) = find_condition(x, rhs, left_c, t0, t0+t, dt, 1.0e-6);
          
        if ((dt > 0 && t_cond > t0+t) || (dt < 0 && t_cond < t0+t)) {
            fmap[i] = vec3(x[0], x[1], t0+t);
        }
        else {
            ++nfailed;
            fmap[i] = vec3(x_cond[0], x_cond[1], t_cond);
        }

#if 0        
        // std::vector<vec3> steps;
        // observer obs(steps, x, t0);
        
        try {
            integrate_adaptive(stepper, rhs, x, t0, t0+t, eps);
        }
        catch(std::runtime_error& e) {
            if (verbose) {
                std::cerr << "integrate_adaptive threw exception: "
                          << e.what() << '\n';
            }
            ++nfailed;
        }
        fmap[i] = x;
#endif

#if 0        
        vec3 x0(p[0], p[1], t0);
        streamline sl(x0, t0);
        streamline::no_stop go_on;
        sl.record = false;
        sl.reltol = sl.abstol = eps;
        sl.stepsz = 0;
        try {
            streamline::state state = sl.advance(*nrrd_vfs[thread], t0+t, go_on);
            fmap[i] = sl(fwd ? sl.t_max() : sl.t_min());
            if (verbose && state != streamline::OK) {
                std::cout 
                    << "integration from " << x0 << " failed: "
                    << "return code: " << state << std::endl;
                ++nfailed;
            }
        }
        catch( std::runtime_error& e) {
            stopped[i] = true;
        }
#endif
    }
    }
    progress.end();
    
    std::cout << nfailed << " incomplete trajectories\n";
    
    for (int i=0; i<nthreads; ++i) {
        delete nrrd_vfs[i];
    }
    
    if (out_name.empty()) {
        out_name = spurt::filename::remove_extension(in_name) + "-fmap.nrrd";
    }
    fmap.save_as_nrrd(out_name);
    
    return 0;
}

