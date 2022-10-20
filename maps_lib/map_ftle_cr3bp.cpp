// #define NO_TBB
#define VERBOSE_OUTPUT 0

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

#include <util/timer.hpp>
#include <misc/meta_utils.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <data/raster.hpp>
#include <format/filename.hpp>
#include <vtk/vtk_io_helper.hpp>
#include <vtk/vtk_data_helper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/dopri5.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef NO_TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <tbb/mutex.h>
#endif

#ifndef NO_TBB
tbb::mutex output_mutex;
tbb::mutex progress_mutex;
tbb::mutex orbit_mutex;
tbb::atomic<size_t> progress_counter;
#else
size_t progress_counter;
#endif

void write_to_ostream(std::ostream& os, const std::string& str) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(output_mutex);
    #endif
        os << str << '\n';
    }
}

void update_progress(spurt::ProgressDisplay& progress) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(progress_mutex);
    #endif
        progress.update(progress_counter);
    }
}

template< typename T >
void append_orbit(std::vector< std::vector<T> >& out, const std::vector<T>& in) {
    {
    #ifndef  NO_TBB
        tbb::mutex::scoped_lock lock(orbit_mutex);
    #endif
        out.push_back(in);
    }
}

//[ vexcl_state_types

template< typename T >
struct value_traits {
    typedef T                                                 value_type;
    typedef std::array<value_type, 2>                         array_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   column_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   vector_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   state_type;
    typedef Eigen::Matrix<value_type, 2, 2>                   matrix_type;

    // sampling grid for Poincare section
    typedef spurt::raster_grid<2, value_type>                grid_type;
    typedef typename grid_type::point_type                    position_type;
    typedef nvis::bounding_box< position_type >               bounds_type;
    typedef typename grid_type::coord_type                    coordinates_type;
    typedef typename grid_type::point_type                    point_type;
    typedef spurt::raster_data<column_type, 2, value_type>   vector_raster_type;
};
//]

const double Jupiter_Europa_mu = 2.528017705e-5;
const double Earth_Moon_mu = 0.012150571430596;
const double invalid_double = std::numeric_limits<double>::max();
constexpr double crossing_eps = 1.0e-6;

//[ vexcl_system
std::string filename;
std::array<size_t, 2> res = {{64, 64}}; // flow map sampling resolution
bool        verbose = false;
double      dt = 1.0e-2;        // initial integration time step
double      t_max = 10.0;       // propagation duration
double      eps = 1.0e-6;       // precision
double      C = 2.96;           // Jacobi constant (low'ish energy in EM)
double      mu = Earth_Moon_mu; // gravity constant
double      xmin=0.4, xmax=1.1, xdmin=-2.5, xdmax=2.5; // sampling bounds
bool        double_prec = true;
size_t      nb_iterations = 10;
size_t      iteration_stride = 1;
bool        compute_ftle = true;
bool        do_streamlines = false;
bool        do_monitor = true;
bool        do_discrete_orbits = true;

typedef nvis::vec3 world_pos_type;

template<typename T, size_t N>
inline std::array<T, N> array(const nvis::fixed_vector<T, N>& v) {
    std::array<T, N> a;
    for (int i=0; i<N; ++i) a[i]=v[i];
    return a;
}

template<typename T, size_t N>
inline std::array<T, N> array(const std::array<T, N>& a) {
    return a;
}

template<typename T1, typename T2, size_t N>
inline std::array<T2, N> convert(const std::array<T1, N>& a) {
    std::array<T2, N> r;
    std::copy(a.begin(), a.end(), r.begin());
    return r;
}

template<typename value_t>
inline value_t yd(value_t x, value_t xd) {
    value_t ysq = 2*(1-mu)/std::abs(x+mu) + 2*mu/std::abs(x-1+mu) + x*x - xd*xd - C;
    if (ysq<0) throw std::exception();
    return sqrt(ysq);
}

template<typename T, typename Vector = Eigen::Matrix<T, 2, 1> >
struct raster_wrapper {
    typedef T                                      scalar_type;
    typedef spurt::raster_grid<2, T>              grid_type;
    typedef spurt::raster_data<Vector, 2, T>      raster_type;
    typedef Vector                                 value_type;
    typedef typename grid_type::coord_type         coord_type;

    raster_wrapper(const raster_type& r)
        : _raster(r), _res(r.grid().resolution()) {}

    const value_type& operator()(long int i, long int j) const {
        coord_type c(i,j);
        c+=_res;
        return _raster.operator()(c[0]%_res[0], c[1]%_res[1]);
    }

    const raster_type& _raster;
    const coord_type _res;
};


template<typename T, size_t N, typename Vector=Eigen::Matrix<T, 2, 1>, typename Matrix=Eigen::Matrix<T, 2, 2> >
void gradient(const spurt::raster_data<Vector, 2, T>& fmap,
              const std::array<T, N>& kernel,
              std::vector<Matrix>& grad,
              std::vector<T>& ftle)
{
    typedef Matrix mat_t;

    raster_wrapper<T, Vector> raster(fmap);
    grad.resize(fmap.size());
    ftle.resize(fmap.size());

    auto res = fmap.grid().resolution();
    auto spc = fmap.grid().spacing();
    int shift = -static_cast<long int>(N/2);

    size_t npoints = fmap.grid().size();

#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif

    spurt::ProgressDisplay progress(true);

    progress_counter = 0;
    progress.fraction_on();
    progress.start(npoints, "FTLE computation");
#ifndef NO_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, npoints),
                      [&](tbb::blocked_range<size_t> r) {
    for (size_t l=r.begin(); l!=r.end(); ++l) {
#else
    for (size_t l=0 ; l<npoints; ++l) {
#endif
        ++progress_counter;
        update_progress(progress);
        mat_t m = mat_t::Zero();
        nvis::ivec2 coords = fmap.grid().coordinates(l);
        const size_t& i = coords[0];
        const size_t& j = coords[1];
        bool failed = false;
        for (int n=0; n<N; ++n) {
            long int ii = i+shift+n;
            long int jj = j+shift+n;
            if (ii <= 0 || ii >= res[0] || jj <= 0 || jj >= res[1]) {
                ftle[l] = 0;
                grad[l] = m;
                failed = true;
                break;
            }
            m.col(0)[0] += kernel[n]*raster(ii, j)[0];
            m.col(0)[1] += kernel[n]*raster(ii, j)[1];
            m.col(1)[0] += kernel[n]*raster(i, jj)[0];
            m.col(1)[1] += kernel[n]*raster(i, jj)[1];
        }
        if (failed) continue;
        m.col(0) /= spc[0];
        m.col(1) /= spc[1];
        grad[l] = m;
        if (false && verbose) {
           std::cout << "grad[" << fmap.grid().index(i,j) << "]=" << m << '\n';
        }

        Eigen::JacobiSVD<mat_t,Eigen::NoQRPreconditioner> svd(m);
        auto val = 2*log(svd.singularValues()[0]);
        if (std::isnan(val) || std::isinf(val)) val = 0;
        ftle[fmap.grid().index(i,j)]=val;
    }
#ifndef NO_TBB
    });
#endif
    progress.end();
    std::cout << progress << '\n';
}

//]

template<typename T, typename Value_>
void write(const std::string& filename, const std::vector<Value_>& data,
           size_t valsize, const typename value_traits<T>::grid_type& grid)
{
    typedef value_traits<T> types;
    typedef typename types::value_type value_t;
    std::ofstream outf(filename.c_str(), std::ofstream::binary);
    std::ostringstream os;
    os << "NRRD0001\n"
       << "# Complete NRRD file format specification at:\n"
       << "# http://teem.sourceforge.net/nrrd/format.html\n"
       << "type: " << spurt::type2string<value_t>::type_name() << "\n"
       << "dimension: " << (valsize > 1 ? 3 : 2) << '\n';
    os << "sizes:";
    if (valsize > 1) os << " " << valsize;
    for (int i=0; i<grid.dim; ++i) os << " " << grid.resolution()[i];
    // os << "\nspacings:";
    // if (valsize > 1) os << " nan";
    // for (int i=0; i<grid.dim; ++i) os << " " << grid.spacing()[i];
    // os << "\naxis mins:";
    // if (valsize > 1) os << " nan";
    // for (int i=0; i<grid.dim; ++i) os << " " << grid.bounds().min()[i];
    os << "\nspace dimension: 3"
       << "\nspace origin: (0,0,0)"
       << "\nspace directions: (1,0,0) (0,0,1)";
    os  << "\nendian: little\n"
        << "encoding: raw\n\n";
    outf.write(os.str().c_str(), os.str().size());

    if (verbose) {
        std::cout << "NRRD header=\n" << os.str() << '\n';
    }
    size_t size=data.size()*sizeof(Value_);
    if (verbose) {
        std::cout << "sizeof(Value_)=" << sizeof(Value_) << '\n'
                  << "size=" << size << '\n';
    }
    outf.write((char*)&data[0], size);
    outf.close();
}

std::string human_readable_size(size_t bytes) {
    const std::string words[] = {"b", "kb", "mb", "gb", "tb", "pb", "eb"};
    std::vector<size_t> sub;
    while (bytes > 0) {
        size_t q = bytes / 1024;
        size_t r = bytes % 1024;
        sub.push_back(r);
        bytes = q;
    }
    std::ostringstream os;
    while (!sub.empty()) {
        size_t n = sub.back();
        if (n != 0) os << n << " " << words[sub.size()-1] << " ";
        sub.pop_back();
    }
    return os.str();
}

template<typename T=double>
struct dopri45_wrapper {
    typedef T value_type;
    typedef nvis::fixed_vector<T, 4> state_type; // vector type required by nvis::dopri5
    typedef nvis::fixed_vector<T, 2> position_type;
    typedef nvis::dopri5<state_type> ode_solver_type;
    typedef typename ode_solver_type::step step_type;

    struct cr3bp {
        value_type m_C, m_mu;
        cr3bp(value_type C, value_type mu) : m_C(C), m_mu(mu) {}
        cr3bp(const cr3bp& other) : m_C(other.m_C), m_mu(other.m_mu) {}

        state_type operator()(value_type t, const state_type& state) const {
            const value_type& x  = state[0];
            const value_type& y  = state[1];
            const value_type& dx = state[2];
            const value_type& dy = state[3];

            value_type d = sqrt((x+m_mu)*(x+m_mu) + y*y);
            value_type r = sqrt((x+m_mu-1)*(x+m_mu-1) + y*y);
            value_type d3 = d*d*d;
            value_type r3 = r*r*r;

            return state_type(dx, dy,
                2*dy + x + (m_mu-1)*(x+m_mu)/d3 - m_mu*(x+m_mu-1)/r3,
                -2*dx + y + (m_mu-1)*y/d3 - m_mu*y/r3
            );
        }
    };

    static
    void integrate(state_type& last_state, value_type& last_time,
                   position_type& crossing_pos, value_type& crossing_time,
                   const cr3bp& rhs, state_type& y, value_type t0,
                   value_type hinit, value_type hmax, value_type eps,
                   std::vector<world_pos_type>& steps,
                   bool forward=true) {
        ode_solver_type solver;
        cr3bp my_rhs(rhs);
        if (eps > 0) {
            solver.reltol = solver.abstol = eps;
        }
        if (hmax > 0) {
            solver.h_max = hmax;
        }
        value_type duration = (forward ? 1.0e9 : -1.0e-9); // some reasonably large value to allow for a sensible integration length to complete a full iteration of the map
        solver.t = t0;
        solver.t_max = t0 + duration;
        solver.y = y;
        step_type step;
        state_type save = y;
        int iter=0;
        steps.clear();
        try {
            while (true) {
                ++iter;
                auto result = solver.do_step(my_rhs, step);
                if (result == ode_solver_type::OK) {
                    value_type y0 = step.y0()[1];
                    value_type y1 = step.y1()[1];
                    value_type t0 = step.t0();
                    value_type t1 = step.t1();
                    value_type yy;
                    value_type tt;
                    if (do_streamlines || do_monitor) {
                        state_type tmp = step.y1();
                        steps.push_back(world_pos_type(tmp[0], tmp[1], tmp[2]));
                    }
                    if ( (forward && y1<0) || (!forward && y1>0) ) continue;
                    if (y0 * y1 < 0) {
#if VERBOSE_OUTPUT > 1
                        std::ostringstream os;
                        os << "\nzero crossing detected: y0=" << y0 << ", y1=" << y1 << ", t0=" << t0 << ", t1=" << t1 << '\n';
#endif
                        // cross y=0 section in prescribed direction
                        while (std::abs(t1-t0) > crossing_eps) {
                            // bisection method
                            // (1-u)*y0 + u*y1 = 0
                            // y0 + u*(y1-y0) = 0
                            // u = -y0/(y1 - y0)
                            // value_type u = -y0/(y1-y0);

                            // binary search method
                            value_type u = 0.5;
                            tt = (1-u)*t0 + u*t1;
                            yy = step.y(tt, 1);
                            if (yy*y0 > 0) {
                                // y and y0 have same sign.
                                // zero crossing between t and t1
                                y0 = yy;
                                t0 = tt;
                            }
                            else if (yy*y1 > 0) {
                                // y and y1 have same sign
                                // zero crossing between t0 and t
                                y1 = yy;
                                t1 = tt;
                            }
                            else break;
#if VERBOSE_OUTPUT > 1
                            os << "u=" << u << ", y0=" << y0 << ", y1=" << y1 << ", t0=" << t0 << ", t1=" << t1 << '\n';
                            write_to_ostream(std::cerr, os.str());
                            os.clear();
                            os.str("");
#endif
                            if (std::abs(yy) < crossing_eps) break;
                        }
#if VERBOSE_OUTPUT > 1
                        write_to_ostream(std::cerr, os.str());
#endif

#if VERBOSE_OUTPUT > 0
                        {
                            std::ostringstream os;
                            os << "integrate: went from " << save << " to " << step.y1() << " via " << step.y(t0) << " in " << iter << " iterations\n";
                            write_to_ostream(std::cerr, os.str());
                        }
#endif
                        last_state = step.y1();
                        last_time = step.t1();
                        state_type tmp = step.y(tt);
                        crossing_pos = position_type(tmp[0], tmp[2]);
                        crossing_time = tt;
                        return;
                    }
                }
                else if (result == ode_solver_type::T_MAX_REACHED ||
                         result == ode_solver_type::STIFFNESS_DETECTED ||
                         result == ode_solver_type::STEPSIZE_UNDERFLOW) {
                    throw std::runtime_error("Unable to complete an iteration of the map");
                }
            }
        }
        catch(std::runtime_error& e) {
            std::ostringstream os;
#if VERBOSE_OUTPUT > 0
            os << "dopri45_wrapper::integrate: exception caught: " << e.what();
            write_to_ostream(std::cerr, os.str());
#endif
            throw nvis::invalid_position_exception(os.str());
        }
    }

    static bool did_intersect(const step_type& step, bool forward=true) {
        const value_type& y0 = step.y0()[1];
        const value_type& y1 = step.y1()[1];
        if (y0*y1 > 0) return false;
        if (forward)
            return (y0<0);
        else
            return (y0>0);
    }

    static void find_intersection(position_type& xdx, value_type& t,
                                  const step_type& step) {
        // cross y=0 section in prescribed direction
        value_type y0 = step.y0()[1];
        value_type y1 = step.y1()[1];
        value_type t0 = step.t0();
        value_type t1 = step.t1();
        value_type _y;
        value_type _t;
        while (std::abs(t1-t0) > crossing_eps) {
            // binary search method
            _t = (t0 + t1)/2;
            _y = step.y(_t, 1);
            if (_y*y0 > 0) {
                // _y and y0 have same sign.
                // zero crossing between _t and t1
                y0 = _y;
                t0 = _t;
            }
            else if (_y*y1 > 0) {
                // _y and y1 have same sign
                // zero crossing between t0 and _t
                y1 = _y;
                t1 = _t;
            }
            else break;

    #if VERBOSE_OUTPUT > 1
            os << ", y0=" << y0 << ", y1=" << y1 << ", t0=" << t0 << ", t1=" << t1 << '\n';
            write_to_ostream(std::cerr, os.str());
            os.clear();
            os.str("");
    #endif
            if (std::abs(_y) < crossing_eps) break;
        }
    #if VERBOSE_OUTPUT > 1
        write_to_ostream(std::cerr, os.str());
    #endif

        xdx[0] = step.y(_t, 0);
        xdx[1] = step.y(_t, 2);
        t = _t;
    }

    static
    void iterate(state_type& last_state, value_type& last_time,
                 std::vector<position_type>& returns,
                 std::vector<value_type>& times,
                 std::vector<world_pos_type>& steps,
                 const cr3bp& rhs, const state_type& y, value_type t0,
                 size_t n_returns,
                 value_type hinit, value_type hmax, value_type eps,
                 bool forward=true) {
        ode_solver_type solver;
        cr3bp my_rhs(rhs);
        if (eps > 0) {
            solver.reltol = solver.abstol = eps;
        }
        if (hmax > 0) {
            solver.h_max = hmax;
        }
        value_type duration = (forward ? 1.0e9 : -1.0e-9); // some reasonably large value to allow for a sensible integration length to complete a full iteration of the map
        solver.t = t0;
        solver.t_max = t0 + duration;
        solver.y = y;
        step_type step;
        state_type save = y;
        int iter=0;
        steps.clear();
        returns.clear();
        times.clear();
        try {
            while (returns.size() < n_returns) {
                ++iter;
                auto result = solver.do_step(my_rhs, step);
                if (result == ode_solver_type::OK) {
                    if (do_streamlines || do_monitor) {
                        state_type tmp = step.y1();
                        steps.push_back(world_pos_type(tmp[0], tmp[1], tmp[2]));
                    }
                    if (!did_intersect(step, forward)) continue;
#if VERBOSE_OUTPUT > 1
                    std::ostringstream os;
                    os << "\nzero crossing detected: y0=" << y0 << ", y1=" << y1 << ", t0=" << t0 << ", t1=" << t1 << ", iter=" << iter << '\n';
#endif
                    position_type pp;
                    value_type tt;
                    find_intersection(pp, tt, step);
                    returns.push_back(pp);
                    times.push_back(tt);
                }
                else if (result == ode_solver_type::T_MAX_REACHED ||
                         result == ode_solver_type::STIFFNESS_DETECTED ||
                         result == ode_solver_type::STEPSIZE_UNDERFLOW) {
                    throw std::runtime_error("Unable to complete this iteration of the map");
                }
            }
        }
        catch(std::runtime_error& e) {
            std::ostringstream os;
#if VERBOSE_OUTPUT > 0
            os << "dopri45_wrapper::integrate: exception caught: " << e.what();
            write_to_ostream(std::cerr, os.str());
#endif
            throw nvis::invalid_position_exception(os.str());
        }
        last_state = step.y1();
        last_time = step.t1();
    }
};

void export_streamlines(const std::vector< std::vector< world_pos_type > >& lines,
                        const std::string& filename) {
#if VERBOSE_OUTPUT > 0
    std::cout << lines.size() << " streamlines on output\n";
#endif

#if VERBOSE_OUTPUT > 1
    for (int i=0; i<lines.size(); ++i) {
        std::cout << "streamline #" << i << " contains " << lines[i].size() << " vertices\n";
    }
#endif
    VTK_SMART(vtkPolyData) pdata = vtk_utils::make_polylines(lines, 0);
    vtk_utils::saveVTK(pdata, filename);
}

template<typename T = double>
void run(const size_t niterations) {
    typedef T                                                 value_type;
    typedef dopri45_wrapper<value_type>                       wrapper_type;
    typedef typename wrapper_type::state_type                 state_type;
    typedef typename wrapper_type::ode_solver_type            solver_type;
    typedef typename wrapper_type::step_type                  step_type;
    typedef typename wrapper_type::cr3bp                      rhs_type;
    typedef nvis::fixed_vector<value_type, 2>                 position_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   column_type;
    typedef Eigen::Matrix<value_type, 2, 2>                   matrix_type;
    typedef spurt::raster_grid<2, value_type>                grid_type;
    typedef spurt::raster_data<position_type, 2, value_type> vector_raster_type;
    typedef nvis::bounding_box<position_type>                 bounds_type;
    typedef typename grid_type::coord_type                    coordinates_type;

    // set up number of system, time step and integration time
    const size_t n_total = res[0]*res[1];

    bounds_type bounds(position_type(xmin, xdmin), position_type(xmax, xdmax));
    grid_type grid(coordinates_type(res[0], res[1]), bounds);
    const auto spacing = grid.spacing();

    std::cout << "grid spacing = " << spacing << '\n';
    std::cout << "grid bounds = " << bounds << '\n';

    // initialize the state of the CR3BP system
    std::map<size_t, size_t> glob2loc;
    std::vector<state_type>  valid_states;
    std::vector< std::vector<world_pos_type> > continuous_orbits;
    std::vector< std::vector<position_type> > discrete_orbits;
    size_t id=0;
    for (size_t i=0; i<n_total; ++i) {
        try {
            position_type p = grid(grid.coordinates(i));
            value_type dydt = yd(p[0], p[1]);
            glob2loc[i] = id++;
            valid_states.push_back(state_type(p[0], 0, p[1], dydt));
            discrete_orbits.push_back(std::vector<position_type>());
            discrete_orbits.back().push_back(p);
            continuous_orbits.push_back(std::vector<world_pos_type>());
            continuous_orbits.back().push_back(world_pos_type(p[0], 0, p[1]));
        }
        catch(...) { continue; }
    }

    const size_t n_valid = valid_states.size();
    if (verbose) {
        std::cout << "There are " << n_valid << " valid initial conditions\n";
    }

    rhs_type rhs(C, mu);
    std::vector<bool> stopped(n_valid, false);

    spurt::ProgressDisplay progress(true);

    size_t last_iteration = 0;

    while (last_iteration < nb_iterations) {
        size_t to_do = std::min(iteration_stride, nb_iterations-last_iteration);
        progress.fraction_on();
        progress_counter = 0;
        if (verbose) std::cout << "Computing iterations " << last_iteration+1 << " to " << last_iteration + to_do << "...\n";
        progress.begin(n_valid, "Integrated orbits");
#ifndef NO_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n_valid),
                          [&](tbb::blocked_range<size_t> r) {
            for (size_t n=r.begin(); n!=r.end(); ++n) {
#else
            for (size_t n=0; n<n_valid; ++n) {
#endif
                ++progress_counter;
                update_progress(progress);
                if (stopped[n]) continue;
                state_type y0 = valid_states[n];
                std::vector<world_pos_type> steps;
                try {
                    state_type last;
                    std::vector<position_type> returns;
                    std::vector<value_type> times;
                    value_type last_time;
                    wrapper_type::iterate(last, last_time, returns, times, steps,
                                          rhs, y0, 0, to_do, 0, 0, eps, true);
                    valid_states[n] = last;
                    if (do_streamlines) {
                        std::copy(steps.begin(), steps.end(),
                                  std::back_inserter(continuous_orbits[n]));
                        // append_orbit(continuous_orbits, steps);
                    }
                    if (do_discrete_orbits) {
                        std::copy(returns.begin(), returns.end(),
                                  std::back_inserter(discrete_orbits[n]));
                    }
                }
                catch(...) {
                    stopped[n] = true;
                    if (do_monitor) {
                        std::copy(steps.begin(), steps.end(),
                                  std::back_inserter(continuous_orbits[n]));
                    }
                }
            }
#ifndef NO_TBB
        });
#endif
        progress.end();
        std::cout << progress << '\n';
        progress.start(n_total, "Filling array");
        progress_counter = 0;
        std::vector<position_type> fmap(n_total, position_type(0,0));
#ifndef NO_TBB
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n_total),
                          [&](tbb::blocked_range<size_t> r) {
            for (size_t n=r.begin(); n!=r.end(); ++n) {
#else
            for (size_t n=0; n<n_total; ++n) {
#endif
                ++progress_counter;
                update_progress(progress);
                auto what = glob2loc.find(n);
                if (what != glob2loc.end()) {
                    size_t i = what->second;
                    // if (stopped[i]) continue;
                    fmap[n] = discrete_orbits[i].back();
                }
            }
#ifndef NO_TBB
        });
#endif
        progress.end();
        std::cout << progress << '\n';

        std::string iter_str;
        {
            std::ostringstream oss;
            oss << std::setw(5) << std::setfill('0') << (last_iteration + to_do);
            iter_str = oss.str();
        }

        std::string name = filename + "-i=" + iter_str + "-" + std::to_string(res[0]) + 'x' + std::to_string(res[1]);
        write<value_type, position_type>(name + "-fmap.nrrd", fmap, 2, grid );
        if (verbose) {
            std::cout << name + "-fmap.nrrd has been exported\n";
        }

        vector_raster_type raster(grid, fmap);
        std::vector<matrix_type> grad;
        std::vector<value_type> ftle;
        std::array<value_type, 3> kernel = {-0.5, 0, 0.5};
        gradient<value_type, 3, position_type, matrix_type>(raster, kernel, grad, ftle);
        write<value_type, matrix_type>( name + "-grad.nrrd", grad, 4, grid );
        if (verbose) {
            std::cout << name + "-grad.nrrd has been exported\n";
        }
        write<value_type, value_type>( name + "-ftle.nrrd", ftle, 1, grid );
        if (verbose) {
            std::cout << name + "-ftle.nrrd has been exported\n";
        }

        if (do_streamlines) {
            export_streamlines(continuous_orbits, name + "-streamlines.vtp");
        }
        else if (do_monitor && !continuous_orbits.empty()) {
            export_streamlines(continuous_orbits, name + "-failed_streamlines.vtp");
        }

        if (do_discrete_orbits) {
            auto bounds = grid.bounds();
            auto size = bounds.size();
            size_t res[2] = { (size_t)(1024*size[0]), (size_t)(1024*size[1]) };
            float* vals = (float*)calloc(res[0]*res[1], sizeof(float));
            srand48(time(0));
            for (size_t i=0; i<discrete_orbits.size(); ++i) {
                double d = drand48();
                const std::vector<position_type>& orb = discrete_orbits[i];
                for (size_t j=0; j<orb.size(); ++j) {
                    if (!bounds.inside(orb[j])) continue;
                    size_t x = (size_t)((orb[j][0]-bounds.min()[0])/size[0]*res[0]);
                    size_t y = (size_t)((orb[j][1]-bounds.min()[1])/size[1]*res[1]);
                    if (x >= res[0] || y >= res[1]) continue;
                    vals[x+y*res[0]] = d;
                }
            }
            spurt::nrrd_utils::writeNrrd(vals, name + "-orbits.nrrd", nrrdTypeFloat, 2, &res[0], false);

            std::vector<position_type> all_pos;
            std::vector<value_type> values;
            for (size_t i=0; i<discrete_orbits.size(); ++i) {
                const std::vector<position_type>& orb = discrete_orbits[i];
                for (size_t j=0; j<orb.size(); ++j) {
                    if (!bounds.inside(orb[j])) continue;
                    all_pos.push_back(orb[j]);
                    values.push_back(i);
                }
            }
            VTK_SMART(vtkPolyData) pdata = vtk_utils::make_points(all_pos);
            vtk_utils::add_vertices(pdata);
            vtk_utils::add_scalars(pdata, values, true, "orbit index");
            vtk_utils::saveVTK(pdata, name + "-orbits.vtp");
        }
        last_iteration += to_do;
    }
}

int main( int argc , const char **argv )
{
    using namespace std;

    namespace xcl = spurt::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map of circular restricted 3-body problem along with FTLE of return map");

    std::array<double, 4> _bounds;
    std::fill(_bounds.begin(), _bounds.end(), invalid_double);

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("output", filename, "Output filename", required_group);
        parser.add_value("mu", mu, mu, "Gravity parameter", optional_group);
        parser.add_value("C", C, C, "Jacobi constant", optional_group);
        parser.add_tuple<4>("bounds", _bounds, "Bounds", optional_group);
        parser.add_value("T", t_max, t_max, "Integration length", optional_group);
        parser.add_value("iter", nb_iterations, nb_iterations, "Number of map iterations", optional_group);
        parser.add_value("delta", iteration_stride, iteration_stride, "Iteration stride", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_value("double", double_prec, double_prec, "Use double precision", optional_group);
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
        parser.add_value("orbit", do_streamlines, do_streamlines, "Export orbits", optional_group);
        parser.add_value("plot", do_discrete_orbits, do_discrete_orbits, "Record Poincare plot", optional_group);
        parser.add_value("debug", do_monitor, do_monitor, "Save failed orbits", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }

    if (std::none_of(_bounds.begin(), _bounds.end(), [&](double d){ return d==invalid_double;})) {
        std::cout << "bounds are valid in input\n";

        std::copy(_bounds.begin(), _bounds.end(), std::ostream_iterator<double>(std::cout, ","));
        std::cout << '\n';

        xmin  = _bounds[0];
        xmax  = _bounds[1];
        xdmin = _bounds[2];
        xdmax = _bounds[3];
    }

    filename=spurt::filename::remove_extension(filename);

    if (double_prec) run< double >(nb_iterations);
    else run< float >(nb_iterations);

    return 0;
}
