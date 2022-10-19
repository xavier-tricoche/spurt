/*
 * Copyright 2012 Karsten Ahnert
 * Copyright 2013 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#define VEXCL_SHOW_KERNELS
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
#include <VTK/vtk_io_helper.hpp>
#include <VTK/vtk_data_helper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/dopri5.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <vexcl/vexcl.hpp>
#include <vexcl/devlist.hpp>

#include <boost/numeric/odeint.hpp>
//[ vexcl_includes
#include <boost/numeric/odeint/external/vexcl/vexcl.hpp>
//]

#include <opencl/cl_error_codes.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef NO_TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
// #include <tbb/mutex.h>
#endif

#include <mutex>

namespace odeint = boost::numeric::odeint;

#ifndef NO_TBB
std::mutex output_mutex;
std::mutex progress_mutex;
std::mutex orbit_mutex;
std::atomic<size_t> progress_counter;
#else
size_t progress_counter;
#endif

void write_to_ostream(std::ostream& os, const std::string& str) {
    {
    #ifndef  NO_TBB
        std::scoped_lock lock(output_mutex);
    #endif
        os << str << '\n';
    }
}

void update_progress(xavier::ProgressDisplay& progress) {
    {
    #ifndef  NO_TBB
        std::scoped_lock lock(progress_mutex);
    #endif
        progress.update(progress_counter);
    }
}

template< typename T >
void append_orbit(std::vector< std::vector<T> >& out, const std::vector<T>& in) {
    {
    #ifndef  NO_TBB
        std::scoped_lock lock(orbit_mutex);
    #endif
        out.push_back(in);
    }
}

//[ vexcl_state_types

template< typename T >
struct value_traits {
    typedef T                                                 value_type;
    // array of scalars
    typedef vex::vector< value_type >                         cl_scalar_type;
    // array of 4D vectors for integration (x, y, xdot, ydot)
    typedef vex::multivector< value_type, 4>                  cl_state_type;
    typedef std::array<value_type, 2>                         array_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   column_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   vector_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   state_type;
    typedef Eigen::Matrix<value_type, 2, 2>                   matrix_type;

    // sampling grid for Poincare section
    typedef xavier::raster_grid<2, value_type>                grid_type;
    typedef typename grid_type::point_type                    position_type;
    typedef nvis::bounding_box< position_type >               bounds_type;
    typedef typename grid_type::coord_type                    coordinates_type;
    typedef typename grid_type::point_type                    point_type;
    typedef xavier::raster_data<column_type, 2, value_type>   vector_raster_type;
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
bool        use_gpu = true;
bool        double_prec = true;
int         niter = 10;
int         delta_iter = 1;
bool        compute_ftle = true;
bool        do_streamlines = false;
bool        do_monitor = true;

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

// VEX_FUNCTION(double, rhs, (double, x)(double, y)(double, xd)(double, yd)

template<typename T>
struct pcr3bp_rhs
{
    typedef value_traits< T > types;
    typedef typename types::value_type value_t;
    typedef typename types::cl_state_type cl_state_type;

    const value_t& C;
    const value_t& mu;
    pcr3bp_rhs( const value_t& _mu, const value_t& _C )
        : mu(_mu), C(_C) {}

    void operator()( const cl_state_type &x , cl_state_type &dxdt , value_t t ) const
    {
        using vex::tag;
        auto d = vex::make_temp<1>(sqrt((tag<1>(x(0))+tag<5>(mu))*(tag<1>(x(0))+tag<5>(mu)) + tag<2>(x(1))*tag<2>(x(1))));
        auto r = vex::make_temp<2>(sqrt((tag<1>(x(0))-1+tag<5>(mu))*(tag<1>(x(0))-1+tag<5>(mu)) + tag<2>(x(1))*tag<2>(x(1))));
        auto d3 = vex::make_temp<3>(d*d*d);
        auto r3 = vex::make_temp<4>(r*r*r);
        dxdt = std::make_tuple(
            tag<3>(x(2)), // dxdt = xdot
            tag<4>(x(3)), // dydt = ydot
            // 2 * yd  + x - (mu - 1) * (x + mu) / d3 - mu * (x + mu - 1) / r3;
            2*tag<4>(x(3)) + tag<1>(x(0)) +
                (tag<5>(mu)-1)*(tag<1>(x(0))+tag<5>(mu))/d3 - tag<5>(mu)*(tag<1>(x(0))+tag<5>(mu)-1)/r3,
            tag<2>(x(1)) - 2*tag<3>(x(2)) +
                (tag<5>(mu)-1)*tag<2>(x(1))/d3 - tag<5>(mu)*tag<2>(x(1))/r3
        );
    }
};

template<typename T, typename Vector = Eigen::Matrix<T, 2, 1> >
struct raster_wrapper {
    typedef T                                      scalar_type;
    typedef xavier::raster_grid<2, T>              grid_type;
    typedef xavier::raster_data<Vector, 2, T>      raster_type;
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
void gradient(const xavier::raster_data<Vector, 2, T>& fmap,
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

    xavier::ProgressDisplay progress(true);

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
       << "type: " << xavier::type2string<value_t>::type_name() << "\n"
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

#ifndef NO_POW
            value_type d3 = std::pow((x + m_mu    )*(x + m_mu    ) + y*y, 1.5);
            value_type r3 = std::pow((x + m_mu - 1)*(x + m_mu - 1) + y*y, 1.5);
#else
            value_type d = sqrt((x + m_mu    )*(x + m_mu    ) + y*y);
            value_type r = sqrt((x + m_mu - 1)*(x + m_mu - 1) + y*y);
            value_type d3 = d*d*d;
            value_type r3 = r*r*r;
#endif
            return state_type(
                dx,
                dy,
                2*dy + x    + (m_mu - 1)*(x + m_mu)/d3 - m_mu*(x + m_mu - 1)/r3,
                   y - 2*dx + (m_mu - 1)*     y    /d3 - m_mu*     y        /r3
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
};

void export_streamlines(const std::vector< std::vector< world_pos_type > >& lines,
                        const std::string& filename) {
    std::cout << lines.size() << " streamlines on output\n";
    for (int i=0; i<lines.size(); ++i) {
        std::cout << "streamline #" << i << " contains " << lines[i].size() << " vertices\n";
    }
    std::vector<nvis::ivec2> dummy;
    VTK_SMART(vtkPolyData) pdata = vtk_utils::make_polylines(lines, dummy, 0);
    vtk_utils::saveVTK(pdata, filename);
}

template<typename T = double>
void run_cpu_only(const size_t niterations) {
    typedef T                                                 value_type;
    typedef dopri45_wrapper<value_type>                       wrapper_type;
    typedef typename wrapper_type::state_type                 state_type;
    typedef typename wrapper_type::ode_solver_type            solver_type;
    typedef typename wrapper_type::step_type                  step_type;
    typedef typename wrapper_type::cr3bp                      rhs_type;
    typedef nvis::fixed_vector<value_type, 2>                 position_type;
    typedef Eigen::Matrix<value_type, 2, 1>                   column_type;
    typedef Eigen::Matrix<value_type, 2, 2>                   matrix_type;
    typedef xavier::raster_grid<2, value_type>                grid_type;
    typedef xavier::raster_data<position_type, 2, value_type> vector_raster_type;
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
    std::vector<position_type> map_returns;
    size_t id=0;
    for (size_t i=0; i<n_total; ++i) {
        try {
            position_type p = grid(grid.coordinates(i));
            value_type dydt = yd(p[0], p[1]);
            glob2loc[i] = id++;
            valid_states.push_back(state_type(p[0], 0, p[1], dydt));
            map_returns.push_back(p);
        }
        catch(...) { continue; }
    }

    const size_t n_valid = valid_states.size();
    if (verbose) {
        std::cout << "There are " << n_valid << " valid initial conditions\n";
    }

    rhs_type rhs(C, mu);
    std::vector<bool> stopped(n_valid, false);

    xavier::ProgressDisplay progress(true);

    size_t last_iteration = 0;
    std::vector< std::vector<world_pos_type> > orbits;

    for (size_t iteration=1; iteration <= niter; ++iteration) {
        progress.fraction_on();
        progress_counter = 0;
        if (verbose) std::cout << "Integrating iteration " << iteration << "...\n";
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
                position_type xing;
                value_type xing_time, last_time;
                wrapper_type::integrate(last, last_time, xing, xing_time,
                                        rhs, y0, 0, 0, 0, eps, steps, true);
                valid_states[n] = last;
                map_returns[n] = xing;
                if (do_streamlines) {
                    append_orbit(orbits, steps);
                }
            }
            catch(...) {
                stopped[n] = true;
                if (do_monitor) {
                    append_orbit(orbits, steps);
                }
            }
        }
#ifndef NO_TBB
        });
#endif
        progress.end();
        std::cout << progress << '\n';

        if (iteration - last_iteration == delta_iter) {
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
                    fmap[n] = map_returns[i];
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
                oss << std::setw(5) << std::setfill('0') << iteration;
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
            last_iteration = iteration;

            if (do_streamlines) {
                export_streamlines(orbits, name + "-streamlines.vtp");
            }
            else if (do_monitor && !orbits.empty()) {
                export_streamlines(orbits, name + "-failed_streamlines.vtp");
            }
        }
    }
}

template <typename T, typename Filter_>
void run_gpu(const Filter_ filter) {
    typedef value_traits< T > types;
    typedef typename types::value_type           value_t;
    typedef typename types::bounds_type          bounds_t;
    typedef typename types::position_type        pos_t;
    typedef typename types::cl_state_type        cl_state_t;
    typedef typename types::coordinates_type     coord_t;
    typedef typename types::matrix_type          mat_t;
    typedef typename types::vector_raster_type   raster_t;
    typedef typename types::grid_type            grid_t;
    typedef typename types::array_type           array_t;
    typedef typename types::vector_type          vector_t;

    typedef std::array<value_t, 3> ic_t; // initial conditions

    // set up number of system, time step and integration time
    const size_t Ntotal = res[0]*res[1];

    bounds_t bounds(pos_t(xmin, xdmin), pos_t(xmax, xdmax));
    grid_t grid(coord_t(res[0], res[1]), bounds);
    const auto spacing = grid.spacing();

    std::cout << "grid spacing = " << spacing << '\n';
    std::cout << "grid bounds = " << bounds << '\n';

    // initialize the state of the CR3BP system
    std::map<size_t, size_t> glob2loc;
    size_t id=0;

    std::vector<ic_t> valid_ics;
    for (size_t i=0; i<Ntotal; ++i) {
        try {
            pos_t p = grid(grid.coordinates(i));
            value_t dydt = yd(p[0], p[1]);
            glob2loc[i] = id++;
            valid_ics.push_back(ic_t({p[0], p[1], dydt}));
        }
        catch(...) { continue; }
    }

    if (verbose) {
        std::cout << "There are " << valid_ics.size() << " valid ICs\n";
    }

    const size_t Nactive = valid_ics.size();

    vex::Context ctx(filter);

    // create a stepper
    odeint::runge_kutta_dopri5<cl_state_t, value_t> stepper;

    std::cout << "starting integration...\n";
    typedef typename std::vector<value_t>::const_iterator cst_it;

    std::vector<value_t> fmap;
    size_t first=0, size=Nactive;
    nvis::timer timer;

    value_t _eps = static_cast<value_t>(eps);
    size_t global_size = 4*Nactive*sizeof(value_t);
    size_t available_mem = 0;
    for (size_t i=0; i<ctx.size(); ++i) {
        vex::backend::device d = ctx.device(i);
        available_mem += d.getInfo< CL_DEVICE_MAX_MEM_ALLOC_SIZE >();
    }
    ctx.finish();
    size_t nrounds = 1;

    if (verbose) {
        std::cout << "size of problem: " << human_readable_size(global_size) << '\n';
        std::cout << "available memory on selected devices: " << human_readable_size(available_mem) << '\n';
    }

    // determine maximum problem size supported by available devices through trial and error
    while (true) {
        if (verbose) std::cout << "#1 fmap.size()=" << fmap.size() << '\n';
        if (!first) timer.restart();
        vex::Context ctx_loc( filter );
        if (verbose) {
            std::cout << "size=" << size << "\t computing block " << first/size+1
                      << " from " << (Nactive % size ? Nactive/size+1 : Nactive/size) << '\n';
            std::cout << "Required memory for next allocation: "
                      << human_readable_size(4*size*sizeof(value_t)) << '\n';
        }
        size_t tight_size=std::min(size, Nactive-first);
        std::vector<value_t> cpu_x(4*tight_size);
        for (int i=0; i<tight_size; ++i) { // column first storage
            cpu_x[             i] = valid_ics[first+i][0]; // x
            cpu_x[  tight_size+i] = 0;                     // y
            cpu_x[2*tight_size+i] = valid_ics[first+i][1]; // dxdt
            cpu_x[3*tight_size+i] = valid_ics[first+i][2]; // dydt
        }
        try {
            if (verbose) std::cout << "Current context=\n" << ctx_loc << '\n';
            if (verbose) std::cout << "resizing X to " << tight_size << '\n';
            if (verbose) std::cout << "\t (size=" << human_readable_size(tight_size * 4 * sizeof(value_t)) << ")\n";
            cl_state_t X(ctx_loc, tight_size);
            if (verbose) std::cout << "Copying CPU vector to OpenCL\n";
            vex::copy(cpu_x, X);
            if (verbose) std::cout << "Integrating...\n";
            // solve the system
            integrate_const(make_controlled(_eps, _eps, stepper),
                            pcr3bp_rhs<value_t>(mu, C), X,
                            static_cast<value_t>(0),
                            static_cast<value_t>(t_max),
                            static_cast<value_t>(dt));
            //]
            if (verbose) std::cout << "Integration terminated. Copying OpenCL vector to CPU\n";
            vex::copy(X, cpu_x);
            std::for_each(cpu_x.begin(), cpu_x.end(), [&](value_t v) {
                if (std::isnan(v) || std::isinf(v)) {
                    std::cerr << "value is NaN!!\n";
                }
            });
            for (size_t i=0; i<tight_size; ++i) {
                fmap.push_back(cpu_x[i             ]); // x
                fmap.push_back(cpu_x[i+  tight_size]); // y
                fmap.push_back(cpu_x[i+2*tight_size]); // xdot
            }
            if (verbose) std::cout << "#2 fmap.size()=" << fmap.size() << '\n';
            first+=size;
            if (first>=Nactive) break;
        }
        catch(cl::Error& e) {
            if (verbose) {
                std::cout << "exception caught: " << e.what() << '\n';
                std::cout << "OpenCL error code=" << getErrorString(e.err()) << '\n';
            }
            size/=2;
            first=0;
            fmap.clear();
        }
        ctx_loc.finish();
    }

    std::cout << "integration completed in " << timer.elapsed() << "s.\n";
    // convert flow map from column-major to row-major
    std::vector<value_t> complete_fmap(Ntotal*3, 0);

    std::cout << "map contents:\n";
    for (auto it=glob2loc.begin(); it!=glob2loc.end(); ++it) {
        std::cout << it->first << " -> " << it->second << '\n';
    }

    for (size_t i=0; i<Ntotal; ++i) {
        auto what = glob2loc.find(i);
        if (what != glob2loc.end()) {
            // valid flow map result
            size_t j = what->second;
            complete_fmap[3*i  ] = fmap[3*j  ];
            complete_fmap[3*i+1] = fmap[3*j+1];
            complete_fmap[3*i+2] = fmap[3*j+2];
        }
    }
    write<T, value_t>( filename + "-fmap.nrrd", complete_fmap, 3, grid );
    if (verbose) {
        std::cout << filename + "-fmap.nrrd has been exported\n";
    }

    raster_t raster(grid, reinterpret_cast<const std::vector<vector_t>&>(complete_fmap));
    // for (int i=0; i<Ntotal; ++i) {
    //     raster[i][0]=fmap[3*i];
    //     raster[i][1]=fmap[3*i+1];
    //     raster[i][2]=fmap[3*i+2];
    // }
    std::vector<mat_t> grad;
    std::vector<value_t> ftle;
    std::array<value_t, 3> kernel = {-0.5, 0, 0.5};
    gradient<T, 3>(raster, kernel, grad, ftle);
    write<T, mat_t>( filename + "-grad.nrrd", grad, 6, grid );
    if (verbose) {
        std::cout << filename + "-grad.nrrd has been exported\n";
    }
    write<T, value_t>( filename + "-ftle.nrrd", ftle, 1, grid );
    if (verbose) {
        std::cout << filename + "-ftle.nrrd has been exported\n";
    }
}

int main( int argc , const char **argv )
{
    using namespace std;
    using namespace odeint;

    namespace xcl = xavier::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map of ABC flow using OpenCL");

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
        parser.add_value("iter", niter, niter, "Number of map iterations", optional_group);
        parser.add_value("delta", delta_iter, delta_iter, "Iteration stride", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_value("double", double_prec, double_prec, "Use double precision", optional_group);
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
        parser.add_value("gpu", use_gpu, use_gpu, "Use GPU for computation", optional_group);
        parser.add_value("orbit", do_streamlines, do_streamlines, "Export orbits", optional_group);
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

    filename=xavier::filename::remove_extension(filename);

    //[ vexcl_main
    // setup the opencl context
    if (use_gpu) {
        vex::Context ctx_gpu(vex::Filter::GPU);
        if (!ctx_gpu) {
            std::cerr << "WARNING: no GPU available. Switching to CPU device.\n";
            use_gpu = false;
        }
        ctx_gpu.finish();
        if (use_gpu && double_prec) {
            // check if double precision is supported by GPU
            vex::Context ctx_gpu_dbl(vex::Filter::DoublePrecision && vex::Filter::GPU);
            if (!ctx_gpu_dbl) {
                std::cerr << "WARNING: Chosen GPU device does not support double precision\n";
                double_prec = false;
            }
            ctx_gpu_dbl.finish();
        }
    }

    if (double_prec) {
        if (use_gpu) run_gpu< double >(vex::Filter::DoublePrecision && vex::Filter::GPU);
        else run_cpu_only< double >(niter);
    }
    else {
        if (use_gpu) run_gpu< float >(vex::Filter::GPU);
        else run_cpu_only< float >(niter);
    }

    return 0;
}
