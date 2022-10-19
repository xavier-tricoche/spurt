/*
 * Copyright 2012 Karsten Ahnert
 * Copyright 2013 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#include <util/timer.hpp>
#include <misc/meta_utils.hpp>
#include <misc/option_parse.hpp>
#include <data/raster.hpp>
#include <format/filename.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

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

namespace odeint = boost::numeric::odeint;

//[ vexcl_state_types

template< typename T >
struct value_traits {
	typedef T                                                 value_type;
	typedef vex::vector< value_type >                         cl_array_type; // 1D array
	typedef vex::multivector< value_type, 3 >                 cl_state_type; // Nx3 2D array
	typedef std::array<value_type, 3>                         array_type;
	typedef Eigen::Matrix<value_type, 3, 1>                   column_type;
	typedef Eigen::Matrix<value_type, 3, 1>                   vector_type;
	typedef Eigen::Matrix<value_type, 3, 1>                   state_type;
	typedef Eigen::Matrix<value_type, 3, 3>                   matrix_type;
	typedef xavier::raster_grid<3, value_type>                grid_type;
	typedef typename grid_type::point_type                    position_type;
	typedef nvis::bounding_box< position_type >               bounds_type;
	typedef typename grid_type::coord_type                    coordinates_type;
	typedef xavier::raster_data<vector_type, 3, value_type>   vector_raster_type;
};
//]

constexpr double _PI_=3.1415926535897932384626433;

//[ vexcl_system
std::array<double, 3> abc {sqrt(3.0), sqrt(2.0), 1.0};

std::string  filename;
size_t       res=64;
bool         verbose = false;
double       dt = 1.0e-2;
double       t_max = 10.0;
double       eps = 1.0e-6;
bool         use_gpu = true;
double       span = 2*_PI_;
bool         double_prec = true;


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

template<typename T>
struct ABC_rhs
{
    typedef value_traits< T > types;
    typedef typename types::value_type value_t;
    typedef typename types::cl_state_type cl_state_type;

    const value_t& A;
    const value_t& B;
    const value_t& C;
    ABC_rhs( const value_t& _A, const value_t& _B, const value_t& _C )
    	: A(_A), B(_B), C(_C) {}

    void operator()( const cl_state_type &x , cl_state_type &dxdt , value_t t ) const
    {
        dxdt(0) = A*sin(x(2)) + C*cos(x(1));
        dxdt(1) = B*sin(x(0)) + A*cos(x(2));
        dxdt(2) = C*sin(x(1)) + B*cos(x(0));
    }
};

template<typename T>
struct raster_wrapper {
    typedef value_traits< T > types;
    typedef typename types::vector_raster_type raster_t;
    typedef typename types::vector_type        vec_t;
    typedef typename types::coordinates_type   coord_t;

    raster_wrapper(const raster_t& r)
        : _raster(r), _res(r.grid().resolution()) {}

    const vec_t& operator()(long int i, long int j, long int k) const {
        coord_t c(i,j,k);
        c+=_res;
        return _raster(c[0]%_res[0], c[1]%_res[1], c[2]%_res[2]);
    }

    const raster_t& _raster;
    const coord_t _res;
};


template<typename T, size_t N>
void gradient(const typename value_traits<T>::vector_raster_type& fmap,
              const std::array<typename value_traits<T>::value_type, N>& kernel,
              std::vector<typename value_traits<T>::matrix_type>& grad,
              std::vector<typename value_traits<T>::value_type>& ftle)
{

	typedef value_traits<T> types;
	typedef typename types::matrix_type mat_t;

    raster_wrapper<T> raster(fmap);
    grad.resize(fmap.size());
    ftle.resize(fmap.size());

    auto res = fmap.grid().resolution();
    auto spc = fmap.grid().spacing();
    int shift = -static_cast<long int>(N/2);

    size_t npoints = fmap.grid().size();

#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif

    #pragma omp parallel
    for (size_t l=0 ; l<npoints; ++l) {
        mat_t m = mat_t::Zero();
        nvis::ivec3 coords = fmap.grid().coordinates(l);
        const size_t& i = coords[0];
        const size_t& j = coords[1];
        const size_t& k = coords[2];
        long int ii, jj, kk;
        // enforce periodicity of phase space
        for (int n=0; n<N; ++n) {
            ii = i+shift+n;
            if (ii >= res[0])
                ii -= res[0];
            else if (ii <= 0)
                ii += res[0];
            jj = j+shift+n;
            if (jj >= res[1])
                jj -= res[1];
            else if (jj <= 0)
                jj += res[1];
            kk = k+shift+n;
            if (kk >= res[2])
                kk -= res[2];
            else if (kk <= 0)
                kk += res[2];

           m.col(0) += kernel[n]*raster(ii, j, k);
           m.col(1) += kernel[n]*raster(i, jj, k);
           m.col(2) += kernel[n]*raster(i, j, kk);
        }
        m.col(0) /= spc[0];
        m.col(1) /= spc[1];
        m.col(2) /= spc[2];
        grad[l] = m;
        if (false && verbose) {
           std::cout << "grad[" << fmap.grid().index(i,j,k) << "]=" << m << '\n';
        }

        Eigen::JacobiSVD<mat_t,Eigen::NoQRPreconditioner> svd(m);
        ftle[l]=2*log(svd.singularValues()[0]);
    }
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
       << "dimension: " << (valsize > 1 ? 4 : 3) << '\n';
	os << "sizes:";
	if (valsize > 1) os << " " << valsize;
    for (int i=0; i<3; ++i) os << " " << grid.resolution()[i];
    os << "\nspacings:";
    if (valsize > 1) os << " nan";
    for (int i=0; i<3; ++i) os << " " << grid.spacing()[i];
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

template <typename T, typename Filter_>
void run(const Filter_ filter) {
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

    // set up number of system, time step and integration time
    const size_t n = res*res*res;
    const value_t step = span/res; // [0, span)

    bounds_t bounds(pos_t(0), pos_t((res-1)*step));
    grid_t grid(coord_t(res, res, res), bounds);

    // initialize the state of the ABC system
    std::vector< pos_t > ic(n);
    for (int i=0; i<n; ++i) {
        int u=i%res;
        int v=(i/res)%res;
        int w=i/(res*res);
        ic[i][0]=step*static_cast<value_t>(u);
        ic[i][1]=step*static_cast<value_t>(v);
        ic[i][2]=step*static_cast<value_t>(w);
    }

    vex::Context ctx(filter);

    // create a stepper
    odeint::runge_kutta_dopri5<cl_state_t, value_t> stepper;

    //std::cout << "starting integration...\n";
    typedef typename std::vector<value_t>::const_iterator cst_it;

    std::vector<value_t> fmap;
    size_t first=0, size=n;
    nvis::timer timer;

    value_t _eps = static_cast<value_t>(eps);
    array_t _abc = convert<double, value_t>(abc);
    size_t global_size = 3*n*sizeof(value_t);
    size_t available_mem = 0;
    for (size_t i=0; i<ctx.size(); ++i) {
        vex::backend::device d = ctx.device(i);
        available_mem += d.getInfo< CL_DEVICE_MAX_MEM_ALLOC_SIZE >();
    }
    ctx.finish();
    size_t nrounds = global_size/available_mem;
    if (nrounds*available_mem < global_size) ++nrounds;
    size = n/nrounds + 1;

    if (verbose) {
        std::cout << "size of problem: " << human_readable_size(global_size) << '\n';
        std::cout << "available memory on selected devices: " << human_readable_size(available_mem) << '\n';
        std::cout << nrounds << " rounds will be necessary to complete computation\n";
    }

    // determine maximum problem size supported by available devices through trial and error
    while (true) {
        if (verbose) std::cout << "#1 fmap.size()=" << fmap.size() << '\n';
        if (!first) timer.restart();
        vex::Context ctx_loc( filter );
        if (verbose) {
            std::cout << "size=" << size << "\t computing block " << first/size+1
                      << " from " << (n % size ? n/size+1 : n/size) << '\n';
			std::cout << "Required memory for next allocation: "
				      << human_readable_size(3*size*sizeof(value_t)) << '\n';
        }
        size_t actual_size=std::min(size, n-first);
        std::vector<value_t> cpu_x(3*actual_size);
        for (int i=0; i<actual_size; ++i) {
            cpu_x[              i]=ic[first+i][0];
            cpu_x[  actual_size+i]=ic[first+i][1];
            cpu_x[2*actual_size+i]=ic[first+i][2];
        }
        try {
			if (verbose) std::cout << "Current context=\n" << ctx_loc << '\n';
			if (verbose) std::cout << "resizing X to " << actual_size << '\n';
			if (verbose) std::cout << "\t (size=" << human_readable_size(actual_size * 3 * sizeof(value_t)) << ")\n";
            cl_state_t X(ctx_loc, actual_size);
			if (verbose) std::cout << "Copying CPU vector to OpenCL\n";
            vex::copy(cpu_x, X);
			if (verbose) std::cout << "Integrating...\n";
            // solve the system
            integrate_const(make_controlled(_eps, _eps, stepper),
                            ABC_rhs<T>(_abc[0], _abc[1], _abc[2]), X,
                            static_cast<value_t>(0),
                            static_cast<value_t>(t_max),
                            static_cast<value_t>(dt));
            //]
			if (verbose) std::cout << "Integration terminated. Copying OpenCL vector to CPU\n";
            vex::copy(X, cpu_x);
            fmap.insert(fmap.end(), cpu_x.begin(), cpu_x.end());
            if (verbose) std::cout << "#2 fmap.size()=" << fmap.size() << '\n';
            first+=size;
            if (first>=n) break;
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
    write<T, value_t>( filename + "-fmap.nrrd", fmap, 3, grid );
    if (verbose) {
        std::cout << filename + "-fmap.nrrd has been exported\n";
    }

    raster_t raster(grid);
    for (int i=0; i<n; ++i) {
        raster[i][0]=fmap[i];
        raster[i][1]=fmap[i+n];
        raster[i][2]=fmap[i+2*n];
    }
    std::vector<mat_t> grad;
    std::vector<value_t> ftle;
    std::array<value_t, 3> kernel = {-0.5, 0, 0.5};
    std::array<value_t, 2> kernel_border = { -0.5, 0.5 };
    gradient<T, 3>(raster, kernel, grad, ftle);
    write<T, mat_t>( filename + "-grad.nrrd", grad, 9, grid );
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

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("output", filename, "Output filename", required_group);
        parser.add_tuple<3>("ABC", abc, abc, "ABC constants", optional_group);
        parser.add_value("T", t_max, t_max, "Integration length", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
		parser.add_value("double", double_prec, double_prec, "Use double precision", optional_group);
        parser.add_value("res", res, res, "Sampling resolution", optional_group);
        parser.add_value("gpu", use_gpu, use_gpu, "Use only GPU for computation", optional_group);
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
        if (use_gpu) run< double >( vex::Filter::DoublePrecision && vex::Filter::GPU );
        else run< double >( vex::Filter::DoublePrecision );
    }
    else {
        if (use_gpu) run< float >(vex::Filter::GPU);
        else run< float >(vex::Filter::Any);
    }

    return 0;
}
