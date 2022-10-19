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

#include <boost/numeric/odeint.hpp>
//[ vexcl_includes
#include <boost/numeric/odeint/external/vexcl/vexcl.hpp>
//]

#include <opencl/cl_error_codes.hpp>

namespace odeint = boost::numeric::odeint;

typedef double value_t;

//[ vexcl_state_types
typedef vex::vector< value_t >               vector_type;
typedef vex::multivector< value_t, 3 >        state_type;
typedef std::array<value_t, 3>                      pos3;
typedef Eigen::Matrix<value_t, 3, 1>                vec3;
typedef Eigen::Matrix<value_t, 3, 3>                mat3;
typedef xavier::raster_grid<3, double>         grid_type;
typedef xavier::raster_data<vec3, 3, double> raster_type;
typedef grid_type::coord_type                 coord_type;
//]

constexpr value_t _PI_=3.1415926535897932384626433;

//[ vexcl_system
value_t A = sqrt(3.0);
value_t B = sqrt(2.0);
value_t C = 1.0;
size_t res=64;
std::string filename;
bool verbose=false;

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

struct catseye_rhs {
    const value_t& c;
    catseye_rhs(const value_t& _c) : c(_c) {}
    
    void operator()(const state_type& x, state_type& dxdt, value_t t) const {
        // value_t denom=1/(c*cosh(x(1))+sqrt(c*c-1)*cos(x(0)));
        dxdt(0)=c*sinh(x(1))/(c*cosh(x(1))+sqrt(c*c-1)*cos(x(0)));
        dxdt(1)=sqrt(c*c-1)*sin(x(0))/(c*cosh(x(1))+sqrt(c*c-1)*cos(x(0)));
        dxdt(2)=1/(c*cosh(x(1))+sqrt(c*c-1)*cos(x(0)));
    } 
};

struct raster_wrapper {
    raster_wrapper(const xavier::raster3d<vec3>& r) 
        : _raster(r), _res(r.grid().resolution()) {}
    
    const vec3& operator()(long int i, long int j, long int k) const {
        nvis::fixed_vector<long int, 3> c(i,j,k);
        c+=_res;
        return _raster(c[0]%_res[0], c[1]%_res[1], c[2]%_res[2]);
    }
    
    const xavier::raster3d<vec3>& _raster;
    const xavier::raster3d<vec3>::grid_type::coord_type _res;
};
    
    
template<size_t N>
void gradient(const xavier::raster3d<vec3>& fmap,
              const std::array<value_t, N>& kernel,
              std::vector<mat3>& grad,
              std::vector<value_t>& ftle) {
                    
    raster_wrapper raster(fmap);    
    grad.resize(fmap.size());
    ftle.resize(fmap.size());
        
    auto res=fmap.grid().resolution();
    auto spc=fmap.grid().spacing();
    int shift=-static_cast<long int>(N/2);
    
    for (int i=0; i<res[0]; ++i) {
        for (int j=0; j<res[1]; ++j) {
            for (int k=0; k<res[2]; ++k) {
               mat3 m=mat3::Zero();
               for (int n=0; n<N; ++n) {
                   m.col(0)+=kernel[n]*raster(i+shift+n, j,         k);
                   m.col(1)+=kernel[n]*raster(i,         j+shift+n, k);
                   m.col(2)+=kernel[n]*raster(i,         j,         k+shift+n);
               }
               m.col(0)/=spc[0];
               m.col(1)/=spc[1];
               m.col(2)/=spc[2];
               grad[fmap.grid().index(i,j,k)]=m;
               if (false && verbose) {
                   std::cout << "grad[" << fmap.grid().index(i,j,k) << "]=" << m << '\n';
               }
               
               Eigen::JacobiSVD<mat3,Eigen::NoQRPreconditioner> svd(m);
               ftle[fmap.grid().index(i,j,k)]=2*log(svd.singularValues()[0]);
            }
        }
    }
    
    if (verbose) {
        // srand48(time());
        for (int i=0; i<20; ++i) {
            std::cout << "grad[" << i << "]=" << grad[i] << '\n';
        }
    }
}

    
//]

template<typename Value_>
void write(const std::string& filename, const std::vector<Value_>& data,
           size_t valsize, const grid_type& grid) {
    std::ofstream outf(filename.c_str(), std::ofstream::binary);
    std::ostringstream os;
    os << "NRRD0001\n"
       << "# Complete NRRD file format specification at:\n"
       << "# http://teem.sourceforge.net/nrrd/format.html\n"
       << "type: " << xavier::type2string<value_t>::type_name() << "\n"
       << "dimension: ";
    if (valsize>1) os << " " << valsize;
    for (int i=0; i<3; ++i)  os << " " << grid.resolution()[i];
    os << "\nspacings:";
    if (valsize>1) os << " nan";
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

int main( int argc , const char **argv )
{
    using namespace std;
    using namespace odeint;
    
    namespace xcl = xavier::command_line;
    
    std::array<value_t, 3> abc {sqrt(3), sqrt(2), 1};
    value_t dt = 1.0e-2;
    value_t t_max = 10.0;
    value_t eps=1.0e-6;
    value_t ymax=2.5;
    bool use_cpu=true;
    value_t span=2*_PI_;
    
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map of ABC flow using OpenCL");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("output", filename, "Output filename", required_group);
        parser.add_value("ymax", ymax, ymax, "Half height of sampling domain", optional_group);
        parser.add_value("T", t_max, t_max, "Integration length", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_value("res", res, res, "Sampling resolution", optional_group);
        parser.add_value("cpu", use_cpu, use_cpu, "Use host", optional_group);
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
    vex::Context ctx(vex::Filter::DoublePrecision);
    if (verbose) std::cout << "Context=" << ctx << '\n';
    
    // set up number of system, time step and integration time
    const size_t n = res*res*res;
    const value_t step=span/res; // [0, span)
    
    nvis::bbox3 bounds;
    bounds.min()[0]=0;
    bounds.max()[0]=span;
    bounds.min()[1]=bounds.min()[2]=-ymax;
    bounds.max()[1]=bounds.max()[2]=ymax;
    
    xavier::rgrid3d grid(coord_type(res, res, res), bounds);

    // initialize the state of the cat's eye system
    std::vector< vec3 > ic(n);
    for (int i=0; i<n; ++i) {
        int u=i%res;
        int v=(i/res)%res;
        int w=i/(res*res);
        ic[i][0]=step*static_cast<value_t>(u);
        ic[i][1]=step*static_cast<value_t>(v);
        ic[i][2]=step*static_cast<value_t>(w);
    }
    state_type X(ctx, n);

    // create a stepper
    runge_kutta_dopri5<state_type> stepper;

    std::cout << "starting integration...\n";
    typedef std::vector<value_t>::const_iterator cst_it;

    std::vector<value_t> fmap;
    size_t first=0, size=n;
    nvis::timer timer;
    while (true) {
        //std::cout << "#1 fmap.size()=" << fmap.size() << '\n';
        if (!first) timer.restart();
        vex::Context ctx_loc( vex::Filter::DoublePrecision );
        if (verbose) {
            std::cout << "size=" << size << "\t computing block " << first/size+1
                      << " from " << (n % size ? n/size+1 : n/size) << '\n'; 
        }
        size_t actual_size=std::min(size, n-first);
        std::vector<value_t> cpu_x(3*actual_size);
        for (int i=0; i<actual_size; ++i) {
            cpu_x[              i]=ic[first+i][0];
            cpu_x[  actual_size+i]=ic[first+i][1];
            cpu_x[2*actual_size+i]=ic[first+i][2];
        }
        try {
            X.resize(ctx_loc, actual_size);
            vex::copy(cpu_x, X);
            // solve the system
            integrate_const(make_controlled(eps, eps, stepper), 
                            catseye_rhs(2), X, 
                            static_cast<value_t>(0), t_max, dt);
            vex::copy(X, cpu_x);
            fmap.insert(fmap.end(), cpu_x.begin(), cpu_x.end());
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
            continue;
        }
    }

    std::cout << "integration completed in " << timer.elapsed() << "s.\n";
    write(filename+"-fmap.nrrd", fmap, 3, grid);
    if (verbose) {
        std::cout << filename+"-fmap.nrrd has been exported\n";
    }
    
    return 0;
}
