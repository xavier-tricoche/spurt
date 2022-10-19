#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>
#include <misc/time_helper.hpp>

// celltree
#include <vtk/vtk_field.hpp>



typedef spurt::vtk_field<float>  field_type;
typedef field_type::bounds_type   bounds_type;
typedef field_type::point_type    point_type;
typedef field_type::vector_type   vector_type;
typedef field_type::tensor_type   tensor_type;

// parameters
spurt::ivec3   res;
std::string     in_name, out_name;
bounds_type     bounds;

std::string me;
void usage( const std::string& msg="")
{
    if (!msg.empty()) {
        std::cerr << "ERROR: " << msg << std::endl;
    }
    std::cout 
        << "Usage  : " << me << " [parameters] [options]\n"
        << '\n'
        << "Synopsis: resample given dataset using cell tree and native interpolation"
        << '\n'
        << "Parameters:\n"
        << " -i  | --input <string>       Dataset info file\n"
        << " -o  | --output <string>      Output file name\n"
        << " -r  | --resolution <int> x3  Resampling resolution in X, Y, and Z\n"
        << '\n'
        << "Options:\n"
        << " -b  | --bounds <float> x6    Bounding box of region (in world coordinates).\n"
        << "                              Default: all\n"
        << std::endl;
    
    exit(!msg.empty());
}

int main(int argc, char* argv[]) {
    
    res = spurt::ivec3(0);
    
    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                usage("missing input");
            }
            in_name = argv[++i];
        } 
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                usage("missing output");
            }
            out_name = argv[++i];
        }        
        else if (arg == "-h" || arg == "--help") {
            usage();
        }
        else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-3) {
                usage("missing resolution information");
            }
            res[0] = atoi(argv[++i]);
            res[1] = atoi(argv[++i]);
            res[2] = atoi(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-6) {
                usage("missing bounds information");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
            bounds.min()[2] = atof(argv[++i]);
            bounds.max()[2] = atof(argv[++i]);
        }
        else {
            usage("invalid argument: " + arg);
        }
    }
    
    if (in_name.empty() || out_name.empty()) {
        usage("missing input or output file name");
    }
    else if (*std::min_element(res.begin(), res.end()) <= 0) {
        usage("missing / invalid resolution information");
    }
    
    // user reader appropriate for input file type
    spurt::timer timer;
    std::shared_ptr<field_type> field_ptr;
    
    try {
        field_ptr.reset(new field_type(in_name, true, false, false));
    }
    catch( std::exception& e) {
        std::cout << "exception caught while importing " << in_name << ": " << e.what() << '\n';
    }
    field_type& field = *field_ptr;
    
    const bounds_type b = field.bounds();
    
    if (spurt::any(bounds.min() >= bounds.max())) bounds = b;
    std::cout << "bounds = " << bounds << '\n';
    
    point_type step = bounds.size() / point_type(res - spurt::ivec3(1));
    
    size_t number_of_threads = 1;
#ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
#endif
    
    size_t number_of_samples = res[0]*res[1]*res[2];
    
    float *data = (float*)calloc(3*number_of_samples, sizeof(float));
    
    timer.start();
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (size_t n=0 ; n<number_of_samples ; ++n) {
            int i = n % res[0];
            int j = n / res[0];
            int k = j / res[1];
            j = j % res[1];
            point_type p = bounds.min() + step*spurt::fvec3(i, j, k);
            vector_type* v = (vector_type *)(data + 3*n);
            field(p, *v);
            
            int thread = 0;
#if _OPENMP
            thread = omp_get_thread_num();
#endif
            if (!thread) {
                double elapsed = timer.elapsed();
                std::cout << "\r" << n << " / " << number_of_samples
                          << " (" << 100.*(double)n/(double)number_of_samples << "\%) in "
                          << elapsed << " s. (" << (double)n/elapsed << " Hz)         \r"
                          << std::flush;
            }
        }
    }
    std::cout << '\n';
    
    // NRRD file storage
    spurt::nrrd_params<float, 4> params;
    params.mins()[0] = 0;
    params.mins()[1] = bounds.min()[0];
    params.mins()[2] = bounds.min()[1];
    params.mins()[3] = bounds.min()[2];
    params.spacings()[0] = 1;
    params.spacings()[1] = step[0];
    params.spacings()[2] = step[1];
    params.spacings()[3] = step[2];
    params.sizes()[0] = 3;
    params.sizes()[1] = res[0];
    params.sizes()[2] = res[1];
    params.sizes()[3] = res[2];
    params.labels()[0] = "Vx;Vy;Vz";
    params.labels()[1] = "x";
    params.labels()[2] = "y";
    params.labels()[3] = "z";
    
    spurt::writeNrrdFromParams(data, out_name, params);
    std::cout << out_name << " exported\n";
    return 0;
}
