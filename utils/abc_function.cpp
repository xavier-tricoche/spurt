#include <array>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>
#include <utils/functions.hpp>
#include <math/types.hpp>
#include <misc/cxxopts.hpp>


using namespace spurt;

int main(int argc, const char* argv[]) {
    cxxopts::Options options("abc_function", "Sample the ABC function over a 3d raster grid");
    options.add_options()
        ("r,resolution", "Sampling resolution", cxxopts::value<size_t>())
        ("o,output", "Output filename", cxxopts::value<std::string>())
        ("h,help", "Print usage");
    
    auto result = options.parse(argc, argv);

    if (!result.count("resolution") || !result.count("output") || result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    
    size_t n = result["resolution"].as<size_t>();
    std::string filename = result["output"].as<std::string>();
    
    filename = filename::remove_extension(filename);
    
    double step = spurt::TWO_PI/(double)(n-1);
    
    const double __nan__ = std::numeric_limits<double>::quiet_NaN();
    
    functions3d::ABC_field abc;
    
    std::vector<double> f(3*n*n*n, 0);
    for (size_t k=0; k<n; ++k) {
        double z = -spurt::PI + k*step;
        for (size_t j=0; j<n; ++j) {
            double y = -spurt::PI + j*step;
            for (size_t i=0; i<n; ++i) {
                double x = -spurt::PI + i*step;
                vec3 tmp = abc.evaluate(vec3(x,y,z));
                size_t id = 3*(i + n*(j + n*k));
                f[id] = tmp[0];
                f[id+1] = tmp[1];
                f[id+2] = tmp[2];
            }
        }
    }
    
    std::array<size_t, 4> dims({3, n, n, n});
    std::array<double, 4> spcs({__nan__, step, step, step});
    std::array<double, 4> mins({__nan__, -spurt::PI, -spurt::PI, -spurt::PI});
    nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&f[0]), filename + ".nrrd", dims, spcs, mins);
    
    srand48(time(0));
    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << n*n*n << " rows - \"x y z fx(x,y,z) fy(f,y,z) fz(x,y,z)\"\n";
    for (size_t i=0; i<n*n; ++i) {
        double x = -spurt::PI + spurt::TWO_PI*drand48();
        double y = -spurt::PI + spurt::TWO_PI*drand48();
        double z = -spurt::PI + spurt::TWO_PI*drand48();
        vec3 f = abc.evaluate(vec3({x,y,z}));
        output << x << " " << y << " " << z << " " << f[0] << " " << f[1] << " " << f[2] << '\n';
    }
    output.close();
    
    return 0;
}