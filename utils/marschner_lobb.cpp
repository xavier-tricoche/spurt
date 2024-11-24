#include <array>
#include <fstream>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <matio.h>
#include <format/filename.hpp>

#include <math/types.hpp>
#include <utils/functions.hpp>
#include <misc/cxxopts.hpp>

using namespace spurt;

int main(int argc, const char* argv[]) {
    
    cxxopts::Options options("marschner_lobb", "Sample Marschner-Lobb function on regular lattice grid and create a random set of samples");
    options.add_options()
        ("r,resolution", "Sampling resolution", cxxopts::value<size_t>())
        ("o,output", "Output files basename", cxxopts::value<std::string>())
        ("h,help", "Print usage information");
    
    auto result = options.parse(argc, argv);

    if (!result.count("resolution") || !result.count("output") || result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    
    size_t n = result["resolution"].as<size_t>();
    std::string filename = result["output"].as<std::string>();
    
    filename = spurt::filename::remove_extension(filename);
    
    double step = 1/(double)(n-1);
    
    std::vector<double> z(n*n, 0);
    for (size_t j=0; j<n; ++j) {
        double y = -1 + j*step;
        for (size_t i=0; i<n; ++i) {
            double x = -1 + i*step;
            z[i+j*n] = spurt::functions2d::z_marschner_lobb(vec2(x,y));
        }
    }
    
    spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&z[0]), filename + ".nrrd", svec2(n), vec2(step), vec2(-1));
    
    std::vector<std::array<double, 3>> vertices(n*n);
    srand48(time(0));
    for (size_t i=0; i<n*n; ++i) {
        double x = -1 + drand48();
        double y = -1 + drand48();
        double f = spurt::functions2d::z_marschner_lobb(vec2(x, y));
        vertices[i][0] = x;
        vertices[i][1] = y;
        vertices[i][2] = f;
    }
    
    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << n*n << " rows\n";
    std::for_each(vertices.begin(), vertices.end(), [&](const std::array<double, 3>& p) {
       output << p[0] << " " << p[1] << " " << p[2] << '\n'; 
    });
    output.close();
    
    return 0;
}