#include <array>
#include <fstream>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>
#include <math/types.hpp>

#include <utils/functions.hpp>


using namespace spurt;

int main(int argc, char* argv[]) {
    
    size_t n = std::atoi(argv[1]);
    std::string filename = argv[2];
    
    filename = spurt::filename::remove_extension(filename);
    
    double step = 15/(double)(n-1);
    
    std::vector<double> z(n*n, 0);
    std::vector<double> re(n*n, 0);
    std::vector<double> im(n*n, 0);
    for (size_t j=0; j<n; ++j) {
        double y = -15 + j*step;
        for (size_t i=0; i<n; ++i) {
            double x = -15 + i*step;
            const vec2 p(x,y);
            z[i+j*n] = functions2d::sinc(p);
            re[i+j*n] = functions2d::complex_sinc(p).real();
            im[i+j*n] = functions2d::complex_sinc(p).imag();
        }
    }
    
    spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&z[0]),  filename + ".nrrd",      svec2(n), vec2(step), vec2(-15));
    spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&re[0]), filename + "-real.nrrd", svec2(n), vec2(step), vec2(-15));
    spurt::nrrd_utils::writeNrrdFromContainers(reinterpret_cast<double *>(&im[0]), filename + "-imag.nrrd", svec2(n), vec2(step), vec2(-15));
    
    std::vector<vec3> vertices(n*n), real(n*n), imag(n*n);
    srand48(time(0));
    for (size_t i=0; i<n*n; ++i) {
        double x = -15 + 15*drand48();
        double y = -15 + 15*drand48();
        double f = 5*functions2d::sinc(vec2(x, y));
        vertices[i] = vec3(x,y,f);
    }
    
    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << n*n << " rows\n";
    std::for_each(vertices.begin(), vertices.end(), [&](const vec3& p) {
       output << p << '\n'; 
    });
    output.close();
    
    srand48(time(0));
    for (size_t i=0; i<n*n; ++i) {
        double x = -5 + 10.*drand48();
        double y = -5 + 10.*drand48();
        std::complex<double> f = spurt::functions2d::complex_sinc(vec2(x,y));
        real[i] = vec3(x, y, f.real());
        imag[i] = vec3(x, y, f.imag());
    }
    
    output.open(filename + "-real.xyz", std::ios::out);
    output << "#" << n*n << " rows\n";
    std::for_each(real.begin(), real.end(), [&](const vec3& p) {
       output << p << '\n'; 
    });
    output.close();
    
    output.open(filename + "-imag.xyz", std::ios::out);
    output << "#" << n*n << " rows\n";
    std::for_each(imag.begin(), imag.end(), [&](const vec3& p) {
       output << p << '\n'; 
    });
    output.close();
    
    return 0;
}