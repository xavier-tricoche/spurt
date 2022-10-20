#include <array>
#include <fstream>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <matio.h>
#include <format/filename.hpp>


constexpr double PI = 3.14159265358979323846264338;

double z_marschner_lobb(double x, double y) {
    double r = sqrt(x*x + y*y);
    return 2/PI*asin(cos(12*PI*cos(PI*r/2))/4);
}

int main(int argc, char* argv[]) {
    
    size_t n = std::atoi(argv[1]);
    std::string filename = argv[2];
    
    filename = spurt::filename::remove_extension(filename);
    
    double step = 1/(double)(n-1);
    
    std::vector<double> z(n*n, 0);
    for (size_t j=0; j<n; ++j) {
        double y = -1 + j*step;
        for (size_t i=0; i<n; ++i) {
            double x = -1 + i*step;
            z[i+j*n] = z_marschner_lobb(x,y);
        }
    }
    
    std::array<size_t, 2> dims({n, n});
    std::array<double, 2> spcs({step, step});
    std::array<double, 2> mins({-1, -1});
    spurt::writeNrrdFromContainers(reinterpret_cast<double *>(&z[0]), filename + ".nrrd", dims, spcs, mins);
    
    std::vector<std::array<double, 3>> vertices(n*n);
    srand48(time(0));
    for (size_t i=0; i<n*n; ++i) {
        double x = -1 + drand48();
        double y = -1 + drand48();
        double f = z_marschner_lobb(x, y);
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