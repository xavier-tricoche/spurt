#include <array>
#include <fstream>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>

#include <boost/math/special_functions/sinc.hpp>


typedef std::complex<double> complex_t;

constexpr double PI = 3.14159265358979323846264338;

double sinc(double x, double y) {
    double r = sqrt(x*x + y*y);
    if (r==0) return 1;
    else return sin(r)/r;
}

complex_t complex_sinc(double x, double y) {
    complex_t z(x, y);
    return boost::math::sinc_pi(z);
}

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
            z[i+j*n] = sinc(x,y);
            re[i+j*n] = std::real(complex_sinc(x,y));
            im[i+j*n] = std::imag(complex_sinc(x,y));
        }
    }
    
    std::array<size_t, 2> dims({n, n});
    std::array<double, 2> spcs({step, step});
    std::array<double, 2> mins({-15, -15});
    spurt::writeNrrdFromContainers(reinterpret_cast<double *>(&z[0]), filename + ".nrrd", dims, spcs, mins);
    spurt::writeNrrdFromContainers(reinterpret_cast<double *>(&re[0]), filename + "-real.nrrd", dims, spcs, mins);
    spurt::writeNrrdFromContainers(reinterpret_cast<double *>(&im[0]), filename + "-imag.nrrd", dims, spcs, mins);
    
    std::vector<std::array<double, 3>> vertices(n*n), real(n*n), imag(n*n);
    srand48(time(0));
    for (size_t i=0; i<n*n; ++i) {
        double x = -15 + 15*drand48();
        double y = -15 + 15*drand48();
        double f = 5*sinc(x, y);
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
    
    srand48(time(0));
    for (size_t i=0; i<n*n; ++i) {
        double x = -5 + 10.*drand48();
        double y = -5 + 10.*drand48();
        real[i][0] = x;
        real[i][1] = y;
        real[i][2] = std::real(complex_sinc(x,y));
        imag[i][0] = x;
        imag[i][1] = y;
        imag[i][2] = std::imag(complex_sinc(x,y));
    }
    
    output.open(filename + "-real.xyz", std::ios::out);
    output << "#" << n*n << " rows\n";
    std::for_each(real.begin(), real.end(), [&](const std::array<double, 3>& p) {
       output << p[0] << " " << p[1] << " " << p[2] << '\n'; 
    });
    output.close();
    
    output.open(filename + "-imag.xyz", std::ios::out);
    output << "#" << n*n << " rows\n";
    std::for_each(imag.begin(), imag.end(), [&](const std::array<double, 3>& p) {
       output << p[0] << " " << p[1] << " " << p[2] << '\n'; 
    });
    output.close();
    
    return 0;
}