#include <array>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>

typedef double value_t;

constexpr value_t PI=3.14159265358979323844;
constexpr value_t TWO_PI=6.28318530717958647688;

typedef std::array<double, 3> state_t;

struct ABC_field {
    ABC_field(value_t a=sqrt(3), value_t b=sqrt(2), value_t c=1) : A(a), B(b), C(c) {}
        
    state_t evaluate(const state_t& x) const {
        state_t dxdt;
        dxdt[0] = A*sin(x[2]) + C*cos(x[1]);
        dxdt[1] = B*sin(x[0]) + A*cos(x[2]);
        dxdt[2] = C*sin(x[1]) + B*cos(x[0]);
        return dxdt;
    }
    
    value_t A, B, C;
};

int main(int argc, char* argv[]) {
    size_t n = std::atoi(argv[1]);
    std::string filename = argv[2];
    
    filename = spurt::filename::remove_extension(filename);
    
    double step = TWO_PI/(double)(n-1);
    
    const double __nan__ = std::numeric_limits<double>::quiet_NaN();
    
    ABC_field abc;
    
    std::vector<double> f(3*n*n*n, 0);
    for (size_t k=0; k<n; ++k) {
        double z = -PI + k*step;
        for (size_t j=0; j<n; ++j) {
            double y = -PI + j*step;
            for (size_t i=0; i<n; ++i) {
                double x = -PI + i*step;
                state_t tmp = abc.evaluate(state_t({x,y,z}));
                size_t id = 3*(i + n*(j + n*k));
                f[id] = tmp[0];
                f[id+1] = tmp[1];
                f[id+2] = tmp[2];
            }
        }
    }
    
    std::array<size_t, 4> dims({3, n, n, n});
    std::array<double, 4> spcs({__nan__, step, step, step});
    std::array<double, 4> mins({__nan__, -PI, -PI, -PI});
    spurt::writeNrrdFromContainers(reinterpret_cast<double *>(&f[0]), filename + ".nrrd", dims, spcs, mins);
    
    srand48(time(0));
    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << n*n*n << " rows - \"x y z fx(x,y,z) fy(f,y,z) fz(x,y,z)\"\n";
    for (size_t i=0; i<n*n; ++i) {
        double x = -PI + TWO_PI*drand48();
        double y = -PI + TWO_PI*drand48();
        double z = -PI + TWO_PI*drand48();
        state_t f = abc.evaluate(state_t({x,y,z}));
        output << x << " " << y << " " << z << " " << f[0] << " " << f[1] << " " << f[2] << '\n';
    }
    output.close();
    
    return 0;
}