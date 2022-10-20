#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include <vexcl/devlist.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/multivector.hpp>
#include <vexcl/mba.hpp>

#include <vtk/vtk_data_helper.hpp>
#include <vtkDataSetWriter.h>

// https://github.com/ddemidov/mba
#include <mba/mba.hpp>

constexpr double _PI_=3.1415926535897932384626;

template <typename T = double>
inline std::array<T, 2> make_array(T x, T y) {
    std::array<T, 2> p = {{x, y}};
    return p;
}

int main(int argc, char *argv[]) {
    const size_t n = argc < 2 ? 1024 * 1024 : std::stoi(argv[1]);
    const size_t m = 100;
    

    // setup the opencl context
    vex::Context ctx(vex::Filter::DoublePrecision);
    std::cout << "Context=" << ctx << '\n';

    vex::profiler<> prof(ctx);

    prof.tic_cpu("generate data");
    std::default_random_engine rng(0);
    std::uniform_real_distribution<double> rnd(-20, 20);

    std::vector< std::array<double, 2> > p(n);
    std::generate(p.begin(), p.end(),
            [&]() {
                return make_array(rnd(rng), rnd(rng));
            });

    std::vector<double> v(n);
    std::transform(p.begin(), p.end(), v.begin(),
            [](const std::array<double,2> &c) {
                double x = c[0];
                double y = c[1];
                double r = sqrt(x*x +y*y);
                return sin(r)/r;
            });

    std::vector< double > x(n);
    std::vector< double > y(n);
    std::vector< double > z(n);
    std::vector< std::array<double, 3> > pos(n);

    std::generate(x.begin(), x.end(), [&]() { return rnd(rng); });
    std::generate(y.begin(), y.end(), [&]() { return rnd(rng); });
    
    x[0] = y[0] = 0.5;
    prof.toc("generate data");

    prof.tic_cpu("GPU");
    {
        prof.tic_cpu("setup");
        vex::mba<2> surf(
                ctx, make_array(-20.01,-20.01), make_array(20.01,20.01),
                p, v, make_array<size_t>(2, 2)
                );
        prof.toc("setup");

        vex::vector<double> X(ctx, x);
        vex::vector<double> Y(ctx, y);
        vex::vector<double> Z(ctx, n);

        prof.tic_cl("interpolate");
        for(size_t i = 0; i < m; ++i)
            Z = surf(X, Y);
        prof.toc("interpolate");
        std::cout << "surf(0.5, 0.5) = " << Z[0] << std::endl;
        
        vex::copy(Z, z);
    }
    prof.toc("GPU");
    
    for (int i=0; i<n; ++i) {
        pos[i] = {{x[i], y[i], z[i]}};
    }
    vtkSmartPointer<vtkPolyData> pd=vtk_utils::make_points(pos);
    vtk_utils::add_vertices(pd);
    vtkSmartPointer<vtkDataSetWriter> writer=vtkSmartPointer<vtkDataSetWriter>::New();
    writer->SetInputData(pd);
    std::string name("blah.vtk");
    if (argc>2) name=argv[2];
    writer->SetFileName(name.c_str());
    writer->Write();
    std::cout << "saved point cloud to " << name << '\n';
    

    prof.tic_cpu("CPU");
    {
        prof.tic_cpu("setup");
        mba::MBA<2> surf(
                make_array(-20.01, -20.01), make_array(20.01, 20.01),
                make_array<size_t>(2, 2), p, v
                );
        prof.toc("setup");

        std::vector<double> z(n);

        prof.tic_cpu("interpolate");
        for(size_t j = 0; j < m; ++j)
#pragma omp parallel for
            for(size_t i = 0; i < n; ++i)
                z[i] = surf(make_array(x[i], y[i]));
        prof.toc("interpolate");
        std::cout << "surf(0.5, 0.5) = " << z[0] << std::endl;
    }
    prof.toc("CPU");

    std::cout << prof << std::endl;
}
