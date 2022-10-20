#include <vtk/vtk_utils.hpp>
#include <misc/option_parse.hpp>
#include <Eigen/Core>
#include <limits>
#include <iostream>
#include <math/bounding_box.hpp>


struct Vector : public Eigen::Matrix<double, 3, 1> {
    typedef Eigen::Matrix<double, 3, 1> base_t;
    
    Vector() : base_t() {}
    
    Vector(double a, double b, double c) : base_t() {
        base_t::operator()(0) = a;
        base_t::operator()(1) = b;
        base_t::operator()(2) = c;
    }
    
    Vector(const base_t& other) : base_t(other) {}
    
    static size_t size() { return 3; }
    double operator[](int i) const { return base_t::operator()(i); }
    double& operator[](int i) { return base_t::operator()(i); }
};

typedef Eigen::Matrix<double, 3, 3> matrix_t;
typedef Vector vector_t;

int main(int argc, const char* argv[]) {
    
    size_t N=1000; // number of points
    std::vector<vector_t> points;
    std::vector<vector_t> vectors;
    std::vector<double> bounds(6);
    bool verbose = false;
    
    std::fill(bounds.begin(), bounds.end(), 0.);
    
    namespace xcl = spurt::command_line;
    
    xcl::option_traits 
        required_group(true, false, "Required Options"), 
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");
        
    xcl::option_parser parser(argv[0],
        "Testing pointwise Jacobian computation in VTK unstructured grids");
        
    try {
        // parser.use_default_symbol();
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("number", N, N, "Number of points", optional_group);
        parser.add_tuple<6>("bounds", bounds, "Sampling bounds", optional_group);
        parser.add_flag("verbose", verbose, "Verbose", optional_group);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    nvis::bbox3 domain;
    if (*std::min_element(bounds.begin(), bounds.end()) == 
        *std::max_element(bounds.begin(), bounds.end())) {
        domain.min() = nvis::vec3(-1, -1, -1);
        domain.max() = nvis::vec3(1, 1, 1);  
    }
    else {
        domain.min() = nvis::vec3(bounds[0], bounds[2], bounds[4]);
        domain.max() = nvis::vec3(bounds[1], bounds[3], bounds[5]);
    }
    
    if (verbose) {
        std::cout << "sampling domain: " << domain << '\n';
    }
    
    points.resize(N);
    vectors.resize(N);
    
    srand((unsigned int)time(0));
    matrix_t A = matrix_t::Random();
    std::cout << "A=" << A << '\n';
    Eigen::Vector3d blah = Eigen::Vector3d::Random();
    Eigen::Vector3d blih = A*blah;
    
    srand48(time(0));
    for (size_t i=0; i<N; ++i) {
        double x = domain.min()[0] + drand48()*(domain.max()[0]-domain.min()[0]);
        double y = domain.min()[1] + drand48()*(domain.max()[1]-domain.min()[1]);
        double z = domain.min()[2] + drand48()*(domain.max()[2]-domain.min()[2]);
        
        vector_t p(x,y,z);
        vector_t v = Eigen::Vector3d(A*p);
        points[i] = p;
        vectors[i] = v;
    }
    
    vtkSmartPointer<vtkPolyData> pdata = vtk_utils::make_points(points);
    vtk_utils::add_vectors(pdata, vectors);
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtk_utils::create_mesh3d(pdata);
    
    
    vtk_utils::add_jacobian(grid);
    vtkDoubleArray* tensors = vtkDoubleArray::SafeDownCast(grid->GetPointData()->GetTensors());
    double* t;
    
    double total_error=0;
    double max_error=std::numeric_limits<double>::min();
    for (size_t i=0; i<N; ++i) {
        t = tensors->GetTuple(i);
        double err=0;
        for (int r=0; r<3; ++r) {
            for (int c=0; c<3; ++c) {
                double delta = t[c+3*r]-A(r,c);
                err += delta*delta;
            }
        }
        total_error += err;
        if (err > max_error) max_error = err;
    }
    
    std::cout << "mean error=" << total_error/A.norm()/(double)N << '\n';
    std::cout << "max error=" << max_error/A.norm() << '\n';
    
    return 0;
}
