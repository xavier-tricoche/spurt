#include <vector>
#include <initializer_list>
#include <algorithm>

#define TAPKEE_USE_LGPL_COVERTREE
// #define TAPKEE_WITH_VIENNACL
#include <tapkee/tapkee.hpp>
#include <tapkee/defines/methods.hpp>
#include <tapkee/exceptions.hpp>
#include <tapkee/callbacks/precomputed_callbacks.hpp>
#include <tapkee/callbacks/eigen_callbacks.hpp>
#include <tapkee/neighbors/neighbors.hpp>

#include <misc/option_parse.hpp>
#include <graph/connectivity.hpp>

// VTK
#include <vtkColorTransferFunction.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

// VTK helper functions and macros
#include <vtk/vtk_utils.hpp>
#include <graphics/colors.hpp>

#include <misc/progress.hpp>

typedef nvis::fixed_vector<double, 3> vec3;

void swiss_roll(std::vector<vec3>& points,
                std::vector<double>& theta, int N,
                bool regular_sampling=false) {
    std::vector<double> height(N);
    theta.resize(N);
    points.resize(N);

    if (!regular_sampling) {
        srand48(time(NULL));
        std::for_each(theta.begin(), theta.end(), [&](double& t)
            {
                t=3*M_PI/2*(1+2*drand48());
            });
        std::for_each(height.begin(), height.end(), [&](double& h)
            {
                h=drand48()-0.5;
            });
    }
    else {
        int nrows=N/10;

        double dtheta=3.*M_PI/2.*4./static_cast<double>(nrows);
        double t=3.*M_PI/2.;
        double step=t*dtheta;
        std::fill(&theta[0], &theta[10], t);
        for (int i=0; i<nrows; ++i) {
            t+=step/t;
            std::fill(&theta[10*i], &theta[10*(i+1)], t);
        }
        for (int i=0; i<10; ++i) {
            double u=static_cast<double>(i)/10.;
            for (int n=0; n<nrows; ++n) {
                height[10*n+i]=-0.5+u;
            }
        }
    }
    for (int i=0;i<N;++i) {
        points[i][0]=theta[i]*cos(theta[i]);
        points[i][1]=10.*height[i];
        points[i][2]=theta[i]*sin(theta[i]);
    }
}

std::string& lower_case(std::string& str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

bool check_string(std::string arg,
                  const std::initializer_list<std::string>& arg_list) {
    std::vector<std::string> argv(arg_list);
    lower_case(arg);
    std::for_each(argv.begin(), argv.end(), [&](std::string& s)
        { lower_case(s);});

    if (argv.size()==1) {
        return arg==argv[0];
    }
    else {
        std::string str1(argv[0]);
        std::string str2(argv[0]);
        std::string str3(1, argv[0][0]);
        for (size_t i=1; i<argv.size(); ++i) {
            str1+=argv[i];
            str2+='_'+argv[i];
            str3+=argv[i][0];
        }
        return arg==str1 || arg==str2 || arg==str3;
    }
}

void check_connectivity(std::vector<int>& cc_ids,
                        const tapkee::tapkee_internal::Neighbors& neighbors) {
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph;
    for (size_t i=0; i<neighbors.size(); ++i) {
        for (size_t j=0; j<neighbors[i].size(); ++j) {
            boost::add_edge(i, neighbors[i][j], graph);
        }
    }

    size_t ncc=spurt::connected_components(cc_ids, graph);
    std::cout << "There are " << ncc << " connected components\n";
}

int npoints=100;
int n_neighbors=3;
double sphere_radius=0.1;
int verbose=0;
bool show_points=true;
bool show_edges=true;
bool regular_sampling=false;
nvis::vec3 bg_color;
std::string method_name="lle";
std::string neighbors_method_name="brute";

struct distance_callback {
    distance_callback(const std::vector<vec3>& _points) : points(_points) {}

    tapkee::ScalarType distance(tapkee::IndexType a, tapkee::IndexType b) const {
        return nvis::norm(points[a]-points[b]);
    }

    tapkee::ScalarType distance(const vec3& a, const vec3& b) const {
        return nvis::norm(a-b);
    }

    const std::vector<vec3>& points;
};

struct kernel_callback {
    kernel_callback(const std::vector<vec3>& _points) : points(_points) {}

    tapkee::ScalarType kernel(tapkee::IndexType a, tapkee::IndexType b) const {
        return nvis::inner(points[a], points[b]);
    }

    tapkee::ScalarType kernel(const vec3& a, const vec3& b) const {
        return nvis::inner(a, b);
    }

    const std::vector<vec3>& points;
};

bool init(int argc, const char* argv[]) {
    namespace xcl=spurt::command_line;

    xcl::option_traits
        required_group(true, false, "Required parameters"),
        positional_group(true, true, "Positional parameters"),
        optional_group(false, false, "Optional parameters");

    xcl::option_parser parser(argv[0],
        "Test Tapkee implementation using classical swiss roll");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("n", npoints, "Number of points",
                         optional_group);
        parser.add_value("regular", regular_sampling,
                         "Regular manifold sampling", optional_group);
        parser.add_value("neighbors", n_neighbors,
                         "Number of neighbors used to determine connectivity",
                         optional_group);
        parser.add_value("nmethod", neighbors_method_name,
                         "Method used to compute neighborhood graph",
                         optional_group);
        parser.add_value("method", method_name,
                         "Dimensionality reduction method", optional_group);
        parser.add_value("points", show_points,
                        "Show vertices", optional_group);
        parser.add_value("edges", show_edges,
                        "Show edges", optional_group);
        parser.add_value("sphere", sphere_radius,
                        "Radius of sphere glyphs", optional_group);
        parser.add_value("verbose", verbose,
                         "verbose level", optional_group);
        parser.add_tuple<3>("bg", bg_color, bg_color,
                            "Background color", optional_group);
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        return false;
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        return false;
    }

    return true;
}

int main(int argc, const char* argv[]) {
    init(argc, argv);

    tapkee::DimensionReductionMethod _method;
    if (check_string(method_name, {"locally", "linear", "embedding"})) {
        _method=tapkee::KernelLocallyLinearEmbedding;
    }
    else if (check_string(method_name, {"isomap"})) {
        _method=tapkee::Isomap;
    }
    else if (check_string(method_name, {"multi", "dimensional", "scaling"})) {
        _method=tapkee::MultidimensionalScaling;
    }
    else if (check_string(method_name, {"kernel", "pca"})) {
        _method=tapkee::KernelPCA;
    }
    else if (check_string(method_name, {"laplacian", "eigenmaps"})) {
        _method=tapkee::LaplacianEigenmaps;
    }
    else if (check_string(method_name,
             {"linear", "local", "tangent", "space", "alignment"})) {
        _method=tapkee::LinearLocalTangentSpaceAlignment;
    }
    else {
        std::cout << "Unrecognized method name: " << method_name << '\n';
        return 1;
    }

    tapkee::NeighborsMethod _nmethod=tapkee::Brute;
    if (check_string(neighbors_method_name, {"cover", "tree"})) {
        _nmethod=tapkee::CoverTree;
    }
    else if (check_string(neighbors_method_name, {"vantage", "point", "tree"})) {
        _nmethod=tapkee::VpTree;
    }
    else if (!check_string(neighbors_method_name, {"brute"})) {
        std::cout << "Unrecognized neighbors method name: " << neighbors_method_name << '\n';
    }

    std::vector<vec3> points;
    std::vector<double> angles;
    swiss_roll(points, angles, npoints, regular_sampling);

    tapkee::TapkeeOutput output;

    distance_callback dcb(points);
    kernel_callback kcb(points);

    if (verbose>0) {
        std::cout << "Computing dimensionality reduction of "
            << npoints << " points using " << neighbors_method_name
            << " method to compute neighboring graph and "
            << method_name << " to compute the manifold coordinates\n";
    }

    try {
        nvis::timer _timer;
        output=tapkee::initialize()
            .withParameters((tapkee::method=_method,
                             tapkee::target_dimension=2,
                             tapkee::num_neighbors=n_neighbors,
                             tapkee::neighbors_method=_nmethod,
                             tapkee::check_connectivity=true))
            .withDistance(dcb)
            .withKernel(kcb)
            .embedRange(points.begin(), points.end());
        if (verbose>0) {
            std::cout << "dimensionality reduction took: " << _timer.elapsed()
                << " s.\n";
        }
    }
    catch(std::exception& e) {
        std::cerr << "exception caught: " << e.what() << '\n';
        exit(1);
    }

    std::vector<int> cc_ids;
    // check_connectivity(cc_ids, output.neighbors);

    Eigen::MatrixXd mins=output.embedding.colwise().minCoeff();
    Eigen::MatrixXd maxs=output.embedding.colwise().maxCoeff();

    if (verbose>0) {
        std::cout << "manifold coordinates bounding box: ["
                << mins(0,0) << ", " << maxs(0,0) << "] x ["
                << mins(0,1) << ", " << maxs(0,1) << "]\n";
    }

    VTK_SMART(PolyData) pd=vtk_utils::make_points(points);
    std::vector<double> values(npoints);
    for (int i=0; i<npoints; ++i) {
        values[i]=output.embedding(i,1);
    }
    pd=vtk_utils::add_scalars(pd, values);

    VTK_CREATE(vtkColorTransferFunction, ctf);
    {
        double min=mins(0,0);
        double max=maxs(0,0);
        ctf->AddRGBPoint(min, 0, 0, 1);
        ctf->AddRGBPoint(0.75*min+0.25*max, 0.5, 0.5, 1);
        ctf->AddRGBPoint(0.5*min+0.5*max, 1, 1, 1);
        ctf->AddRGBPoint(0.25*min+0.75*max, 1, 1, 0.5);
        ctf->AddRGBPoint(max, 1, 1, 0);
    }

    VTK_SMART(PolyData) spheres=vtk_utils::make_spheres(pd, 0.1);

    VTK_CREATE(vtkPolyDataMapper, point_mapper);
    point_mapper->SetInputData(spheres);
    point_mapper->ScalarVisibilityOn();
    point_mapper->SetLookupTable(ctf);
    VTK_CREATE(vtkActor, point_actor);
    point_actor->SetMapper(point_mapper);
    std::vector<int> edges;
    // for (size_t i=0; i<output.neighbors.size(); ++i) {
        const std::vector<tapkee::IndexType>& neighbors=output.neighbors[i];
        for (size_t j=0; j<neighbors.size(); ++j) {
            edges.push_back(i);
            edges.push_back(neighbors[j]);
        }
    }
    vtk_utils::add_edges_from_numbers(pd, edges);

    VTK_CREATE(vtkPolyDataMapper, edge_mapper);
    edge_mapper->SetInputData(pd);
    edge_mapper->ScalarVisibilityOn();
    edge_mapper->SetLookupTable(ctf);
    VTK_CREATE(vtkActor, edge_actor);
    edge_actor->SetMapper(edge_mapper);
    edge_actor->GetProperty()->EdgeVisibilityOn();

    // Initialize visualization pipeline and interactor
    VTK_CREATE(vtkRenderer, renderer);
    renderer->AddActor(point_actor);
    renderer->AddActor(edge_actor);

    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(800, 800);

    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    renderer->ResetCamera();
    interactor->Initialize();
    window->Render();

    std::vector<nvis::vec2> coords(npoints);
    for (int i=0; i<npoints; ++i) {
        coords[i][0]=output.embedding(i,0);
        coords[i][1]=output.embedding(i,1);
    }
    VTK_SMART(PolyData) pd2d=vtk_utils::make_points(coords);
    pd2d=vtk_utils::add_scalars(pd2d, values);
    vtk_utils::add_edges_from_numbers(pd2d, edges);
    VTK_CREATE(vtkPolyDataMapper, mapper2d);
    mapper2d->SetInputData(pd2d);
    mapper2d->ScalarVisibilityOn();
    mapper2d->SetLookupTable(ctf);
    VTK_CREATE(vtkActor, actor2d);
    actor2d->SetMapper(mapper2d);
    actor2d->GetProperty()->EdgeVisibilityOn();

    VTK_SMART(PolyData) spheres2d=vtk_utils::make_spheres(pd2d, 0.0001);
    VTK_CREATE(vtkPolyDataMapper, point2d_mapper);
    point2d_mapper->SetInputData(spheres2d);
    point2d_mapper->ScalarVisibilityOn();
    point2d_mapper->SetLookupTable(ctf);
    VTK_CREATE(vtkActor, point2d_actor);
    point2d_actor->SetMapper(point2d_mapper);

    VTK_CREATE(vtkRenderer, ren2d);
    ren2d->AddActor(actor2d);
    ren2d->AddActor(point2d_actor);
    VTK_CREATE(vtkRenderWindow, win2d);
    win2d->AddRenderer(ren2d);
    win2d->SetSize(800, 800);
    VTK_CREATE(vtkRenderWindowInteractor, interactor2d);
    interactor2d->SetRenderWindow(win2d);
    ren2d->ResetCamera();
    interactor2d->Initialize();
    win2d->Render();

    interactor->Start();


    return 0;
}
