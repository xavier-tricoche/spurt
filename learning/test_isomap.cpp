#include <vector>

#include <misc/option_parse.hpp>
#include <learning/isomap.hpp>
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

#include <math/types.hpp>

using namespace spurt;

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

bool show_points=true;
bool show_edges=true;
bool regular_sampling=false;
int npoints=100;
int n_neighbors=8;
int verbose=0;
int color_dim=0;
double search_radius=1;
double sphere_radius=0.1;
spurt::vec3 bg_color;

bool init(int argc, const char* argv[]) {
    namespace xcl=spurt::command_line;
    
    xcl::option_traits 
        required_group(true, false, "Required parameters"),
        positional_group(true, true, "Positional parameters"),
        optional_group(false, false, "Optional parameters");
        
    xcl::option_parser parser(argv[0], 
        "Test isomap implementation using classical swiss roll");
    
    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("n", npoints, "Number of points", 
                         optional_group);
        parser.add_value("regular", regular_sampling, 
                         "Regular manifold sampling", optional_group);
        parser.add_value("radius", search_radius, 
                         "Radius of neighbor search region",
                         optional_group);
        parser.add_value("neighbors", n_neighbors, 
                         "Number of neighbors used to determine connectivity",
                         optional_group);
        parser.add_value("points", show_points,
                        "Show vertices", optional_group);
        parser.add_value("edges", show_edges,
                        "Show edges", optional_group);
        parser.add_value("sphere", sphere_radius,
                        "Radius of sphere glyphs", optional_group);
        parser.add_value("verbose", verbose,
                         "verbose level", optional_group);
        parser.add_value("coldim", color_dim,
                         "Dimension used for manifold color coding", 
                         optional_group);
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
    
    std::vector<vec3> points;
    std::vector<double> angles;
    swiss_roll(points, angles, npoints, regular_sampling);
    
    typedef boost::property<boost::vertex_color_t, bool> vertex_prop_t;
    typedef boost::property<boost::edge_weight_t, double> edge_prop_t;
    typedef spurt::graph<vertex_prop_t, 
                          boost::property<boost::edge_weight_t, double, 
                          boost::property<boost::edge_weight2_t, double> > > graph_t;
    typedef spurt::Isomap<double, graph_t> isomap_t;
    typedef isomap_t::vector_t vector_t;
    typedef isomap_t::row_vector_t row_vector_t;
    
    graph_t graph;
    spurt::compute_knn_connectivity<graph_t, vec3, 3>(graph, n_neighbors, points);
    
    isomap_t::matrix_t coordinates;
    bool failed = false;
    try {
        spurt::timer t;
        isomap_t isomap(graph);
        isomap.embed();
        if (verbose>0) {
            std::cout << "total computation time for isomap=" 
                << t.elapsed() << " seconds.\n";
        }
    
        vector_t evals=isomap.eigenvalues(npoints);
#ifdef _SHOW_EMBEDDING_ERROR_
        for (int i=0 ; i<10 ; ++i) {
            std::cout << "eigenvalue #" << i+1 << ": " << evals(npoints-i-1) << '\n';
        }
        
        for (int dim=1; dim<10; ++dim) {
            std::cout << "Error norm for dim=" 
                      << dim << ": " << isomap.l2_error(dim) << '\n';
        }
#endif
        coordinates=isomap.coordinates(2); // (Nxdim)
    }
    catch (std::exception& e) {
        std::cout << "Exception caught: " << e.what() << '\n';
        std::cout << "Exiting.\n";
        failed = true;
    }

    row_vector_t mins, maxs, spans;
    if (!failed) {
        mins = coordinates.colwise().minCoeff(); // (1xdim)
        maxs = coordinates.colwise().maxCoeff(); // (1xdim)
        spans = maxs-mins; // (1xdim)
    }
    VTK_SMART(vtkPolyData) pd = vtk_utils::make_points(points);
    
    /*
    typedef spurt::fixed_vector<unsigned char, 3> color_t;
    const vec3 red = vec3(255,0,0); // (1,0)
    const vec3 blue = vec3(0,0,255); // (0,1)
    const vec3 white = vec3(255,255,255); // (1,1)
    const vec3 black = vec3(0,0,0); // (0,0)
    std::vector<color_t> colors(npoints);
    for (int i=0; i<npoints; ++i) {
        row_vector_t loc = (coordinates.row(i)-mins).array()/spans.array();
        double u=loc(0);
        double v=loc(1);
        if (u<0 || u>1 || v<0 || v>1) {
            std::cout << "ERROR: wrong local coordinates: " << loc << '\n';
        }
        colors[i] = color_t((1.-u)*(1.-v)*black + 
                            u*(1.-v)*red + 
                            u*v*white + 
                            (1.-u)*v*blue);
    }
    pd = vtk_utils::add_colors(pd, colors);
    */
    VTK_CREATE(vtkColorTransferFunction, ctf);
    std::vector<double> values(npoints);
    if (!failed) {
        for (int i=0; i<npoints; ++i) {
            values[i] = coordinates(i, color_dim);
        }
        pd = vtk_utils::add_scalars(pd, values);
        double min=mins(0, color_dim);
        double max=maxs(0, color_dim);
        ctf->AddRGBPoint(min, 0, 0, 1);
        ctf->AddRGBPoint(0.75*min+0.25*max, 0.5, 0.5, 1);
        ctf->AddRGBPoint(0.5*min+0.5*max, 1, 1, 1);
        ctf->AddRGBPoint(0.25*min+0.75*max, 1, 1, 0.5);
        ctf->AddRGBPoint(max, 1, 1, 0);
    }
    
    VTK_SMART(vtkPolyData) spheres = vtk_utils::make_spheres(pd, 0.1);
    
    VTK_CREATE(vtkPolyDataMapper, point_mapper);
    point_mapper->SetInputData(spheres);
    if (!failed) {
        point_mapper->ScalarVisibilityOn();
        point_mapper->SetLookupTable(ctf);
    }
    VTK_CREATE(vtkActor, point_actor);
    point_actor->SetMapper(point_mapper);
    
    typedef boost::graph_traits<graph_t>::edge_iterator edge_iter_t;
    edge_iter_t eit, eend;
    std::vector<int> edges;
    for (boost::tie(eit, eend)=boost::edges(graph); eit!=eend; ++eit) {
        int i=boost::source(*eit, graph);
        int j=boost::target(*eit, graph);
        edges.push_back(i);
        edges.push_back(j);
    }
    vtk_utils::add_edges_from_numbers(pd, edges);
    
    VTK_CREATE(vtkPolyDataMapper, edge_mapper);
    edge_mapper->SetInputData(pd);
    edge_mapper->ScalarVisibilityOn();
    if (!failed)
        edge_mapper->SetLookupTable(ctf);
    VTK_CREATE(vtkActor, edge_actor);
    edge_actor->SetMapper(edge_mapper);
    
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
    
    // visualization pipeline for manifold coordinates
    if (!failed) {
        std::vector<spurt::vec2> coords(npoints);
        for (int i=0; i<npoints; ++i) {
            coords[i][0]=coordinates(i, 0);
            coords[i][1]=coordinates(i, 1);
        }
        VTK_SMART(vtkPolyData) pd2d=vtk_utils::make_points(coords);
        pd2d=vtk_utils::add_scalars(pd2d, values);
        vtk_utils::add_edges_from_numbers(pd2d, edges);
        VTK_CREATE(vtkPolyDataMapper, mapper2d);
        mapper2d->SetInputData(pd2d);
        mapper2d->ScalarVisibilityOn();
        mapper2d->SetLookupTable(ctf);
        VTK_CREATE(vtkActor, actor2d);
        actor2d->SetMapper(mapper2d);
        actor2d->GetProperty()->EdgeVisibilityOn();
    
        VTK_SMART(vtkPolyData) spheres2d=vtk_utils::make_spheres(pd2d, 0.0001);
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
    }
    
    interactor->Start();
    
    return 0;
}