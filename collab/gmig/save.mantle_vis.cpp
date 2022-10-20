#include <string>
#include <iostream>
#include <map>

#include "core_mantle.hpp"
#include "core_mantle_io.hpp"
#include "utils.hpp"

// spurt's utilities
#include <format/format.hpp>
#include <misc/option_parse.hpp>

// VTK
#include <vtkAxesActor.h>
#include <vtkCaptionActor2D.h>
#include <vtkColorTransferFunction.h>
#include <vtkDelaunay2D.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

// VTK helper functions and macros
#include <VTK/vtk_utils.hpp>
#include <graphics/colors.hpp>

// Boost graph library
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/utility.hpp>

// Qt
#include <QMenu>
#include <QFileDialog>
#include "QVTKInteractor.h"
#include "mantle_vis_control.hpp"
#include "mantle_vis_renderer.hpp"

// miscellaneous
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

std::string home_directory;

using namespace spurt::gmig;

typedef float  scalar_t;
typedef int    integer_t;

const scalar_t _core_radius_ =3486; // in km
const scalar_t _earth_radius_=6371; // in km

const scalar_t _infinity_=std::numeric_limits<scalar_t>::max();

const integer_t invalid_index=static_cast<integer_t>(-1);

typedef std::pair<integer_t, scalar_t>  edge_t;

typedef core_mantle::Vertex<scalar_t,integer_t>  vertex_t;
typedef vertex_t::pos_t    pos_t;
typedef vertex_t::bbox_t   bbox_t;
typedef core_mantle::VertexDistance<scalar_t, integer_t> distance_t;

// * a sector is comprised of sections
// * a section is a set of vertices contained on a given meridian
// * a layer is a set of points within a section that have same index and 
//   same type (+1/-1)
// * an interface is the surface formed by matching layers across sections
// * a front is the layer most recently added to an interface

typedef core_mantle::Section<scalar_t, integer_t>    section_t;
typedef std::vector<integer_t>                       layer_t;
typedef core_mantle::SectionID<scalar_t, integer_t>  section_id_t;
typedef std::vector<section_id_t>                    sector_t;
typedef core_mantle::LayerID<scalar_t, integer_t>    layer_id_t;
typedef core_mantle::SectorID<scalar_t, integer_t>   sector_id_t;
typedef std::vector<layer_id_t>                      interface_t;

typedef core_mantle::SetContainer<section_id_t, section_t, integer_t>  
    section_container;
typedef core_mantle::SetContainer<layer_id_t, layer_t, integer_t>      
    layer_container;
typedef core_mantle::SetContainer<sector_id_t, sector_t, section_id_t> 
    sector_container;

typedef section_container::iterator   section_iterator;
typedef layer_container::iterator     layer_iterator;
typedef sector_container::iterator    sector_iterator;

std::vector<vertex_t>    vertices_;
section_container        sections_;
layer_container          layers_;
sector_container         sectors_;
std::vector<interface_t> interfaces_;

typedef std::list<point_t> neighborhood_t;
std::vector<neighborhood_t> neighbors_;

typedef spurt::point_locator<scalar_t, integer_t, 3> locator_t;
typedef locator_t::point_type point_t;
typedef locator_t::coord_type coord_t;
typedef nvis::lexicographical_order order_t;

typedef LessDistFrom<coord_t, point_t> less_dist_from;
typedef LessPosTolerance<scalar_t, coord_t, order_t> less_pos_tol;

std::vector<integer_t> type_plus_one, type_minus_one;
std::vector< std::vector<integer_t> > connected_plus, connected_minus;

template <typename EdgeWeightMap>
struct neighborhood_filter {
  neighborhood_filter(EdgeWeightMap threshold, 
                      scalar_t m_threshold=params::max_search_radius)
          : m_weight(weight), m_threshold(threshold) {}
                      
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return get(m_weight, e) < m_threshold;
  }
  
  EdgeWeightMap m_weight;
  scalar_t m_threshold;
};

namespace boost {
typedef property<vertex_color_t, bool> vertex_prop_t;
typedef property<edge_weight_t, double> edge_prop_t;
typedef adjacency_list<vecS, vecS, undirectedS, vertex_prop_t, edge_prop_t> graph_t;
typedef graph_traits<graph_t> graph_traits_t;
typedef graph_traits_t::vertex_descriptor graph_vertex_t;
typedef graph_traits_t::edge_descriptor graph_edge_t;
typedef graph_t::vertex_iterator graph_vertex_iterator;
typedef graph_t::edge_iterator graph_edge_iterator;
template<typename Tag_>
using graph_map_t=property_map<graph_t, Tag_>;
typedef graph_map_t<edge_weight_t>::type edge_map_t;
typedef neighborhood_filter<edge_map_t> graph_filter_t;
typedef filtered_graph<graph_t, graph_filter_t> filtered_graph_t;
typedef graph_traits<filtered_graph_t> filtered_graph_traits_t;
template<typename Tag_>
using filter_graph_map_t=property_map<filtered_graph_t, Tag_>;
}

typedef Counter<integer_t> counter_t;
typedef nvis::fixed_vector<scalar_t, 3> color_t;

namespace params {
// Program parameters and options
std::string filename;       // input file name
int  verbose=0;             // output wordiness level
bool world_coords=false;    // map geocoordinates to spatial coordinates 
bool show_earth=false;      // display earth surface
bool show_axes=false;       // display coordinate axes
bool show_points=true;      // display vertices
bool show_edges=false;      // display edges
bool show_progress=false;   // show progress bar during computation
bool do_type_plus_one = true;   // process type +1 vertices
bool do_type_minus_one = false; // process type -1 vertices
scalar_t search_radius=3.;      // radius of neighbor search region
scalar_t max_search_radius=10.; // maximum size of search area
scalar_t sphere_radius=0.5;  // radius of sphere glyphs in point depiction
integer_t n_neighbors=3;     // number of neighbors used for connectivity
color_t bg_color=color_t(1,1,1); // background color
}

struct VertexDistance {
    scalar_t operator()(integer_t vid1, integer_t vid2) const {
        const pos_t& p1=vertices_[vid1].position;
        const pos_t& p2=vertices_[vid2].position;
        // ignore longitudinal distance
        pos_t subv1(p1[0], 0, p1[2]);
        pos_t subv2(p2[0], 0, p2[2]);
        return nvis::norm(subv1-subv2); 
    }
};

void ControlWindow::slot_points_check(bool checked) {
    show_points = checked;
    // update
}

void ControlWindow::slot_edges_checked(bool checked) {
    show_edges = checked;
    // update
}

void ControlWindow::slot_axes_checked(bool checked) {
    show_axes = checked;
    // update
}

void ControlWindow::slot_radius_slider(int r) {
    search_radius = r;
    // update
}

void ControlWindow::slot_radius_spinbox(double r) {
    search_radius = r;
    // update
}

void ControlWindow::slot_size_spinbox(int r) {
    n_neighbors = r;
    // update
}

void ControlWindow::camera_out() {
    QString fn = 
        QFileDialog::getSaveFileName(this, 
                                     QString("Select an output camera file"),
                                     QString(home_directory.c_str()));
    std::string filename = fn.toStdString();
    vtk_utils::export_camera_settings(filename, renderer);
}

void ControlWindow::camera_in() {
    QString fn = 
        QFileDialog::getOpenFileName(this,
                                     QString("Select an input camera file"),
                                     QString(home_directory.c_str()));
    std::string filename = fn.toStdString();
    vtk_utils::import_camera_settings(filename, renderer);
}

// int to double conversion slot (slider to spinbox signal)
void ControlWindow::slot_slider_to_spinbox(int value) {
    this->ui->radiusSpin->setValue(static_cast<double>(value)/1000.);
}

// double to int conversion slot (spinbox to slider signal)
void ControlWindow::slot_spinbox_to_slider(double value) {
    this->ui->radiusSlider->setValue(static_cast<int>(1000*value));
}

void ControlWindow::load_file(const std::string& data_name) {
    import_vertices(data_name);
    internal_main();
}

bbox_t import_vertices(const std::string& data_name) {
    // load vertices from file
    std::vector<vertex_t> __vertices;
    bbox_t bounds = 
        core_mantle::read_text(__vertices, params::filename, false);
    
    // filter out redundancy
    vertices_.clear();
    type_plus_one.clear();
    type_minus_one.clear();
    std::set<pos_t, less_pos_tol> unique_pos;
    integer_t n_points = 0;
    for (integer_t i=0; i<__vertices.size(); ++i) {
        const vertex_t v = __vertices[i];
        const pos_t& x=v.position;
        if (unique_pos.find(x) == unique_pos.end()) {
            unique_pos.insert(x);
            vertices_.push_back(v);
            if (v.type>0) type_plus_one.push_back(n_points++);
            else type_minus_one.push_back(n_points++);
        } // else skip redundant vertex
    }
    
    if (params::verbose) {
        std::cout 
            << __vertices.size() << " vertices imported from file\n"
            << vertices_.size() << " unique vertices\n"
            << type_plus_one.size() << " type +1\n"
            << type_minus_one.size() << " type -1\n"
            << "bounding box:\n" << bounds << '\n';
    }
    
    return bounds;
}

void compute_connectivity(boost::graph& graph, 
                          const std::vector<integer_t>& verts) {
    spurt::locator_t locator;
    for (integer_t i=0 ; i<verts.size() ; ++i) {
        integer_t id = verts[i];
        const vertex_t& v = vertices_[id];
        const pos_t& x = v.position;
        locator.insert(point_t(x, i));
    }
    
    Counter<size_t> neigh_sizes;
    progress.start(locator.size());
    size_t n=0;
    std::for_each(locator.begin(), locator.end(), 
                  [&] (const locator_t::point_type& p) {
        neighborhood_t neighbors;
        integer_t i=p.data();
        const pos_t& x=p.coordinate();
        locator.find_within_range(neighbors, x, params::max_search_radius);
        neighbors.sort(less_dist_from(x));
        neigh_sizes.increment(neighbors.size()-1);
        if (neighbors.empty()) {
            std::cout << "ERROR: no neighbor found at vertex #"
                     << i << "'s position\n";
        } 
        else if (neighbors.front().data() != i) {
            std::cout << "ERROR: closest point from vertex #"
                << i << "'s position is not vertex #"
                << i << " itself!" << '\n';
            std::cout << "position #" << i << "=" << x << '\n';
            std::cout << "closest point #" << neighbors.front().data()
                << "=" << neighbors.front().coordinate() << '\n';
            std::cout << "distance=" 
                << nvis::norm(neighbors.front().coordinate()-x)
                << '\n';
        }
        neighborhood_t::iterator it=neighbors.begin();
        for (++it; it!=neighbors.end() ; ++it) {
            scalar_t dist = nvis::norm(x-it->coordinate());
            if (dist > params::max_search_radius) break;
            boost::add_edge(i, it->data(), 
                            nvis::norm(x-it->coordinate()), 
                            graph);
        }
        progress.update(n++);
    });
    progress.end();
    
    if (params::verbose) {
        std::cout << boost::num_edges(graph) << " edges created overall\n";
        if (params::verbose > 1) {
            std::cout << "Neighborhood size distribution\n";
            std::cout << neigh_sizes << '\n';
        }
    }
}

filtered_graph_t 
filter_connectivity(const boost::graph& in) {
    graph_filter_t filter(boost::get(edge_weight, in), params::search_radius);
    return filtered_graph_t(in, filter);
}

vtkSmartPointer<vtkColorTransferFunction>
create_color_transfer_function(const std::vector<scalar_t>& values) {
    vtkSmartPointer<vtkColorTransferFunction> ctf = 
        vtkSmartPointer<vtkColorTransferFunction>::New();
    
    // adaptive color map with spiral color scale
    std::vector<color_t> colors(20);
    spurt::spiral_scale(colors, 20, 0.2);
    scalar_t minval=*std::min_element(values.begin(), values.end());
    scalar_t maxval=*std::max_element(values.begin(), values.end());
    scalar_t dval=(maxval-minval)/19;
    
    for (int i=0; i<20; ++i) {
        scalar_t value=minval + i*dval;
        color_t color=colors[i];
        ctf->AddRGBPoint(value, color[0], color[1], color[2]);
    }
    return ctf;
}

vtkSmartPointer<vtkPolyData> 
update_geometry(vtkSmartPointer<vtkPolyData>& pd,
                const std::vector<integer_t>& indices) {
    
    typedef vtk_utils::vtk_array_traits<scalar_t>::array_type array_type;

    integer_t n_points = indices.size();
    
    VTK_CREATE(array_type, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(n_points);
    size_t counter=0;
    for (int i=0 ; i<n_points ; ++i, ++counter) {
        pos_t x = vertices_[indices[i]].coordinate();
        if (!params::world_coords) {
            // rescaling to improve aspect ratio
            x[1] *= 2;
            x[2] /= 3;
        }
        else {
            scalar_t r=_core_radius_ + x[2];
            scalar_t la=deg2rad(x[0]);
            scalar_t lo=deg2rad(x[1]);
            x = pos_t(r*sin(la)*cos(lo),
                      r*sin(la)*sin(lo),
                      r*cos(la));
            coords->SetTuple(counter, &x[0]);
        }
    }

    VTK_CREATE(vtkPoints, points);
    points->SetData(coords);
    pd->SetPoints(points); // old points, if any, will be discarded
    return pd; // return updated polydata object
}

bool init(int argc, char** argv) {
    namespace xcl=spurt::command_line;
    
    xcl::option_traits 
        required_group("Required parameters", true, false),
        positional_group("Positional parameters", true, true),
        optional_group("Optional parameters", false, false);
        
    xcl::option_parser parser(argv[0], 
        "Extract and visualize core-mantle interfaces");
    
    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", params::filename, "Input file (string)", 
                         positional_group);
        parser.add_value("radius", params::search_radius, 
                         "Radius of neighbor search region",
                         optional_group);
        parser.add_value("neighbors", params::n_neighbors, 
                         "Number of neighbors used to determine connectivity",
                         optional_group);
        parser.add_value("points", params::show_points,
                        "Show vertices", optional_group);
        parser.add_value("edges", params::show_edges,
                        "Show edges", optional_group);
        parser.add_value("positive", params::do_type_plus_one,
                        "Process type +1 vertices", optional_group);
        parser.add_value("negative", params::do_type_minus_one,
                        "Process type -1 vertices", optional_group);
        parser.add_value("sphere", params::sphere_radius,
                        "Radius of sphere glyphs", optional_group);
        parser.add_value("verbose", params::verbose,
                         "verbose level", optional_group);
        parser.add_tuple<3>("bg", params::bg_color, params::bg_color, 
                            "Background color", optional_group);
        parser.add_flag("world", params::world_coords, 
                        "Display data in spherical coordinates", 
                        optional_group);
        parser.add_flag("earth", params::show_earth, "Show earth crust", 
                        optional_group);
        parser.add_flag("axes", params::show_axes, "Show coordinates axes", 
                        optional_group);
        parser.add_flag("progress", params::show_progress,
                        "Show progress bar", optional_group);
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

vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> window;
vtkSmartPointer<vtkRenderWindowInterator> interactor;
vtkSmartPointer<vtkPolyData> polydata_plus, polydata_minus;
vtkSmartPointer<vtkOrientationMarkerWidget> axes_widget;
vtkSmartPointer<vtkColorTransferFunction> color_tf;
vtkSmartPointer<vtkActor> points_actor, edges_actor, faces_actor; 

vtkSmartPointer<vtkOrientationMarkerWidget> coordinate_axes_widget() {
    color_t text_color = 
        spurt::luminosity(params::bg_color) > 0.5 ?
        spurt::black : spurt::white;
    color_t text_bg_color = 0.5*(text_color + params::bg_color);
    vtkSmartPointer<vtkAxesActor> axes=
        vtkSmartPointer<vtkAxesActor>::New();
    axes->SetCylinderRadius(0.05);
    axes->SetShaftTypeToCylinder();             
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(text_color[0], text_color[1], text_color[2]);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(text_color[0], text_color[1], text_color[2]);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(text_color[0], text_color[1], text_color[2]);
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(text_bg_color[0], text_bg_color[1], 
                           text_bg_color[2]);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(text_bg_color[0], text_bg_color[1], 
                           text_bg_color[2]);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(text_bg_color[0], text_bg_color[1], 
                           text_bg_color[2]);
    axes->GetXAxisCaptionActor2D()->SetCaption("latitude");
    axes->GetYAxisCaptionActor2D()->SetCaption("longitude");
    axes->GetZAxisCaptionActor2D()->SetCaption("height");
 
    vtkSmartPointer<vtkOrientationMarkerWidget> widget=
      vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor(bg_color[0], bg_color[1], bg_color[2]);
    widget->SetOrientationMarker(axes);
    widget->SetInteractor( interactor );
    widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    return widget;
}

void connected_components(std::vector< std::vector<integer_t> >& comp,
                          const boost::filtered_graph_t& graph) {
    // compute connected components in graph
    std::vector<integer_t> components(graph.size());
    integer_t num=boost::connected_components(graph, &components[0]);
    if (params::verbose) {
        std::cout << "There are " << num << " connected components\n";
    }
    comp.clear();
    comp.resize(num);
    for (int i=0 ; i<components.size() ; ++i) {
        comp[components[i]].push_back(i);
    }
}

vtkSmartPointer<vtkPolyData> 
assign_values(vtkSmartPointer<vtkPolyData> pd, 
              const std::vector< std::vector<integer_t> >& comp,
              bool update_ctf=false) {
    std::vector<scalar_t> values(pd->GetPoints()->GetNumberOfPoints());
    for (int i=0; i<comp.size() ; ++i) {
        std::for_each(com[i].begin(), comp[i].end(), 
                      [&](integer_t j) { values[j] = i; });
    }
    // Update color mapping
    if (update_ctf) {
        color_tf = create_color_transfer_function(values);
    }
    if (params::verbose) {
        std::cout << values.size() << " values added\n";
    }
    return vtk_utils::add_scalars(pd, values);
}

vtkSmartPointer<vtkActor> display_points(vtkSmartPointer<vtkPolyData> pd) {
    vtkSmartPointer<vtkPolyData> spheres=
        vtk_utils::make_spheres(pd, params::sphere_radius);
    VTK_MAKE_ACTOR(sphere_actor, spheres);
    sphere_actor->GetMapper()->SetLookupTable(color_transfer_function);
    sphere_actor->GetMapper()->ScalarVisibilityOn();
    return sphere_actor;
}

vtkSmartPointer<vtkActor> 
display_edges(vtkSmartPointer<vtkPolyData> pd,
              const boost::filtered_graph_t& graph) {
    std::vector<integer_t> edge_ids;
    typedef boost::filtered_graph_traits_t::edge_iterator edge_iterator_t;
    edge_iterator_t ei, ei_end;
    for (boost::tie(ei, ei_end)=boost::edges(graph); ei != ei_end; ++ei) {
        edge_ids.push_back(boost::source(*ei, graph));
        edge_ids.push_back(boost::target(*ei, graph));
    }
    if (params::verbose) {
        std::cout << edge_ids.size()/2 << " edges created\n";
    }
    vtk_utils::add_edges_from_numbers(pd, edge_ids);
    VTK_MAKE_ACTOR(edge_actor, pd);
    edge_actor->GetMapper()->SetLookupTable(color_transfer_function);
    edge_actor->GetMapper()->ScalarVisibilityOn();
    edge_actor->GetProperty()->SetLineWidth(2);
    return edge_actor;                   
}

int internal_main(const std::string& filename) {
    import_vertices(filename);
    
    srand48(time(NULL));
    ProgressDisplay progress(params::show_progress);
    
    // determine superset of connectivity
    boost::graph_t graph_plus, graph_minus;
    edge_map_t edge_weight_plus, edge_weight_minus;
    edge_weight_plus=get(boost::edge_weight, graph_plus);
    edge_weight_minus=get(boost::edge_weight, graph_minus);
    compute_connectivity(graph_plus, type_plus_one);
    compute_connectivity(graph_minus type_minus_one);
    
    // identify subset of connectivity corresponding to chosen radius
    filtered_graph_t fg_plus = filter_connectivity(graph_plus);
    filtered_graph_t fg_minus = filter_connectivity(graph_minus);
    
    if (params::verbose) {
        std::cout << boost::num_edges(fg_plus) 
                  << " edges created for type +1\n";
        std::cout << boost::num_edges(fg_minus)
                  << " edges created for type -1\n";
    }
    
#if 0
    std::vector<scalar_t> 
        degree_plus(type_plus_one.size()),
        degree_minus(type_minus_one.size());
    std::pair<filtered_graph_t::vertex_iterator, 
              filtered_graph_t::vertex_iterator> iter_range;
    filtered_map_t<vertex_index>::type vertex_id_map =
        boost::get(vertex_index, fg_plus);
    for (iter_range=boost::vertices(fg_plus); 
         iter_range.first != iter_range.end(); ++iter_range.first) {
        auto v = *iter_range.first;
        integer_t i = vertex_id_map[v];
        degree_plus[i] = boost::out_degree(i, fg_plus);
    }
    vertex_id_map=boost::get(vertex_index, fg_minus);
    for (iter_range=boost:vertices(fg_minus); 
         iter_range.first!=iter_range.second; ++iter.first) {
        auto v = *iter_range.first;
        integer_t i = vertex_id_map[v];
        degree_minus[i] = boost::out_degree(i, fg_minus);
    }
#endif
    
    // Create polydata to be visualized
    polydata_plus = update_geometry(polydata_plus, type_plus_one);
    polydata_minus = update_geometry(polydata_minus, type_minus_one);
    polydata_plus = vtk_utils::make_points(positions_plus);
    polydata_minus = vtk_utils::make_points(positions_minus);
    
    // Initialize visualization pipeline and interactor
    renderer = vtkSmartPointer<vtkRenderer>::New();
    window = vtkSmartPointer<vtkRenderWindow>::New();
    window->AddRenderer(renderer);
    window->SetSize(800, 800);
    interactor = vtkRenderWindowInteractor::New();
    interactor->SetRenderWindow(window);
    
    // coordinate axes
    axes_widget = coordinate_axes_widget();
    
    // initialize connected components
    connected_components(connected_plus, fg_plus);
    connected_components(connected_minus, fg_minus);
    
    if (params::verbose) {
        std::cout << "There are " << connected_plus.size() 
                  << " connected components for type +1 and "
                  << connected_minus.size()
                  << " connected components for type -1\n";
    }
    
    // assign connected component id to each vertex
    polydata_plus = assign_values(polydata_plus, connected_plus, true);
    polydata_minus = assign_values(polydata_minus, connected_minus);
    
    if (params::show_edges) {
        std::vector<integer_t> edge_ids;
        typedef boost::graph_traits_t::edge_iterator edge_iterator_t;
        edge_iterator_t ei, ei_end;
        for (boost::tie(ei, ei_end)=boost::edges(graph); ei != ei_end; ++ei) {
            edge_ids.push_back(boost::source(*ei, graph));
            edge_ids.push_back(boost::target(*ei, graph));
        }
        
        if (params::verbose) {
            std::cout << edge_ids.size()/2 << " edges created\n";
        }
        vtk_utils::add_edges_from_numbers(pd, edge_ids);
        VTK_MAKE_ACTOR(edge_actor, pd);
        edge_actor->GetMapper()->SetLookupTable(ctf);
        edge_actor->GetMapper()->ScalarVisibilityOn();
        edge_actor->GetProperty()->SetLineWidth(2);
        renderer->AddActor(edge_actor);
        
#if 0
        VTK_CREATE(vtkTubeFilter, tubef);
        VTK_CONNECT(tubef, pd);
        tubef->SetRadius(1);
        tubef->SetNumberOfSides(8);
        tubef->Update();
        VTK_MAKE_ACTOR(tube_actor, tubef->GetOutput());
        
        // tube_actor->GetProperty()->SetColor(1,0,0);
        tube_actor->GetMapper()->SetLookupTable(ctf);
        tube_actor->GetMapper()->ScalarVisibilityOn();
        renderer->AddActor(tube_actor);
#endif
    }
    
    if (params::show_points) {
        vtkSmartPointer<vtkPolyData> spheres=
            vtk_utils::make_spheres(pd, params::sphere_radius);
        VTK_MAKE_ACTOR(sphere_actor, spheres);
        sphere_actor->GetMapper()->SetLookupTable(ctf);
        sphere_actor->GetMapper()->ScalarVisibilityOn();
        renderer->AddActor(sphere_actor);
    }
    
    if (params::show_earth) {
        std::cout << "showing earth crust\n";
        VTK_CREATE(vtkSphereSource, earth);
        earth->SetRadius(_earth_radius_);
        earth->SetCenter(0,0,0);
        earth->SetThetaResolution(50);
        earth->SetPhiResolution(50);
        earth->Update();
        VTK_MAKE_ACTOR(earth_actor, earth->GetOutput());
        earth_actor->GetProperty()->SetColor(0,0,0);
        earth_actor->GetProperty()->SetOpacity(0.1);
        renderer->AddActor(earth_actor);
    }
    
    // setup interaction loop
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    
    // configure rendering window
    renderer->SetBackground(params::bg_color[0],
                            params::bg_color[1],
                            params::bg_color[2]);
    color_t text_color = 
        spurt::luminosity(params::bg_color) > 0.5 ?
        spurt::black : spurt::white;
    color_t text_bg_color = 0.5*(text_color + params::bg_color);
    
    vtkSmartPointer<vtkAxesActor> axes=
        vtkSmartPointer<vtkAxesActor>::New();
    
    double __cyl_radius=axes->GetCylinderRadius();
    double __cone_radius=axes->GetConeRadius();
    int __cyl_res=axes->GetCylinderResolution();
    std::cout << "cylinder radius=" << __cyl_radius << '\n';
    std::cout << "cone radius=" << __cone_radius << '\n';
    std::cout << "cylinder resolution=" << __cyl_res << '\n';
    
    axes->SetCylinderRadius(0.05);
    axes->SetShaftTypeToCylinder();             
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(text_color[0], text_color[1], text_color[2]);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(text_color[0], text_color[1], text_color[2]);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(text_color[0], text_color[1], text_color[2]);
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(text_bg_color[0], text_bg_color[1], 
                           text_bg_color[2]);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(text_bg_color[0], text_bg_color[1], 
                           text_bg_color[2]);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(text_bg_color[0], text_bg_color[1], 
                           text_bg_color[2]);
    axes->GetXAxisCaptionActor2D()->SetCaption("latitude");
    axes->GetYAxisCaptionActor2D()->SetCaption("longitude");
    axes->GetZAxisCaptionActor2D()->SetCaption("height");
 
    vtkSmartPointer<vtkOrientationMarkerWidget> widget=
      vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor( bg_color[0], bg_color[1], bg_color[2] );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( interactor );
    widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    if (params::show_axes) {
        widget->SetEnabled( 1 );
        widget->InteractiveOn();
    }

    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    window->Render();
    
    interactor->Initialize();
    interactor->Start();
       
    return 0; 
}

int main(int argc, char* argv[])
{
    if (!init(argc, argv)) {
        std::cout << "Invalid command line options. Exiting.\n";
        exit(1);
    }
    
    bbox_t bounds=core_mantle::read_text(vertices_, params::filename, false);
    std::cout << vertices_.size() << " vertices imported across region:\n"
              << bounds << '\n';
    
    srand48(time(NULL));
    ProgressDisplay progress(params::show_progress);
        
    // split set of vertices into type categories
    std::vector<integer_t> type_plus, type_minus;
    for (integer_t i=0; i<vertices_.size(); ++i) {
        if (vertices_[i].type > 0) type_plus.push_back(i);
        else type_minus.push_back(i); 
    }
    
    if (params::verbose) {
        std::cout << "There are " << type_plus.size() 
            << " type +1 vertices and " << type_minus.size()
            << " type -1 vertices\n";
    }
    
    std::map<int, int> old2new, new2old;
    locator_t locator;
    std::set<pos_t, less_pos_tol> unique_pos;
    integer_t n_points=0;
    std::for_each(type_plus.begin(), type_plus.end(), [&] (size_t i) {
        const pos_t& x=vertices_[i].position;
        if (unique_pos.find(x) == unique_pos.end()) {
            unique_pos.insert(x);
            old2new[i] = n_points;
            new2old[n_points] = i;
            locator.insert(point_t(x, n_points++));
        } // else skip redundant vertex
    });
    
    // create graph data structure 
    boost::graph_t graph(unique_pos.size());
    boost::property_map<boost::graph_t, boost::edge_weight_t>::type eweight;
    eweight=get(boost::edge_weight, graph);
    
    if (params::verbose) {
        std::cout << unique_pos.size() << " (" << n_points 
                  << ") points added to kd-tree\n";
    }
    
    Counter<size_t> neigh_sizes;
    progress.start(locator.size());
    size_t n=0;
    std::for_each(locator.t.begin(), locator.t.end(), 
                  [&] (const locator_t::point_type& p) {
        std::list<point_t> neighbors;
        integer_t i=p.data();
        const pos_t& x=p.coordinate();
        locator.find_within_range(neighbors, x, params::search_radius);
        neighbors.sort(less_dist_from(x));
        neigh_sizes.increment(neighbors.size()-1);
        if (neighbors.empty()) {
            std::cout << "ERROR: no neighbor found at vertex #"
                     << i << "'s position\n";
        } 
        else if (neighbors.front().data() != i) {
            std::cout << "ERROR: closest point from vertex #"
                << i << "'s position is not vertex #"
                << i << " itself!" << '\n';
            std::cout << "position #" << i << "=" << x << '\n';
            std::cout << "closest point #" << neighbors.front().data()
                << "=" << neighbors.front().coordinate() << '\n';
            std::cout << "distance=" 
                << nvis::norm(neighbors.front().coordinate()-x)
                << '\n';
        }
        std::list<point_t>::iterator it=neighbors.begin();
        ++it;
        for (int k=0; it!=neighbors.end() && k<params::n_neighbors;
             ++it, ++k) {
            scalar_t dist = nvis::norm(x-it->coordinate());
            if (dist > params::search_radius) break;
            boost::add_edge(i, it->data(), 
                            nvis::norm(x-it->coordinate()), 
                            graph);
        }
        progress.update(n++);
    });
    progress.end();
    
    if (params::verbose) {
        std::cout << boost::num_edges(graph) << " edges created\n";
        // std::cout << "Neighborhood size statistics (for radius="
        //           << params::search_radius << "):\n"
        //           << neigh_sizes << '\n';
    }
       
    std::vector<pos_t> positions(locator.size());
    std::for_each(locator.begin(), locator.end(), 
                  [&](const locator_t::point_type& p)
    {
        if (!params::world_coords) {
            positions[p.data()] = p.coordinate();
            // rescaling to improve aspect ratio
            positions[p.data()][1] *= 2;
            positions[p.data()][2] /= 3;
        }
        else {
            scalar_t r=_core_radius_ + p.coordinate()[2];
            scalar_t la=deg2rad(p.coordinate()[0]);
            scalar_t lo=deg2rad(p.coordinate()[1]);
            pos_t q(r*sin(la)*cos(lo),
                    r*sin(la)*sin(lo),
                    r*cos(la));
            positions[p.data()] = q;
        }
    });

    std::vector<scalar_t> values(locator.size());
    
    std::cout << "Summary:\n"
              << '\t' << vertices_.size() << " vertices\n"
              << '\t' << type_plus.size() << " of type +1\n"
              << '\t' << type_minus.size() << " of type -1\n"
              << '\t' << positions.size() << " positions\n";
    
    // compute connected components in graph
    std::vector<integer_t> components(locator.size());
    integer_t num=boost::connected_components(graph, &components[0]);
    if (params::verbose) {
        std::cout << "There are " << num << " connected components\n";
    }
    
    Counter<integer_t> comp_sizes;
    std::for_each(components.begin(), components.end(), [&](integer_t c) {
        comp_sizes.increment(c);
    });
    
    for (auto it=locator.begin(); it!=locator.end(); ++it) {
        values[it->data()] = comp_sizes[components[it->data()]];
    }
        
    // create rendering engine and render window
    VTK_CREATE(vtkRenderer, renderer);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(800, 800);
    
    vtkSmartPointer<vtkPolyData> pd=vtk_utils::make_points(positions);
    if (params::verbose) {
        pd->ComputeBounds();
        double _bounds[6];
        pd->GetBounds(_bounds);
        std::cout << "bounds of created polydata are: \n"
                  << _bounds[0] << " -> " << _bounds[1] << "\n"
                  << _bounds[2] << " -> " << _bounds[3] << "\n"
                  << _bounds[4] << " -> " << _bounds[5] << '\n';
    }
    vtk_utils::add_scalars(pd, values);
    
    if (params::verbose) {
        std::cout << values.size() << " values added\n";
    }
    
    // create color map
    VTK_CREATE(vtkColorTransferFunction, ctf);
    {
        // adaptive color map with spiral color scale
        std::vector<color_t> colors(20);
        spurt::spiral_scale(colors, 20, 0.2);
        scalar_t minval=*std::min_element(values.begin(), values.end());
        scalar_t maxval=*std::max_element(values.begin(), values.end());
        scalar_t dval=(maxval-minval)/19;
        
        for (int i=0; i<20; ++i) {
            scalar_t value=minval + i*dval;
            color_t color=colors[i];
            ctf->AddRGBPoint(value, color[0], color[1], color[2]);
        }
    }
    
    if (params::show_edges) {
        std::vector<integer_t> edge_ids;
        typedef boost::graph_traits_t::edge_iterator edge_iterator_t;
        edge_iterator_t ei, ei_end;
        for (boost::tie(ei, ei_end)=boost::edges(graph); ei != ei_end; ++ei) {
            edge_ids.push_back(boost::source(*ei, graph));
            edge_ids.push_back(boost::target(*ei, graph));
        }
        
        if (params::verbose) {
            std::cout << edge_ids.size()/2 << " edges created\n";
        }
        vtk_utils::add_edges_from_numbers(pd, edge_ids);
        VTK_MAKE_ACTOR(edge_actor, pd);
        edge_actor->GetMapper()->SetLookupTable(ctf);
        edge_actor->GetMapper()->ScalarVisibilityOn();
        edge_actor->GetProperty()->SetLineWidth(2);
        renderer->AddActor(edge_actor);
        
#if 0
        VTK_CREATE(vtkTubeFilter, tubef);
        VTK_CONNECT(tubef, pd);
        tubef->SetRadius(1);
        tubef->SetNumberOfSides(8);
        tubef->Update();
        VTK_MAKE_ACTOR(tube_actor, tubef->GetOutput());
        
        // tube_actor->GetProperty()->SetColor(1,0,0);
        tube_actor->GetMapper()->SetLookupTable(ctf);
        tube_actor->GetMapper()->ScalarVisibilityOn();
        renderer->AddActor(tube_actor);
#endif
    }
    
    if (params::show_points) {
        vtkSmartPointer<vtkPolyData> spheres=
            vtk_utils::make_spheres(pd, params::sphere_radius);
        VTK_MAKE_ACTOR(sphere_actor, spheres);
        sphere_actor->GetMapper()->SetLookupTable(ctf);
        sphere_actor->GetMapper()->ScalarVisibilityOn();
        renderer->AddActor(sphere_actor);
    }
    
    if (params::show_earth) {
        std::cout << "showing earth crust\n";
        VTK_CREATE(vtkSphereSource, earth);
        earth->SetRadius(_earth_radius_);
        earth->SetCenter(0,0,0);
        earth->SetThetaResolution(50);
        earth->SetPhiResolution(50);
        earth->Update();
        VTK_MAKE_ACTOR(earth_actor, earth->GetOutput());
        earth_actor->GetProperty()->SetColor(0,0,0);
        earth_actor->GetProperty()->SetOpacity(0.1);
        renderer->AddActor(earth_actor);
    }
    
    // setup interaction loop
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    
    // configure rendering window
    renderer->SetBackground(params::bg_color[0],
                            params::bg_color[1],
                            params::bg_color[2]);
    
    vtkSmartPointer<vtkAxesActor> axes=
        vtkSmartPointer<vtkAxesActor>::New();
    
    double __cyl_radius=axes->GetCylinderRadius();
    double __cone_radius=axes->GetConeRadius();
    int __cyl_res=axes->GetCylinderResolution();
    std::cout << "cylinder radius=" << __cyl_radius << '\n';
    std::cout << "cone radius=" << __cone_radius << '\n';
    std::cout << "cylinder resolution=" << __cyl_res << '\n';
    
    axes->SetCylinderRadius(0.05);
    axes->SetShaftTypeToCylinder();
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0.5,0.5,0.5);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0.5,0.5,0.5);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0.5,0.5,0.5);
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetBackgroundColor(0.75,0.75,0.75);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetBackgroundColor(0.75,0.75,0.75);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetBackgroundColor(0.75,0.75,0.75);
 
    vtkSmartPointer<vtkOrientationMarkerWidget> widget=
      vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor( 1, 1, 1 );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( interactor );
    widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    if (params::show_axes) {
        widget->SetEnabled( 1 );
        widget->InteractiveOn();
    }

    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    window->Render();
    
    interactor->Initialize();
    interactor->Start();
       
    return 0;
}