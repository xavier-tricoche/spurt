#include <iostream>
#include <map>
#include <memory>
#include <string>


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
#include <boost/graph/filtered_graph.hpp>
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

// **************************************************************************
// 
//                             Type Definitions
//
// **************************************************************************
typedef float  scalar_t;
typedef int    integer_t;

typedef std::pair<integer_t, scalar_t>  edge_t;

typedef core_mantle::Vertex<scalar_t, integer_t>  vertex_t;
typedef vertex_t::pos_t    pos_t;
typedef vertex_t::bbox_t   bbox_t;

typedef spurt::point_locator<scalar_t, integer_t, 3> locator_t;
typedef locator_t::point_type point_t;
typedef locator_t::coord_type coord_t;
typedef nvis::lexicographical_order order_t;

typedef std::list<point_t> neighborhood_t;

typedef LessDistFrom<coord_t, point_t> less_dist_from;
typedef LessPosTolerance<scalar_t, coord_t, order_t> less_pos_tol;
typedef Counter<integer_t> counter_t;
typedef nvis::fixed_vector<scalar_t, 3> color_t;

template <typename EdgeWeightMap>
struct neighborhood_filter {
    neighborhood_filter(scalar_t threshold=0) {}
    neighborhood_filter(EdgeWeightMap weight, scalar_t threshold=0)
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

struct spurt::gmig::TypeAttributes {
    TypeAttributes() : graph(new boost::graph_t()) {}
    
    std::vector<integer_t> indices;
    std::shared_ptr<boost::graph_t> graph;
    std::shared_ptr<boost::graph_t> subgraph;
    VTK_SMART(PolyData) polydata, filtered_polydata;
    VTK_SMART(Actor) points_actor, edges_actor, faces_actor;
    bool active;
    std::vector< std::vector<integer_t> > cc;
    std::vector<scalar_t> degrees;
    int type;
    std::string type_name() const {
        return (type>0) ? "plus type" : "minus type";
    }
};


// **************************************************************************
// 
//                      Program parameters and options
//
// **************************************************************************
namespace params {
std::string filename;       // input file name
int  verbose=0;             // output wordiness level
bool world_coords=false;    // map geocoordinates to spatial coordinates 
bool show_earth=false;      // display earth surface
bool show_axes=false;       // display coordinate axes
bool show_points=true;      // display vertices
bool show_edges=false;      // display edges
bool show_progress=false;   // show progress bar during computation
bool do_type[2] = { true, true }; // process type +1 (resp. -1) vertices
scalar_t search_radius=3.;      // radius of neighbor search region
scalar_t max_search_radius=10.; // maximum size of search area
scalar_t sphere_radius=0.5;  // radius of sphere glyphs in point depiction
integer_t n_neighbors=3;     // number of neighbors used for connectivity
integer_t min_size = 10;     // min size of connected component for display
color_t bg_color=color_t(1,1,1); // background color
}


// **************************************************************************
// 
//                      Global constants and variables
//
// **************************************************************************
const scalar_t _infinity_=std::numeric_limits<scalar_t>::max();
const integer_t invalid_index=static_cast<integer_t>(-1);
const scalar_t _core_radius_ =3486; // in km
const scalar_t _earth_radius_=6371; // in km
const bool render_now = true;
const bool delay_rendering = false;
const bool create_transfer_function = true;
const bool no_transfer_function = false;

std::vector<vertex_t>    vertices_;
std::vector<neighborhood_t> neighbors_;

ProgressDisplay progress;

TypeAttributes type_attr[2];


// **************************************************************************
// 
//                             VTK operations
//
// **************************************************************************
VTK_SMART(Actor)
get_points_actor(VTK_SMART(PolyData) pd, VTK_SMART(ColorTransferFunction) color_tf) {
    std::cout << "get_points_actor\n";
    VTK_SMART(PolyData) spheres=
        vtk_utils::make_spheres(pd, params::sphere_radius);
    VTK_MAKE_ACTOR(actor, spheres);
    actor->GetMapper()->SetLookupTable(color_tf);
    actor->GetMapper()->ScalarVisibilityOn();
    return actor;
}

VTK_SMART(Actor)
get_edges_actor(VTK_SMART(PolyData) pd, const boost::graph_t& graph, VTK_SMART(ColorTransferFunction) color_tf) {
    std::cout << "get_edges_actor\n";
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
    VTK_MAKE_ACTOR(actor, pd);
    actor->GetMapper()->SetLookupTable(color_tf);
    actor->GetMapper()->ScalarVisibilityOn();
    actor->GetProperty()->SetLineWidth(2);
    return actor;                   
}

VTK_SMART(Actor) get_earth_actor() {
    VTK_CREATE(vtkSphereSource, earth);
    earth->SetRadius(_earth_radius_);
    earth->SetCenter(0,0,0);
    earth->SetThetaResolution(50);
    earth->SetPhiResolution(50);
    earth->Update();
    VTK_MAKE_ACTOR(actor, earth->GetOutput());
    actor->GetProperty()->SetColor(0,0,0);
    actor->GetProperty()->SetOpacity(0.1);
    return actor;
}

VTK_SMART(ColorTransferFunction)
create_color_transfer_function(const std::vector<scalar_t>& values) {    
    // adaptive color map with spiral color scale
    VTK_CREATE(vtkColorTransferFunction, ctf);
    
    std::vector<color_t> colors(20);
    spurt::spiral_scale(colors, 20, 0.2);
    scalar_t dval, minval, maxval;
    if (values.size()) { 
        minval=*std::min_element(values.begin(), values.end());
        maxval=*std::max_element(values.begin(), values.end());
    }
    else {
        minval = 0;
        maxval = 1;
    }
    dval=(maxval-minval)/19;
    
    for (int i=0; i<20; ++i) {
        scalar_t value=minval + i*dval;
        color_t color=colors[i];
        ctf->AddRGBPoint(value, color[0], color[1], color[2]);
    }
    return ctf;
}

VTK_SMART(OrientationMarkerWidget) 
coordinate_axes_widget(QVTKWidget* vtk_widget) {
    color_t text_color = 
        spurt::luminosity(params::bg_color) > 0.5 ?
        spurt::black : spurt::white;
    color_t text_bg_color = 0.5*(text_color + params::bg_color);
    VTK_CREATE(vtkAxesActor, axes);
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
 
    VTK_CREATE(vtkOrientationMarkerWidget, widget);
    widget->SetOutlineColor(params::bg_color[0], params::bg_color[1], params::bg_color[2]);
    widget->SetOrientationMarker(axes);
    widget->SetInteractor( vtk_widget->GetInteractor() );
    widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    return widget;
}

VTK_SMART(PolyData) 
create_point_set(VTK_SMART(PolyData) pd, const std::vector<integer_t>& indices) {
    std::cout << "create_point_set\n";
    typedef vtk_array_traits<scalar_t>::array_type array_type;

    integer_t n_points = indices.size();
    
    VTK_CREATE(array_type, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(n_points);
    for (int i=0 ; i<n_points ; ++i) {
        pos_t x(vertices_[indices[i]].position);
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
        }
        coords->SetTuple(i, &x[0]);
    }

    VTK_CREATE(vtkPoints, points);
    points->SetData(coords);
    if (pd.Get() == NULL) {
        pd = VTK_SMART(PolyData)::New();
    }
    pd->SetPoints(points); // old points, if any, will be discarded
    return pd; // return updated polydata object
}

VTK_SMART(PolyData)
mantle_vis_renderer::assign_component_ids(VTK_SMART(PolyData) pd, const std::vector< std::vector<integer_t> >& comp, bool update_ctf) {
    std::vector<scalar_t> values(pd->GetPoints()->GetNumberOfPoints());
    for (int i=0; i<comp.size() ; ++i) {
        std::for_each(comp[i].begin(), comp[i].end(), 
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

VTK_SMART(PolyData)
mantle_vis_renderer::assign_value_by_component(VTK_SMART(PolyData) pd, const std::vector< std::vector<integer_t> >& comp) {
    std::vector<scalar_t> values(pd->GetPoints()->GetNumberOfPoints());
    std::for_each(comp.begin(), comp.end(), [&](const std::vector<int>& v)
        {
            scalar_t r = static_cast<scalar_t>(drand48());
            std::for_each(v.begin(), v.end(), 
                          [&](integer_t j) { values[j] = r; });
        });
        
    if (params::verbose) {
        std::cout << values.size() << " values added\n";
    }
    return vtk_utils::add_scalars(pd, values);
}

VTK_SMART(PolyData)
filter_polydata(TypeAttributes& ta) {
    std::vector<integer_t> valid_cc;
    for (int i=0 ; i<ta.cc.size() ; ++i) {
        if (ta.cc[i].size() >= params::min_size) valid_cc.push_back(i);
    }
    std::vector<integer_t> valid_verts;
    for (int i=0 ; i<valid_cc.size() ; ++i) {
        int j=valid_cc[i];
        std::copy(ta.cc[j].begin(), ta.cc[j].end(), 
                  std::back_inserter(valid_verts));
    }
    std::map<integer_t, integer_t> old2new;
    typedef std::map<integer_t, integer_t>::iterator iterator;
    for (int i=0 ; i<valid_verts.size() ; ++i) {
        old2new[valid_verts[i]] = i;
    }
    std::vector<integer_t> valid_edges;
    typedef boost::filtered_graph_traits_t::edge_iterator edge_iterator_t;
    edge_iterator_t ei, ei_end;
    for (boost::tie(ei, ei_end)=boost::edges(*ta.subgraph); ei != ei_end; ++ei) {
        integer_t a = boost::source(*ei, *ta.subgraph);
        integer_t b = boost::target(*ei, *ta.subgraph);
        iterator it_a = old2new.find(a);
        iterator it_b = old2new.find(b);
        if (it_a != old2new.end() && it_b != old2new.end()) {
            valid_edges.push_back(it_a->second);
            valid_edges.push_back(it_b->second);
        }
    }
    
    std::vector<pos_t> valid_pos(valid_verts.size());
    for (int i=0 ; i<valid_verts.size() ; ++i) {
        valid_pos[i] = vertices_[old2new[valid_verts[i]]].position;
    }

    if (params::verbose) {
        std::cout << valid_pos.size() << " points left after filtering\n";
        std::cout << valid_edges.size()/2 << " edges left after filtering\n";
    }
    VTK_SMART(PolyData) pd = vtk_utils::make_points(valid_pos);
    vtk_utils::add_edges_from_numbers(pd, valid_edges);
    return pd;
}


// **************************************************************************
// 
//                             Graph operations
//
// **************************************************************************
void compute_connectivity(boost::graph_t& graph, const std::vector<integer_t>& verts) {
    std::cout << "compute_connectivity\n";
    locator_t locator;
    for (integer_t i=0 ; i<verts.size() ; ++i) {
        integer_t id = verts[i];
        const vertex_t& v = vertices_[id];
        const pos_t& x = v.position;
        locator.insert(point_t(x, i));
    }
    
    graph.clear();
    
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

std::shared_ptr<boost::graph_t>
filter_by_size(boost::filtered_graph_t in) {
    std::shared_ptr<boost::graph_t> out(new boost::graph_t());
    
    // ensure that the size of each neighborhood does not exceed prescribed 
    // threshold
    integer_t n_too_many = 0;
    integer_t n_too_few = 0;
    typedef boost::filtered_graph_traits_t::vertex_iterator vertex_iterator;    
    typedef boost::filtered_graph_traits_t::out_edge_iterator edge_iterator;
    typedef boost::filtered_graph_t::edge_descriptor edge_t;
    vertex_iterator vit, vend; 
    edge_iterator eit, eend;
    std::map<scalar_t, edge_t> weight2edge;
    auto weight_map = boost::get(boost::edge_weight, in);
    for (boost::tie(vit, vend)=boost::vertices(in); vit!=vend ; ++vit) {
        weight2edge.clear();
        for (boost::tie(eit, eend)=boost::out_edges(*vit, in); eit!=eend; ++eit) {
            weight2edge[weight_map[*eit]] = *eit;
        }
        integer_t n=0;
        for (auto it=weight2edge.begin() ; 
             it!=weight2edge.end() && n<params::n_neighbors ; ++it, ++n) {
                 const edge_t& e=it->second;
                 boost::add_edge(boost::source(e, in), boost::target(e, in), 
                                 weight_map[e], *out);
        }
        if (weight2edge.size()>params::n_neighbors) ++n_too_many;
        else if (weight2edge.size()<params::n_neighbors) ++n_too_few;
    }
    
    std::cout << n_too_many << " vertices had too many neighbors included\n";
    std::cout << n_too_few << " vertices have too few neighbors\n";
    std::cout << "after pruning: filtered graph contains " << boost::num_edges(*out) << " edges\n";
    return out;
}

std::shared_ptr<boost::graph_t >
filter_connectivity(boost::graph_t& in, scalar_t radius=params::search_radius) {
    std::cout << "filter_connectivity with r=" << radius << "\n";
    boost::graph_filter_t filter(boost::get(boost::edge_weight, in), radius);
    return filter_by_size(boost::filtered_graph_t(in, filter));
}

void compute_vertex_degrees(std::vector<scalar_t>& degrees, const boost::graph_t& graph) {
    typedef boost::graph_t::vertex_iterator viterator;
    
    degrees.resize(boost::num_vertices(graph));
    std::pair<viterator, viterator> iter_range;
    boost::graph_map_t<boost::vertex_index_t>::type vertex_id_map =
        boost::get(boost::vertex_index, graph);
    for (iter_range=boost::vertices(graph); 
         iter_range.first != iter_range.second; ++iter_range.first) {
        auto v = *iter_range.first;
        integer_t i = vertex_id_map[v];
        degrees[i] = boost::out_degree(i, graph);
    }
}

void connected_components(std::vector< std::vector<integer_t> >& comp, const boost::filtered_graph_t& graph) {
    std::cout << "connected_components\n" << std::flush;
    // compute connected components in graph
    std::vector<integer_t> components(boost::num_vertices(graph));
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


// **************************************************************************
// 
//                             Visibility control
//
// **************************************************************************
void mantle_vis_renderer::hide_points(TypeAttributes& attr, bool _draw) {
    if (!this->initialized) return;
    std::cout << "hide_points(...)\n" << std::flush;
    renderer->RemoveActor(attr.points_actor);
    if (_draw) draw();
}

void mantle_vis_renderer::show_points(TypeAttributes& attr, bool _draw) {
    if (!this->initialized) return;
    std::cout << "show_points(...)\n" << std::flush;
    if (!attr.active) return;
    std::cout << "show one set of points\n";
    if (attr.points_actor.Get() == NULL) {
        std::cout << "point set does not exist\n";
        std::cout << "recreate points actor\n";
        attr.points_actor = get_points_actor(attr.polydata, color_tf);
    }
    std::cout << "add actor for one point set\n";
    renderer->AddActor(attr.points_actor);
    if (_draw) draw();
}

void mantle_vis_renderer::hide_points(bool _draw) {
    std::cout << "hide_points\n" << std::flush;
    hide_points(type_attr[0], delay_rendering);
    hide_points(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::show_points(bool _draw) {
    std::cout << "show_points\n" << std::flush;
    show_points(type_attr[0], delay_rendering);
    show_points(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::hide_edges(TypeAttributes& attr, bool _draw) {
    if (!this->initialized) return;
    std::cout << "hide_edge(...)\n" << std::flush;
    renderer->RemoveActor(attr.edges_actor);
    if (_draw) draw();
}

void mantle_vis_renderer::show_edges(TypeAttributes& attr, bool _draw) {
    if (!this->initialized) return;
    std::cout << "show_edges(...)\n" << std::flush;
    if (!attr.active) return;
    std::cout << "show one edge set\n";
    if (attr.edges_actor.Get() == NULL) {
        std::cout << "edge set does not exist\n";
        std::cout << "create a new actor for edge set\n";
        attr.edges_actor = get_edges_actor(attr.polydata, *attr.subgraph, color_tf);
    }
    std::cout << "add actor for edge set\n";
    renderer->AddActor(attr.edges_actor);
    if (_draw) draw();
}

void mantle_vis_renderer::hide_edges(bool _draw) {
    std::cout << "hide_edges\n" << std::flush;
    hide_edges(type_attr[0], delay_rendering);
    hide_edges(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::show_edges(bool _draw) {
    std::cout << "show_edges\n" << std::flush;
    show_edges(type_attr[0], delay_rendering);
    show_edges(type_attr[1], delay_rendering);
    if (_draw) draw();
}


// **************************************************************************
// 
//                             Qt Slots
//
// **************************************************************************
void mantle_vis_renderer::slot_points_check(int state) {
    bool b = (state == Qt::Checked);
    if (b != params::show_points) {
        params::show_points = b;
        if (b) show_points(delay_rendering);
        else hide_points(delay_rendering);
        draw();
    }
}

void mantle_vis_renderer::slot_edges_check(int state) {
    bool b = (state == Qt::Checked);
    if (b != params::show_edges) {
        params::show_edges = b;
        if (b) show_edges(delay_rendering);
        else hide_edges(delay_rendering);
        draw();
    }
}

void mantle_vis_renderer::slot_axes_check(int state) {
    bool b = (state == Qt::Checked);
    if (b != params::show_axes) {
        params::show_axes = b;
        if (params::show_axes) {
            axes_widget->SetEnabled( 1 );
            axes_widget->InteractiveOn();
        }
        else {
            axes_widget->SetEnabled(0);
        }
        draw();
    }
}

void mantle_vis_renderer::slot_world_check(int state) {
    // do nothing
}

bool mantle_vis_renderer::slot_type_check(int state, TypeAttributes& ta) {
    bool b = (state == Qt::Checked);
    if (b != ta.active) {
        ta.active = b;
        if (ta.active) {
            if (params::show_points) show_points(ta, delay_rendering);
            if (params::show_edges) show_edges(ta, delay_rendering);
        }
        else {
            hide_points(ta, delay_rendering);
            hide_edges(ta, delay_rendering);
        }
    }
    return (b != ta.active); // something has changed
}

void mantle_vis_renderer::slot_plus_check(int state) {
    if (slot_type_check(state, type_attr[0])) draw();
    params::do_type[0] = (state == Qt::Checked);
}

void mantle_vis_renderer::slot_minus_check(int state) {
    if (slot_type_check(state, type_attr[1])) draw();
    params::do_type[1] = (state == Qt::Checked);
}

void mantle_vis_renderer::slot_radius_slider(int r) {
    std::cout << "slot_radius_slider(" << r << ")\n";
    double radius = scalar_t(r*params::max_search_radius)/100.;
    if (radius != params::search_radius) {
        params::search_radius = radius;
        // pass on change to spinbox, which will draw
        this->controls->ui->radiusSpin->setValue(radius);
    }
}

void mantle_vis_renderer::slot_radius_spinbox(double r) {
    std::cout << "slot_radius_spinbox(" << r << ")\n";
    if (r != params::search_radius) {
        int radius = static_cast<int>(r/params::max_search_radius*100.);
        this->controls->ui->radiusSlider->setValue(radius);
        params::search_radius = r;
        draw_from_scratch(type_attr[0], delay_rendering, create_transfer_function);
        draw_from_scratch(type_attr[1], delay_rendering, no_transfer_function);
        draw();
    }
}

void mantle_vis_renderer::slot_size_spinbox(int s) {
    std::cout << "slot_size_spinbox(" << s << ")\n";
    if (s != params::n_neighbors) {
        params::n_neighbors=s;
        draw_from_scratch(type_attr[0], delay_rendering, create_transfer_function);
        draw_from_scratch(type_attr[1], delay_rendering, no_transfer_function);
        // type_attr[0].subgraph = filter_connectivity(*type_attr[0].graph);
        // type_attr[1].subgraph = filter_connectivity(*type_attr[1].graph);
        // redraw_points(false);
        // redraw_edges(true);
        draw();
    }
}

void mantle_vis_renderer::slot_min_spinbox(int m) {
    if (m != params::min_size) {
        params::min_size = m;
    }
}

void mantle_vis_renderer::slot_file_action(QAction* action) {
    if (action->text() == "Save Current Frame") {
        save_snapshot();
    }
    else if (action->text() == "Load Camera Settings") {
        camera_in();
    }
    else if (action->text() == "Save Camera Settings") {
        camera_out();
    }
    else if (action->text() == "Quit") {
        this->close();
    }
}


// **************************************************************************
// 
//                             I/O Operations
//
// **************************************************************************

bool init(int argc, char** argv) {
    std::cout << "init\n" << std::flush;
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
        parser.add_value("min", params::min_size, 
                         "Min size of connected component for display",
                         optional_group);
        parser.add_value("points", params::show_points,
                        "Show vertices", optional_group);
        parser.add_value("edges", params::show_edges,
                        "Show edges", optional_group);
        parser.add_value("positive", params::do_type[0],
                        "Process type +1 vertices", optional_group);
        parser.add_value("negative", params::do_type[1],
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

bbox_t import_vertices(const std::string& data_name) {
    std::cout << "import_vertices\n";
    // load vertices from file
    std::vector<vertex_t> __vertices;
    bbox_t bounds = 
        core_mantle::read_text(__vertices, params::filename, false);
    
    // filter out redundancy
    vertices_.clear();
    type_attr[0].indices.clear();
    type_attr[1].indices.clear();
    std::set<pos_t, less_pos_tol> unique_pos;
    integer_t n_points = 0;
    for (integer_t i=0; i<__vertices.size(); ++i) {
        const vertex_t v = __vertices[i];
        const pos_t& x=v.position;
        if (unique_pos.find(x) == unique_pos.end()) {
            unique_pos.insert(x);
            vertices_.push_back(v);
            if (v.type>0) type_attr[0].indices.push_back(n_points++);
            else type_attr[1].indices.push_back(n_points++);
        } // else skip redundant vertex
    }
    
    if (params::verbose) {
        std::cout 
            << __vertices.size() << " vertices imported from file\n"
            << vertices_.size() << " unique vertices\n"
            << type_attr[0].indices.size() << " type +1\n"
            << type_attr[1].indices.size() << " type -1\n"
            << "bounding box:\n" << bounds << '\n';
    }
    
    return bounds;
}

void mantle_vis_renderer::camera_out() {
    QString fn = 
        QFileDialog::getSaveFileName(this, 
                                     QString("Select an output camera file"),
                                     QString(home_directory.c_str()));
    std::string filename = fn.toStdString();
    vtk_utils::export_camera_settings(filename, renderer);
}

void mantle_vis_renderer::camera_in() {
    QString fn = 
        QFileDialog::getOpenFileName(this,
                                     QString("Select an input camera file"),
                                     QString(home_directory.c_str()));
    std::string filename = fn.toStdString();
    vtk_utils::import_camera_settings(filename, renderer);
    draw();
}

void mantle_vis_renderer::load_file(const std::string& data_name) {
    import_vertices(data_name);
    internal_main();
}

void mantle_vis_renderer::save_snapshot() {
    QString fn = 
        QFileDialog::getSaveFileName(this, QString("Select file name to save snapshot"),
                                     QString(home_directory.c_str()).
                                        append("/Desktop/untitled.png"));
    std::string filename = fn.toStdString();
    vtk_utils::save_frame(this->ui->qvtkWidget->GetRenderWindow(), 
                          filename);
}


// **************************************************************************
// 
//                             Drawing
//
// **************************************************************************
void mantle_vis_renderer::redraw_edges(bool _draw) {
    std::cout << "redraw_edges\n" << std::flush;
    // if type is active: actor must be hidden, deleted, redrawn, and shown
    // if type is inactive: actor must be deleted and redrawn
    bool fb_changed = false;
    for (int i=0 ; i<2 ; ++i) {
        if (type_attr[i].active) {
            std::cout << "type " << i << " is active\n";
            std::cout << "hide edges for type " << i << '\n';
            hide_edges(type_attr[i], delay_rendering);
            std::cout << "delete edges for type " << i << '\n';
            type_attr[i].edges_actor = NULL;
            std::cout << "show edges for type " << i << '\n';
            show_edges(type_attr[i], delay_rendering);
            fb_changed = true;
        }
        else {
            std::cout << "type " << i << " is inactive\n";
            std::cout << "delete edges for type " << i << '\n';
            type_attr[i].edges_actor = NULL; // we will draw it once needed
        }
    }
    if (fb_changed && _draw) draw();
}

void mantle_vis_renderer::redraw_points(bool _draw) {
    std::cout << "redraw_points\n" << std::flush;
    // if type is active: points must be hidden, deleted, redrawn, and shown
    // if type is inactive: actor must be deleted and redrawn
    bool fb_changed = false;
    for (int i=0 ; i<2 ; ++i) {
        if (type_attr[i].active) {
            std::cout << "type " << i << " is active\n";
            hide_points(type_attr[i], delay_rendering);
            std::cout << "delete point set " << i << '\n';
            type_attr[i].points_actor = NULL;
            show_points(type_attr[i], delay_rendering);
            fb_changed = true;
        }
        else {
            std::cout << "type " << i << " is inactive\n";
            if (type_attr[i].points_actor.Get() != NULL) {
                std::cout << "point set " << i << " exists\n";
                type_attr[i].points_actor = NULL;
                std::cout << "delete point set " << i << '\n';
                std::cout << "create new actor for point set " << i << '\n';
                type_attr[i].points_actor = get_points_actor(type_attr[i].polydata, color_tf);
            } // else we leave it NULL and will draw it once needed
        }
    }
    if (fb_changed && _draw) draw();
}

void mantle_vis_renderer::draw() {
    if (this->initialized) 
        this->ui->qvtkWidget->GetRenderWindow()->Render();
}

void mantle_vis_renderer::internal_main() {
    std::cout << "internal_main\n" << std::flush;
    
    import_vertices(params::filename);
    
    // create coordinate axes
    axes_widget = coordinate_axes_widget(this->ui->qvtkWidget);
    
    for (int i=0 ; i<2 ; ++i) {
        type_attr[i].active = params::do_type[i];
        type_attr[i].graph = std::shared_ptr<boost::graph_t>(new boost::graph_t());
    
        // determine superset of type-specific vertex connectivity
        boost::edge_map_t edge_weight;
        edge_weight=get(boost::edge_weight, *type_attr[i].graph);
        compute_connectivity(*type_attr[i].graph, type_attr[i].indices);
        draw_from_scratch(type_attr[i], delay_rendering, !i);
    }
    
    if (params::show_axes) {
        axes_widget->SetEnabled( 1 );
        axes_widget->InteractiveOn();
    }

    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    draw();
}

void mantle_vis_renderer::
draw_from_scratch(TypeAttributes& ta, bool _draw, bool create_ctf) {
    std::cout << "draw_from_scratch(" << ta.type_name() << ")\n" << std::flush;
    
    std::cout << "\t1. there are currently " 
        << renderer->GetActors()->GetNumberOfItems() << " actors\n";
    
    hide_points(ta);
    hide_edges(ta);
    
    std::cout << "\t2. there are currently " 
        << renderer->GetActors()->GetNumberOfItems() << " actors\n";
        
    // identify subset of connectivity corresponding to chosen radius
    ta.subgraph = filter_connectivity(*ta.graph);

    if (params::verbose) {
        std::cout << "\t\t" << boost::num_edges(*ta.subgraph) 
                  << " edges selected for " << ta.type_name() << "\n";
    }

    if (false) {
        compute_vertex_degrees(ta.degrees, *ta.subgraph);
    }

    // initialize connected components
    connected_components(ta.cc, *ta.subgraph);

    if (params::verbose) {
        std::cout << "\t\tThere are " << ta.cc.size() 
            << " connected components for " << ta.type_name() << "\n";
    }

    // Create polydata to be visualized
    ta.polydata = VTK_SMART(PolyData)::New();
    ta.polydata = create_point_set(ta.polydata, ta.indices);

    // assign connected component id to each vertex
    /*ta.polydata = assign_component_ids(ta.polydata, ta.cc, create_ctf);*/
    // assign random color per connected component
    ta.polydata = assign_value_by_component(ta.polydata, ta.cc);
    
    if (ta.active) {
        if (params::show_points) show_points(ta, delay_rendering);
    
    std::cout << "\t3. there are currently " 
        << renderer->GetActors()->GetNumberOfItems() << " actors\n";
        if (params::show_edges) show_edges(ta, delay_rendering);

    
        std::cout << "\t4. there are currently " 
            << renderer->GetActors()->GetNumberOfItems() << " actors\n";
    }
    
    std::cout << "\t5. there are currently " 
        << renderer->GetActors()->GetNumberOfItems() << " actors\n";

    if (_draw) {
        renderer->ResetCamera();
        renderer->ResetCameraClippingRange();
        draw();
    }
}


// **************************************************************************
// 
//                             Constructors / Destructors
//
// **************************************************************************
mantle_vis_renderer::mantle_vis_renderer(int argc, char* argv[]) : initialized(false) {
    std::cout << "constructor\n" << std::flush;
        
    if (!init(argc, argv)) {
        std::cout << "Invalid command line options. Exiting.\n";
        exit(1);
    }
    
    type_attr[0].type = 1;
    type_attr[1].type = -1;
    
    srand48(time(NULL));
    
    progress = ProgressDisplay(params::show_progress);
    
    this->ui = new Ui_MainWindow;
    this->ui->setupUi(this);
    
    // use home directory as default location for file browsing
    home_directory = getpwuid(getuid())->pw_dir;
    
    // Initialize visualization pipeline and interactor
    renderer = VTK_SMART(Renderer)::New();
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(800, 800);
                               
    this->ui->qvtkWidget->SetRenderWindow(window);
    
    this->controls = new mantle_vis_control();
    this->controls->ui = new Ui_ControlWidget;
    this->controls->ui->setupUi(this->controls);
    this->controls->ui->pointsCheck->setTristate(false);
    this->controls->ui->edgesCheck->setTristate(false);
    this->controls->ui->axesCheck->setTristate(false);
    this->controls->ui->minusCheck->setTristate(false);
    this->controls->ui->plusCheck->setTristate(false);
    
    // create connections
    connect(this->controls->ui->pointsCheck, SIGNAL(stateChanged(int)),
            this, SLOT(slot_points_check(int)));
    connect(this->controls->ui->edgesCheck, SIGNAL(stateChanged(int)), 
            this, SLOT(slot_edges_check(int)));
    connect(this->controls->ui->axesCheck, SIGNAL(stateChanged(int)),
            this, SLOT(slot_axes_check(int)));
    connect(this->controls->ui->minusCheck, SIGNAL(stateChanged(int)), 
            this, SLOT(slot_minus_check(int)));
    connect(this->controls->ui->plusCheck, SIGNAL(stateChanged(int)), 
            this, SLOT(slot_plus_check(int)));
    connect(this->controls->ui->radiusSpin, SIGNAL(valueChanged(double)), 
            this, SLOT(slot_radius_spinbox(double)));
    connect(this->controls->ui->radiusSlider, SIGNAL(valueChanged(int)), 
            this, SLOT(slot_radius_slider(int)));
    connect(this->controls->ui->sizeSpin, SIGNAL(valueChanged(int)),
            this, SLOT(slot_size_spinbox(int)));
    connect(this->ui->menuFile, SIGNAL(triggered(QAction*)),
            this, SLOT(slot_file_action(QAction*)));
            
    // initialize control values with parameters
    this->controls->ui->pointsCheck->setCheckState(params::show_points ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->edgesCheck->setCheckState(params::show_edges ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->axesCheck->setCheckState(params::show_axes ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->plusCheck->setCheckState(params::do_type[0] ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->minusCheck->setCheckState(params::do_type[1] ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->radiusSpin->setValue(params::search_radius);
    this->controls->ui->radiusSlider->setValue(static_cast<int>(params::search_radius/params::max_search_radius*100.));
    this->controls->ui->sizeSpin->setValue(params::n_neighbors);
        
    // configure rendering window
    renderer->SetBackground(params::bg_color[0],
                            params::bg_color[1],
                            params::bg_color[2]);
                            
    // create default color transfer function
    color_tf = create_color_transfer_function(std::vector<scalar_t>());
    
    this->initialized = true;
    internal_main();
    
    this->controls->show();
}

mantle_vis_control::mantle_vis_control() {
    this->ui = new Ui_ControlWidget;
    this->ui->setupUi(this);
}

mantle_vis_control::~mantle_vis_control() {
    delete this->ui;
}

mantle_vis_renderer::~mantle_vis_renderer() {
    delete this->controls;
    delete this->ui;
}