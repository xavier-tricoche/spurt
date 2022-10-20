#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <string>


#include "core_mantle.hpp"
#include "core_mantle_io.hpp"
#include "utils.hpp"

// spurt's utilities
#include <format/format.hpp>
#include <misc/option_parse.hpp>
#include <learning/isomap.hpp>
#include <graph/connectivity.hpp>

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
#include "QVTKWidget.h"
#include "QVTKInteractor.h"
#include "mantle_vis_control.hpp"
#include "mantle_vis_renderer.hpp"

// Tapkee
// #include <tapkee/tapkee.hpp>
// #include <tapkee/methods.hpp>

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

typedef nvis::fixed_vector<vtkIdType, 3> triangle_t;

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
    return get(m_weight, e)<m_threshold;
  }
  
  EdgeWeightMap m_weight;
  scalar_t m_threshold;
};

namespace boost {
typedef property<vertex_color_t, bool> vertex_prop_t;
typedef property<edge_weight_t, scalar_t> edge_prop_t;
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

typedef spurt::Isomap<scalar_t, boost::graph_t> dim_reduct_t;

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
bool show_cells=false;      // display cells
bool show_progress=false;   // show progress bar during computation
bool do_type[2]={ true, true }; // process type +1 (resp. -1) vertices
scalar_t search_radius=3.;      // radius of neighbor search region
scalar_t max_search_radius=10.; // maximum size of search area
scalar_t sphere_radius=0.5;  // radius of sphere glyphs in point depiction
integer_t n_neighbors=3;     // number of neighbors used for connectivity
integer_t min_size=10;     // min size of connected component for display
integer_t edge_line_width=2; // line width to draw edges
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
const bool delay_rendering=false;
const bool create_transfer_function=true;
const bool no_transfer_function=false;

std::vector<vertex_t>    vertices_;
std::vector<neighborhood_t> neighbors_;

ProgressDisplay progress;

struct Mesh {
    Mesh() 
        : geometry(NULL), points_actor(NULL), edges_actor(NULL), 
          cells_actor(NULL) {}
    Mesh(const Mesh& other) 
        : geometry(other.geometry), points_actor(other.points_actor),
            edges_actor(other.edges_actor), cells_actor(other.cells_actor) {}
    
    VTK_SMART(PolyData) geometry;
    VTK_SMART(Actor) points_actor, edges_actor, cells_actor;
    color_t color;
};

struct spurt::gmig::TypeAttributes {
    TypeAttributes() : graph(new boost::graph_t()) {}
    
    std::vector<pos_t> positions;
    std::shared_ptr<boost::graph_t> graph, subgraph;
    std::vector<Mesh> meshes;
    bool active;
    std::vector< std::vector<integer_t> > cc;
    std::vector<scalar_t> degrees;
    int type;
    std::string type_name() const {
        return (type>0) ? "plus type" : "minus type";
    }
    const pos_t& position(integer_t i) const { 
        return positions[i];
    }
};

TypeAttributes type_attr[2];

vtkPolyData* dummy;


// **************************************************************************
// 
//                             VTK operations
//
// **************************************************************************
VTK_SMART(Actor)
get_points_actor(VTK_SMART(PolyData) pd, const color_t& color) {
    VTK_SMART(PolyData) spheres=
        vtk_utils::make_spheres(pd, params::sphere_radius, 6, 6);
    VTK_MAKE_ACTOR(actor, spheres);
    //actor->GetMapper()->SetLookupTable(color_tf);
    actor->GetMapper()->ScalarVisibilityOff();
    actor->GetProperty()->SetColor(color[0], color[1], color[2]);
    return actor;
}

VTK_SMART(Actor)
get_edges_actor(VTK_SMART(PolyData) pd, const color_t& color) {
    std::cout << "get_edges_actor\n";
    
    VTK_CREATE(vtkExtractEdges, edges);
    edges->SetInputData(pd);
    edges->Update();
    VTK_MAKE_ACTOR(actor, edges->GetOutput());
    //actor->GetMapper()->SetLookupTable(color_tf);
    actor->GetMapper()->ScalarVisibilityOff();
    actor->GetProperty()->SetColor(color[0], color[1], color[2]);
    return actor;                   
}

VTK_SMART(Actor)
get_cells_actor(VTK_SMART(PolyData) pd, const color_t& color) {
    std::cout << "get_cells_actor\n";
    
    VTK_MAKE_ACTOR(actor, pd);
    //actor->GetMapper()->SetLookupTable(color_tf);
    actor->GetMapper()->ScalarVisibilityOff();
    actor->GetProperty()->SetColor(color[0], color[1], color[2]);
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
        minval=0;
        maxval=1;
    }
    dval=(maxval-minval)/19;
    
    for (int i=0; i<20; ++i) {
        scalar_t value=minval+i*dval;
        color_t color=colors[i];
        ctf->AddRGBPoint(value, color[0], color[1], color[2]);
    }
    return ctf;
}

VTK_SMART(OrientationMarkerWidget) 
coordinate_axes_widget(QVTKWidget* vtk_widget) {
    color_t text_color=
        spurt::luminosity(params::bg_color)>0.5 ?
        spurt::black : spurt::white;
    color_t text_bg_color=0.5*(text_color+params::bg_color);
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

inline pos_t& coordinate_transform(pos_t& x) {
    if (!params::world_coords) {
        // rescaling to improve aspect ratio
        x[1] *= 2;
        x[2] /= 3;
    }
    else {
        scalar_t r=_core_radius_+x[2];
        scalar_t la=deg2rad(x[0]);
        scalar_t lo=deg2rad(x[1]);
        x=pos_t(r*sin(la)*cos(lo),
                  r*sin(la)*sin(lo),
                  r*cos(la));
    }
    return x;
}

VTK_SMART(PolyData) 
create_point_set(VTK_SMART(PolyData) pd, const std::vector<integer_t>& indices) {
    std::cout << "create_point_set\n";
    typedef vtk_array_traits<scalar_t>::array_type array_type;

    integer_t n_points=indices.size();
    
    VTK_CREATE(array_type, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(n_points);
    for (int i=0; i<n_points; ++i) {
        pos_t x(vertices_[indices[i]].position);
        x=coordinate_transform(x);
        coords->SetTuple(i, &x[0]);
    }

    VTK_CREATE(vtkPoints, points);
    points->SetData(coords);
    if (pd.Get() == NULL) {
        pd=VTK_SMART(PolyData)::New();
    }
    pd->SetPoints(points); // old points, if any, will be discarded
    return pd; // return updated polydata object
}

VTK_SMART(PolyData)
mantle_vis_renderer::assign_component_ids(VTK_SMART(PolyData) pd, const std::vector< std::vector<integer_t> >& comp, bool update_ctf) {
    std::vector<scalar_t> values(pd->GetPoints()->GetNumberOfPoints());
    for (int i=0; i<comp.size(); ++i) {
        std::for_each(comp[i].begin(), comp[i].end(), 
                      [&](integer_t j) { values[j]=i; });
    }
    // Update color mapping
    if (update_ctf) {
        color_tf=create_color_transfer_function(values);
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
            scalar_t r=static_cast<scalar_t>(drand48());
            std::for_each(v.begin(), v.end(), 
                          [&](integer_t j) { values[j]=r; });
        });
        
    if (params::verbose) {
        std::cout << values.size() << " values added\n";
    }
    return vtk_utils::add_scalars(pd, values);
}

void triangulate(std::vector<triangle_t>&, boost::graph_t&, const std::vector<integer_t>&);

void filter_polydata(TypeAttributes& ta) {
    std::cout << "filter_polydata\n";
    
    ta.meshes.clear();
    scalar_t min=_infinity_, max=0;
    for (integer_t i=0; i<ta.cc.size(); ++i) {
        const std::vector<integer_t>& comp=ta.cc[i];
        integer_t n=comp.size();
        if (n<params::min_size) continue;
        
        Mesh mesh;

        std::map<integer_t, integer_t> old2new;
        typedef std::map<integer_t, integer_t>::iterator iterator;
        std::vector<pos_t> pos(n);
        for (int j=0; j<n; ++j) {
            old2new[comp[j]]=j;
            pos[j]=ta.position(comp[j]);
        }
        mesh.geometry=vtk_utils::make_points(pos);
        
        std::vector<triangle_t> triangles;
        triangulate(triangles, *ta.subgraph, comp);
        VTK_CREATE(vtkCellArray, cells);
        for (int j=0 ; j<triangles.size() ; ++j) {
            triangle_t tri;
            tri[0]=old2new[triangles[j][0]];
            tri[1]=old2new[triangles[j][1]];
            tri[2]=old2new[triangles[j][2]];
            if (tri[0]<0 || tri[0]>=comp.size()) std::cout << "Invalid triangle index\n";
            if (tri[1]<0 || tri[1]>=comp.size()) std::cout << "Invalid triangle index\n";
            if (tri[2]<0 || tri[2]>=comp.size()) std::cout << "Invalid triangle index\n";
            
            cells->InsertNextCell(3, tri.begin());
        }
        mesh.geometry->SetPolys(cells);
        ta.meshes.push_back(mesh);
        if (pos.size()<min) min=pos.size();
        if (pos.size()>max) max=pos.size();
    }
    for (int i=0; i<ta.meshes.size(); ++i) {
        scalar_t n=ta.meshes[i].geometry->GetNumberOfPoints();
        scalar_t u=(n-min)/max;
        //ta.meshes[i].color=color_t(u, u, 1-u);
        ta.meshes[i].color=color_t(drand48(), drand48(), drand48());
    }
}


// **************************************************************************
// 
//                             Graph operations
//
// **************************************************************************
   
// ensure that the size of each vertex neighborhood does not exceed 
// prescribed threshold
std::shared_ptr<boost::graph_t>
filter_by_size(TypeAttributes& ta, boost::filtered_graph_t& in) {
    std::shared_ptr<boost::graph_t> out(new boost::graph_t());
    integer_t n_too_many=0;
    integer_t n_too_few=0;
    typedef boost::filtered_graph_traits_t::vertex_iterator vertex_iterator;    
    typedef boost::filtered_graph_traits_t::out_edge_iterator edge_iterator;
    typedef boost::filtered_graph_t::edge_descriptor edge_t;
    // typedef boost::filtered_graph_map<vertex_index_t>::type index_map_t;
    vertex_iterator vit, vend; 
    edge_iterator eit, eend;
    std::map<scalar_t, edge_t> weight2edge;
    auto weight_map=boost::get(boost::edge_weight, in);
    scalar_t min_length, max_length;
    min_length=std::numeric_limits<scalar_t>::max();
    max_length=std::numeric_limits<scalar_t>::min();
    for (boost::tie(vit, vend)=boost::vertices(in); vit!=vend; ++vit) {
        weight2edge.clear();
        for (boost::tie(eit, eend)=boost::out_edges(*vit, in); eit!=eend; ++eit) {            
            weight2edge[weight_map[*eit]]=*eit;
        }
        integer_t n=0;
        for (auto it=weight2edge.begin(); 
             it!=weight2edge.end() && n<params::n_neighbors; ++it, ++n) {
             const edge_t& e=it->second;
             boost::add_edge(boost::source(e, in), boost::target(e, in), 
                             weight_map[e], *out);
             scalar_t w=weight_map[e];
             if (w<min_length) min_length=w;
             if (w>max_length) max_length=w;
        }
        if (weight2edge.size()>params::n_neighbors) ++n_too_many;
        else if (weight2edge.size()<params::n_neighbors) ++n_too_few;
    }
    
    std::cout << n_too_many << " vertices had too many neighbors included\n";
    std::cout << n_too_few << " vertices have too few neighbors\n";
    std::cout << "after pruning: filtered graph contains " << boost::num_edges(*out) << " edges\n";
    std::cout << "after pruning: edge lengths in [" << min_length
        << ", " << max_length << "]\n";

    return out;
}

// prune edges whose length is above prescribed threshold
std::shared_ptr<boost::graph_t >
filter_connectivity(TypeAttributes& ta, scalar_t radius=params::search_radius) {
    std::cout << "filter_connectivity with r=" << radius << "\n";
    
    boost::graph_t& in=*ta.graph;
    //
    // typedef boost::graph_traits_t::vertex_iterator vertex_iterator;
    // typedef boost::graph_traits_t::out_edge_iterator edge_iterator;
    // typedef boost::graph_t::edge_descriptor edge_t;
    // vertex_iterator vit, vend;
    // edge_iterator eit, eend;
    // std::map<scalar_t, edge_t> weight2edge;
    // auto weight_map=boost::get(boost::edge_weight, in);
    // auto index_map=boost::get(boost::vertex_index, in);
    // for (boost::tie(vit, vend)=boost::vertices(in); vit!=vend; ++vit) {
    //     weight2edge.clear();
    //     for (boost::tie(eit, eend)=boost::out_edges(*vit, in); eit!=eend; ++eit) {
    //         weight2edge[weight_map[*eit]]=*eit;
    //         integer_t a=index_map[boost::source(*eit, in)];
    //         integer_t b=index_map[boost::target(*eit, in)];
    //         if (a==0) {
    //             std::cout << "|edge (" << a << ", " << b << ")|="
    //                 << nvis::norm(ta.position(a)-ta.position(b))
    //                     << " (=" << weight_map[*eit] << "\?\?)\n";
    //         }
    //     }
    // }
    
    
    boost::graph_filter_t filter(boost::get(boost::edge_weight, in), radius);
    boost::filtered_graph_t filtered_by_diameter(in, filter);
    return filter_by_size(ta, filtered_by_diameter);
}

void compute_vertex_degrees(std::vector<scalar_t>& degrees, const boost::graph_t& graph) {
    typedef boost::graph_t::vertex_iterator viterator;
    
    degrees.resize(boost::num_vertices(graph));
    std::pair<viterator, viterator> iter_range;
    boost::graph_map_t<boost::vertex_index_t>::type vertex_id_map =
        boost::get(boost::vertex_index, graph);
    for (iter_range=boost::vertices(graph); 
         iter_range.first!=iter_range.second; ++iter_range.first) {
        auto v=*iter_range.first;
        integer_t i=vertex_id_map[v];
        degrees[i]=boost::out_degree(i, graph);
    }
}


// **************************************************************************
// 
//                             Meshing
//
// **************************************************************************

struct vertex_predicate {
    vertex_predicate(const std::vector<integer_t>& _ids) : ids() {
        std::copy(_ids.begin(), _ids.end(), std::inserter(ids, ids.end()));
    }
    
    template <typename Vertex_>
    bool operator()(const Vertex_& v) const {
        return ids.find(v) != ids.end();
    }
    
    std::set<integer_t> ids;
};

size_t ntris=0;
void triangulate(std::vector<triangle_t>& triangles, boost::graph_t& graph, const std::vector<integer_t>& ids) {
    typedef boost::graph_traits_t::out_edge_iterator edge_iterator;
    
    // std::cout << "triangulate #1: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
        
    boost::graph_t subgraph;
    // std::cout << "extracting subgraph for " << ids.size() << " vertices\n" << std::flush;
    // std::cout << "Those indices are:\n";
    // std::copy(ids.begin(), ids.end(), std::ostream_iterator<integer_t>(std::cout, ", "));
    // std::cout << "\ninput graph contains " << boost::num_vertices(graph) << " vertices and "
    //     << boost::num_edges(graph) << " edges\n";
    spurt::extract_subgraph(subgraph, graph, vertex_predicate(ids));
    // std::cout << "resulting subgraph contains " << boost::num_vertices(subgraph) << " vertices and "
        // << boost::num_edges(subgraph) << " edges\n" << std::flush;
    
    // std::cout << "triangulate #2 " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
    
    dim_reduct_t isomap(subgraph);
    
    // std::cout << "triangulate #3: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
    isomap.embed();
    // std::cout << "triangulate #4: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";

    typename dim_reduct_t::matrix_t coordinates=isomap.coordinates(2); //(Nxdim)
    std::vector<pos_t> coords(ids.size(), 0.);
    for (int i=0; i<ids.size() ; ++i) {
        coords[i][0] = coordinates(i,0);
        coords[i][1] = coordinates(i,1);
    }
    // std::cout << "triangulate #5: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
    VTK_SMART(PolyData) pd=vtk_utils::make_points(coords);
    // std::cout << "triangulate #6: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
    vtk_utils::add_mesh2d(pd);
    // std::cout << "triangulate #7: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
    VTK_SMART(CellArray) tris(pd->GetPolys());
    vtkIdType* tri_ids = tris->GetPointer();
    size_t N = tris->GetNumberOfCells();
    triangles.resize(N);
    int j=0;
    for (int i=0 ; i<N ; ++i) {
        triangles[i][0] = tri_ids[j++];
        triangles[i][1] = tri_ids[j++];
        triangles[i][2] = tri_ids[j++];
        if (i<10) std::cout << "Triangle " << i<< ": " << tri_ids[0] << "," << tri_ids[1] <<","<<tri_ids[2] << "\n";
    }
    // std::cout << "triangulate #8: " << dummy->GetVerts()->GetNumberOfCells() << " verts\n";
    std::cout << "triangulation of " << ids.size() << " vertices produced "
        << triangles.size() << " triangles\n";
    ntris+=triangles.size();
    std::cout << ntris << " total triangles so far\n";
} 

// **************************************************************************
// 
//                             Visibility control
//
// **************************************************************************
// Points
void mantle_vis_renderer::hide_points(TypeAttributes& ta, bool _draw) {
    if (!this->initialized) return;
    std::cout << "hide_points(...)\n" << std::flush;
    for (int i=0; i<ta.meshes.size(); ++i) {
        renderer->RemoveActor(ta.meshes[i].points_actor);
    }
    if (_draw) draw();
}

void mantle_vis_renderer::hide_points(bool _draw) {
    std::cout << "hide_points\n" << std::flush;
    hide_points(type_attr[0], delay_rendering);
    hide_points(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::show_points(TypeAttributes& ta, bool _draw) {
    if (!this->initialized) return;
    std::cout << "show_points(...)\n" << std::flush;
    if (!ta.active) return;
    std::cout << "show one set of points\n";
    for (int i=0; i<ta.meshes.size(); ++i) {
        if (ta.meshes[i].points_actor.Get()==NULL) {
            std::cout << "point set does not exist\n";
            std::cout << "recreate points actor\n";
            ta.meshes[i].points_actor=
                get_points_actor(ta.meshes[i].geometry, ta.meshes[i].color);
        }
        std::cout << "add actor for one point set\n";
        renderer->AddActor(ta.meshes[i].points_actor);
    }
    if (_draw) draw();
}

void mantle_vis_renderer::show_points(bool _draw) {
    std::cout << "show_points\n" << std::flush;
    show_points(type_attr[0], delay_rendering);
    show_points(type_attr[1], delay_rendering);
    if (_draw) draw();
}

// Edges
void mantle_vis_renderer::hide_edges(TypeAttributes& ta, bool _draw) {
    if (!this->initialized) return;
    std::cout << "hide_edge(...)\n" << std::flush;
    for (int i=0; i<ta.meshes.size(); ++i) {
        renderer->RemoveActor(ta.meshes[i].edges_actor);
    }
    std::cout << "edges actors for " << ta.type_name() << " have been removed\n";
   
    if (_draw) draw();
}

void mantle_vis_renderer::hide_edges(bool _draw) {
    std::cout << "hide_edges\n" << std::flush;
    hide_edges(type_attr[0], delay_rendering);
    hide_edges(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::show_edges(TypeAttributes& ta, bool _draw) {
    if (!this->initialized) return;
    std::cout << "show_edges(...)\n" << std::flush;
    if (!ta.active) return;
    std::cout << "show one edge set\n";
    for (int i=0; i<ta.meshes.size(); ++i) {
        if (ta.meshes[i].edges_actor.Get()==NULL) {
            std::cout << "edge set does not exist\n";
            std::cout << "create a new actor for edge set\n";
            ta.meshes[i].edges_actor=get_edges_actor(ta.meshes[i].geometry, ta.meshes[i].color);
        }
        std::cout << "add actor for edge set\n";
        renderer->AddActor(ta.meshes[i].edges_actor);
    }
    if (_draw) draw();
}

void mantle_vis_renderer::show_edges(bool _draw) {
    std::cout << "show_edges\n" << std::flush;
    show_edges(type_attr[0], delay_rendering);
    show_edges(type_attr[1], delay_rendering);
    if (_draw) draw();
}

// Cells
void mantle_vis_renderer::hide_cells(TypeAttributes& ta, bool _draw) {
    if (!this->initialized) return;
    std::cout << "hide_cells(...)\n" << std::flush;
    for (int i=0; i<ta.meshes.size(); ++i) {
        renderer->RemoveActor(ta.meshes[i].cells_actor);
    }
    std::cout << "cells actors for " << ta.type_name() << " have been removed\n";
   
    if (_draw) draw();
}

void mantle_vis_renderer::hide_cells(bool _draw) {
    std::cout << "hide_edges\n" << std::flush;
    hide_edges(type_attr[0], delay_rendering);
    hide_edges(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::show_cells(TypeAttributes& ta, bool _draw) {
    if (!this->initialized) return;
    std::cout << "show_cells(...)\n" << std::flush;
    if (!ta.active) return;
    std::cout << "show one cell set\n";
    for (int i=0; i<ta.meshes.size(); ++i) {
        if (ta.meshes[i].cells_actor.Get()==NULL) {
            std::cout << "edge set does not exist\n";
            std::cout << "create a new actor for edge set\n";
            ta.meshes[i].cells_actor=get_cells_actor(ta.meshes[i].geometry, ta.meshes[i].color);
        }
        std::cout << "add actor for cell set\n";
        renderer->AddActor(ta.meshes[i].cells_actor);
    }
    if (_draw) draw();
}

void mantle_vis_renderer::show_cells(bool _draw) {
    std::cout << "show_cells\n" << std::flush;
    show_cells(type_attr[0], delay_rendering);
    show_cells(type_attr[1], delay_rendering);
    if (_draw) draw();
}

void mantle_vis_renderer::show_hide_control(QAction* action) {
    if (this->controls->isVisible()) {
        this->controls->setVisible(false);
        action->setText(QString("Show control window"));
    }
    else {
        this->controls->setVisible(true);
        action->setText(QString("Hide control window"));
    }
}

// **************************************************************************
// 
//                             Qt Slots
//
// **************************************************************************
void mantle_vis_renderer::slot_points_check(int state) {
    bool b=(state == Qt::Checked);
    if (b!=params::show_points) {
        params::show_points=b;
        if (b) show_points(delay_rendering);
        else hide_points(delay_rendering);
        draw();
    }
}

void mantle_vis_renderer::slot_edges_check(int state) {
    bool b=(state == Qt::Checked);
    if (b!=params::show_edges) {
        params::show_edges=b;
        if (b) show_edges(delay_rendering);
        else hide_edges(delay_rendering);
        draw();
    }
}

void mantle_vis_renderer::slot_cells_check(int state) {
    bool b=(state == Qt::Checked);
    if (b!=params::show_cells) {
        params::show_cells=b;
        if (b) show_cells(delay_rendering);
        else hide_cells(delay_rendering);
        draw();
    }
}

void mantle_vis_renderer::slot_axes_check(int state) {
    bool b=(state == Qt::Checked);
    if (b!=params::show_axes) {
        params::show_axes=b;
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
    bool b=(state == Qt::Checked);
    if (b!=ta.active) {
        ta.active=b;
        if (ta.active) {
            if (params::show_points) show_points(ta, delay_rendering);
            if (params::show_edges) show_edges(ta, delay_rendering);
        }
        else {
            hide_points(ta, delay_rendering);
            hide_edges(ta, delay_rendering);
        }
    }
    return (b!=ta.active); // something has changed
}

void mantle_vis_renderer::slot_plus_check(int state) {
    if (slot_type_check(state, type_attr[0])) draw();
    params::do_type[0]=(state == Qt::Checked);
}

void mantle_vis_renderer::slot_minus_check(int state) {
    if (slot_type_check(state, type_attr[1])) draw();
    params::do_type[1]=(state == Qt::Checked);
}

void mantle_vis_renderer::slot_radius_slider(int r) {
    std::cout << "slot_radius_slider(" << r << ")\n";
    double radius=scalar_t(r*params::max_search_radius)/100.;
    if (radius!=params::search_radius) {
        params::search_radius=radius;
        // pass on change to spinbox, which will draw
        this->controls->ui->radiusSpin->setValue(radius);
    }
}

void mantle_vis_renderer::slot_radius_spinbox(double r) {
    std::cout << "slot_radius_spinbox(" << r << ")\n";
    if (r!=params::search_radius) {
        int radius=static_cast<int>(r/params::max_search_radius*100.);
        this->controls->ui->radiusSlider->setValue(radius);
        params::search_radius=r;
        draw_from_scratch(type_attr[0], delay_rendering, create_transfer_function);
        draw_from_scratch(type_attr[1], delay_rendering, no_transfer_function);
        draw();
    }
}

void mantle_vis_renderer::slot_size_spinbox(int s) {
    std::cout << "slot_size_spinbox(" << s << ")\n";
    if (s!=params::n_neighbors) {
        params::n_neighbors=s;
        draw_from_scratch(type_attr[0], delay_rendering, create_transfer_function);
        draw_from_scratch(type_attr[1], delay_rendering, no_transfer_function);
        // type_attr[0].subgraph=filter_connectivity(*type_attr[0].graph);
        // type_attr[1].subgraph=filter_connectivity(*type_attr[1].graph);
        // redraw_points(false);
        // redraw_edges(true);
        draw();
    }
}

void mantle_vis_renderer::slot_min_spinbox(int m) {
    if (m!=params::min_size) {
        params::min_size=m;
        for (int i=0; i<2; ++i) {  
            if (type_attr[i].active)
                filter_polydata(type_attr[0]);
            redraw_points(type_attr[i], false);
            redraw_edges(type_attr[i], false);
        }
    }
    draw();
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
    else if (action->text() == "Save type +1 patches") {
        save_patches(type_attr[0]);
    }
    else if (action->text() == "Save type -1 patches") {
        save_patches(type_attr[1]);
    }
    else if (action->text() == "Hide control window") {
        show_hide_control(action);
    }
    else if (action->text() == "Show control window") {
        show_hide_control(action);
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
    bbox_t bounds=
        core_mantle::read_text(__vertices, params::filename, false);
    
    // filter out redundancy
    vertices_.clear();
    type_attr[0].positions.clear();
    type_attr[1].positions.clear();
    std::set<pos_t, less_pos_tol> unique_pos;
    integer_t n_points=0;
    for (integer_t i=0; i<__vertices.size(); ++i) {
        const vertex_t v=__vertices[i];
        const pos_t& x=v.position;
        if (unique_pos.find(x) == unique_pos.end()) {
            unique_pos.insert(x);
            vertices_.push_back(v);
            if (v.type>0) type_attr[0].positions.push_back(x);
            else type_attr[1].positions.push_back(x);
            ++n_points;
        } // else skip redundant vertex
    }
    
    if (params::verbose) {
        std::cout 
            << __vertices.size() << " vertices imported from file\n"
            << vertices_.size() << " unique vertices\n"
            << type_attr[0].positions.size() << " type +1\n"
            << type_attr[1].positions.size() << " type -1\n"
            << "bounding box:\n" << bounds << '\n';
    }
    
    return bounds;
}

void mantle_vis_renderer::camera_out() {
    QString fn=
        QFileDialog::getSaveFileName(this, 
                                     QString("Select an output camera file"),
                                     QString(home_directory.c_str()));
    std::string filename=fn.toStdString();
    vtk_utils::export_camera_settings(filename, renderer);
}

void mantle_vis_renderer::camera_in() {
    QString fn=
        QFileDialog::getOpenFileName(this,
                                     QString("Select an input camera file"),
                                     QString(home_directory.c_str()));
    std::string filename=fn.toStdString();
    vtk_utils::import_camera_settings(filename, renderer);
    draw();
}

void mantle_vis_renderer::load_file(const std::string& data_name) {
    import_vertices(data_name);
    internal_main();
}

void mantle_vis_renderer::save_snapshot() {
    QString fn=
        QFileDialog::getSaveFileName(this, QString("Select file name to save snapshot"),
                                     QString(home_directory.c_str()).
                                        append("/Desktop/untitled.png"));
    std::string filename=fn.toStdString();
    vtk_utils::save_frame(this->ui->qvtkWidget->GetRenderWindow(), 
                          filename);
}

void mantle_vis_renderer::save_patches(TypeAttributes& ta) {
    QString fn=
        QFileDialog::getSaveFileName(this, 
                                     QString("Select an output file"),
                                     QString(home_directory.c_str()));
    std::string filename=fn.toStdString();
    std::fstream file(filename.c_str(), std::ios::out);
    
    if (file) {
        // print out header
        file << "# Patches for " << params::filename 
             << ", neighborhood radius=" << params::search_radius
             << ", neighborhood max size=" << params::n_neighbors
             << ", patch min size=" << params::min_size
             << ", " << ta.type_name()
             << '\n';
        for (int i=0; i<ta.cc.size(); ++i) {
            const std::vector<integer_t>& cc=ta.cc[i];
            if (cc.size()<params::min_size) continue;
            file << "patch_size=" << cc.size() << '\n';
            std::copy(cc.begin(), cc.end(), 
            std::ostream_iterator<integer_t>(file, " "));
            file << '\n';
        }
        file.close();
    }
    else {
        std::cerr << "Operation failed. Unable to open " << filename << "\n";
    }
}


// **************************************************************************
// 
//                             Drawing
//
// **************************************************************************
void mantle_vis_renderer::redraw_cells(TypeAttributes& ta, bool _draw) {
    std::cout << "redraw_cells(" << ta.type_name() << ")\n" << std::flush;
    // if type is active: actor must be hidden, deleted, redrawn, and shown
    // if type is inactive: actor must be deleted and redrawn
    bool fb_changed=false;
    if (ta.active) {
        std::cout << ta.type_name() << " is active\n";
        std::cout << "hide cells for " << ta.type_name() << '\n';
        hide_cells(ta, delay_rendering);
        std::cout << "delete cells for " << ta.type_name() << '\n';
        for (int i=0; i<ta.meshes.size(); ++i) {
            ta.meshes[i].cells_actor=NULL;
        }
        std::cout << "show cells for " << ta.type_name() << '\n';
        show_cells(ta, delay_rendering);
        fb_changed=true;
    }
    else {
        std::cout << ta.type_name() << " is inactive\n";
        std::cout << "delete cells for " << ta.type_name() << '\n';
        for (int i=0; i<ta.meshes.size(); ++i) {
            ta.meshes[i].cells_actor=NULL;
        } // we will draw it once needed
    }
    
    if (fb_changed && _draw) draw();
}

void mantle_vis_renderer::redraw_edges(TypeAttributes& ta, bool _draw) {
    std::cout << "redraw_edges(" << ta.type_name() << ")\n" << std::flush;
    // if type is active: actor must be hidden, deleted, redrawn, and shown
    // if type is inactive: actor must be deleted and redrawn
    bool fb_changed=false;
    if (ta.active) {
        std::cout << ta.type_name() << " is active\n";
        std::cout << "hide edges for " << ta.type_name() << '\n';
        hide_edges(ta, delay_rendering);
        std::cout << "delete edges for " << ta.type_name() << '\n';
        for (int i=0; i<ta.meshes.size(); ++i) {
            ta.meshes[i].edges_actor=NULL;
        }
        std::cout << "show edges for " << ta.type_name() << '\n';
        show_edges(ta, delay_rendering);
        fb_changed=true;
    }
    else {
        std::cout << ta.type_name() << " is inactive\n";
        std::cout << "delete edges for " << ta.type_name() << '\n';
        for (int i=0; i<ta.meshes.size(); ++i) {
            ta.meshes[i].edges_actor=NULL;
        } // we will draw it once needed
    }
    
    if (fb_changed && _draw) draw();
}

void mantle_vis_renderer::redraw_points(TypeAttributes& ta, bool _draw) {
    std::cout << "redraw_points\n" << std::flush;
    // if type is active: points must be hidden, deleted, redrawn, and shown
    // if type is inactive: actor must be deleted and redrawn
    bool fb_changed=false;
    if (ta.active) {
        std::cout << ta.type_name() << " is active\n";
        hide_points(ta, delay_rendering);
        std::cout << "delete point set for " << ta.type_name() << '\n';
        for (int i=0; i<ta.meshes.size(); ++i) {
            ta.meshes[i].points_actor=NULL;
        }
        show_points(ta, delay_rendering);
        fb_changed=true;
    }
    else {
        std::cout << ta.type_name() << " is inactive\n";
        for (int i=0; i<ta.meshes.size(); ++i) {
            if (ta.meshes[i].points_actor!=NULL) {
                ta.meshes[i].points_actor=NULL;
                ta.meshes[i].points_actor=
                    get_points_actor(ta.meshes[i].geometry, ta.meshes[i].color);
            }
        }
    }
    
    if (fb_changed && _draw) draw();
}

void mantle_vis_renderer::draw() {
    if (this->initialized) {
        std::cout << "draw()\n";
        vtkActorCollection* actors=renderer->GetActors();/*
        std::cout << actors->GetNumberOfItems() 
            << " actors to render currently, of which "
            << renderer->VisibleActorCount() << " are visible\n";*/
        actors->InitTraversal();
        vtkActor* actor;
        int i=0;
        while ((actor=actors->GetNextActor())) {
            double bounds[6];
            actor->GetBounds(bounds);/*
            std::cout << "actor #" << i++ << " has bounds: "
                << "(" 
                << bounds[0] << ", " << bounds[2] << ", " << bounds[4]
                << ") -> (" 
                << bounds[1] << ", " << bounds[3] << ", " << bounds[5]
                << ")\n";*/
        }
        this->ui->qvtkWidget->GetRenderWindow()->Render();
    }
}

void mantle_vis_renderer::internal_main() {
    std::cout << "internal_main\n" << std::flush;
    
    import_vertices(params::filename);
    
    // create coordinate axes
    axes_widget=coordinate_axes_widget(this->ui->qvtkWidget);
    
    for (int i=0; i<2; ++i) {
        type_attr[i].active=params::do_type[i];
        type_attr[i].graph=std::shared_ptr<boost::graph_t>(new boost::graph_t());
    
        // determine superset of type-specific vertex connectivity
        boost::edge_map_t edge_weight;
        edge_weight=get(boost::edge_weight, *type_attr[i].graph);
        spurt::compute_connectivity<boost::graph_t, pos_t, 3>(*type_attr[i].graph, params::search_radius, type_attr[i].positions);
        draw_from_scratch(type_attr[i], delay_rendering, !i);
    }
    
    if (params::show_axes) {
        axes_widget->SetEnabled( 1 );
        axes_widget->InteractiveOn();
    }

//     std::cout << "calling draw()\n";
    draw();
    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
//     std::cout << "calling draw() again\n";
    draw();
}

void mantle_vis_renderer::
draw_from_scratch(TypeAttributes& ta, bool _draw, bool create_ctf) {
    std::cout << "draw_from_scratch(" << ta.type_name() << ")\n" << std::flush;
    
    std::cout << "\t1. there are currently " 
        << renderer->GetActors()->GetNumberOfItems() << " actors\n";
    
    hide_points(ta, false);
    hide_edges(ta, false);
    
    std::cout << "\t2. there are currently " 
        << renderer->GetActors()->GetNumberOfItems() << " actors\n";
        
    // identify subset of connectivity corresponding to chosen radius
    // and neighborhood size
    ta.subgraph=filter_connectivity(ta);

    if (params::verbose) {
        std::cout << "\t\t" << boost::num_edges(*ta.subgraph) 
                  << " edges selected for " << ta.type_name() << "\n";
    }

    if (false) {
        compute_vertex_degrees(ta.degrees, *ta.subgraph);
    }

    // initialize connected components
    std::vector<integer_t> cc_ids;
    spurt::connected_components(cc_ids, *ta.subgraph);
    ta.cc.clear();
    for (integer_t i=0; i<cc_ids.size(); ++i) {
        integer_t id=cc_ids[i];
        if (id>=ta.cc.size()) ta.cc.resize(id+1);
        ta.cc[id].push_back(i);
    }

    if (params::verbose) {
        std::cout << "\t\tThere are " << ta.cc.size() 
            << " connected components for " << ta.type_name() 
            << " before filtering\n";
    }

    // Create polydata to be visualized
    // ta.polydata=VTK_SMART(PolyData)::New();
    filter_polydata(ta);
    // ta.polydata=create_point_set(ta.polydata, ta.indices);

    // assign connected component id to each vertex
    /*ta.polydata=assign_component_ids(ta.polydata, ta.cc, create_ctf);*/
    // assign random color per connected component
    // ta.polydata=assign_value_by_component(ta.polydata, ta.cc);
    // filtered_polydata(ta.polydata);
    
    
    
    redraw_points(ta, false);
    redraw_edges(ta, false);
    redraw_cells(ta, false);

    if (_draw) {
        renderer->ResetCamera();
        renderer->ResetCameraClippingRange();
        draw();
    }
}


// **************************************************************************
// 
//                             Constructors/Destructors
//
// **************************************************************************
mantle_vis_renderer::mantle_vis_renderer(int argc, char* argv[]) : initialized(false) {
    std::cout << "constructor\n" << std::flush;
        
    if (!init(argc, argv)) {
        std::cout << "Invalid command line options. Exiting.\n";
        exit(1);
    }
    
    dummy=vtkPolyData::New();
    
    type_attr[0].type=1;
    type_attr[1].type=-1;
    
    srand48(time(NULL));
    
    progress=ProgressDisplay(params::show_progress);
    
    this->ui=new Ui_MainWindow;
    this->ui->setupUi(this);
    
    // use home directory as default location for file browsing
    home_directory=getpwuid(getuid())->pw_dir;
    
    // Initialize visualization pipeline and interactor
    renderer=VTK_SMART(Renderer)::New();
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(800, 800);
                               
    this->ui->qvtkWidget->SetRenderWindow(window);
    
    this->controls=new mantle_vis_control();
    this->controls->ui=new Ui_ControlWidget;
    this->controls->ui->setupUi(this->controls);
    this->controls->ui->pointsCheck->setTristate(false);
    this->controls->ui->edgesCheck->setTristate(false);
    this->controls->ui->cellsCheck->setTristate(false);
    this->controls->ui->axesCheck->setTristate(false);
    this->controls->ui->minusCheck->setTristate(false);
    this->controls->ui->plusCheck->setTristate(false);
    
    // create connections
    connect(this->controls->ui->pointsCheck, SIGNAL(stateChanged(int)),
            this, SLOT(slot_points_check(int)));
    connect(this->controls->ui->edgesCheck, SIGNAL(stateChanged(int)), 
            this, SLOT(slot_edges_check(int)));
    connect(this->controls->ui->cellsCheck, SIGNAL(stateChanged(int)), 
            this, SLOT(slot_cells_check(int)));
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
    connect(this->controls->ui->minSpin, SIGNAL(valueChanged(int)),
            this, SLOT(slot_min_spinbox(int)));
    connect(this->ui->menuFile, SIGNAL(triggered(QAction*)),
            this, SLOT(slot_file_action(QAction*)));
    connect(this->ui->menuWindow, SIGNAL(triggered(QAction*)),
            this, SLOT(slot_file_action(QAction*)));
            
    // initialize control values with parameters
    this->controls->ui->pointsCheck->setCheckState(params::show_points ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->edgesCheck->setCheckState(params::show_edges ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->cellsCheck->setCheckState(params::show_cells ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->axesCheck->setCheckState(params::show_axes ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->plusCheck->setCheckState(params::do_type[0] ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->minusCheck->setCheckState(params::do_type[1] ? Qt::Checked : Qt::Unchecked);
    this->controls->ui->radiusSpin->setValue(params::search_radius);
    this->controls->ui->radiusSlider->setValue(static_cast<int>(params::search_radius/params::max_search_radius*100.));
    this->controls->ui->sizeSpin->setValue(params::n_neighbors);
    this->controls->ui->minSpin->setValue(params::min_size);
        
    // configure rendering window
    renderer->SetBackground(params::bg_color[0],
                            params::bg_color[1],
                            params::bg_color[2]);
                            
    // create default color transfer function
    color_tf=create_color_transfer_function(std::vector<scalar_t>());
    
    this->initialized=true;
    internal_main();
    
    this->controls->show();
}

mantle_vis_control::mantle_vis_control() {
    this->ui=new Ui_ControlWidget;
    this->ui->setupUi(this);
}

mantle_vis_control::~mantle_vis_control() {
    delete this->ui;
}

mantle_vis_renderer::~mantle_vis_renderer() {
    delete this->controls;
    delete this->ui;
    dummy->Delete();
}