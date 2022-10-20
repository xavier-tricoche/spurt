#include <string>
#include <iostream>
#include <map>

// Boost
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include "core_mantle.hpp"
#include "core_mantle_io.hpp"
#include "utils.hpp"

// spurt's utilities
#include <format/format.hpp>
#include <misc/option_parse.hpp>

// VTK
#include <vtkColorTransferFunction.h>
#include <vtkDelaunay2D.h>
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

using namespace spurt::gmig;

const float core_radius = 3486; // in km
const float earth_radius = 6371; // in km

typedef float  scalar_t;
typedef int    index_t;
typedef int    index_t;

typedef core_mantle::Vertex<scalar_t,index_t> vertex_t;
typedef vertex_t::pos_t    pos_t;
typedef vertex_t::bbox_t   bbox_t;
typedef vertex_t::index_t  index_t;

using boost::multi_index_container;
using namespace boost::multi_index;

// ---------------------------------------------------------------
// Vertex container supporting efficient search according to all
// relevant pieces of information associated with a Vertex object:
// ---------------------------------------------------------------
// - longitude (member function)
// - latitude (member function)
// - height (member function)
// - lat_section (member)
// - long_section (member)
// - section_index (member function)

// phony types to be used as tags
struct longitude{};
struct latitude{};
struct height{};
struct latitude_section{};
struct longitude_section{};
struct section_index{};
typedef multi_index_container<
    vertex_t,
    indexed_by<
        ordered_non_unique<
            tag<longitude>,
            BOOST_MULTI_INDEX_CONST_MEM_FUN(vertex_t, scalar_t, longitude)
        >,
        ordered_non_unique<
            tag<latitude>, 
            BOOST_MULTI_INDEX_CONST_MEM_FUN(vertex_t, scalar_t, latitude)
        >,
        ordered_non_unique<
            tag<height>, 
            BOOST_MULTI_INDEX_CONST_MEM_FUN(vertex_t, scalar_t, height)
        >,
        ordered_non_unique<
            tag<latitude_section>, 
            BOOST_MULTI_INDEX_MEMBER(vertex_t, index_t, lat_section)
        >,
        ordered_non_unique<
            tag<longitude_section>, 
            BOOST_MULTI_INDEX_MEMBER(vertex_t, index_t, long_section)
        >,
        ordered_non_unique<
            tag<section_index>, 
            BOOST_MULTI_INDEX_CONST_MEM_FUN(vertex_t, index_t, section_index)
        >
    >
> vertex_set_t;
            
typedef boost::multi_index::index<vertex_set_t, longitude>::type
    longitude_map; 
typedef boost::multi_index::index<vertex_set_t, latitude>::type
    latitude_map;
typedef boost::multi_index::index<vertex_set_t, height>::type
    height_map;
typedef boost::multi_index::index<vertex_set_t, latitude_section>::type
    latitude_section_map; 
typedef boost::multi_index::index<vertex_set_t, longitude_section>::type
    longitude_section_map; 
typedef boost::multi_index::index<vertex_set_t, section_index>::type
    sections_map; 

// typedef core_mantle::Section<index_t>             section_t;
typedef std::vector<index_t>                      layer_t;
// typedef std::vector<index_t>                      interface_t;
// typedef core_mantle::SectionID<scalar_t,index_t>  section_id_t;
// typedef core_mantle::LayerID<scalar_t,index_t>    layer_id_t;

// typedef std::map<section_id_t, index_t> section_map_t;
// typedef std::map<layer_id_t, index_t>   layer_map_t;
typedef longitude_map section_map_t;
typedef section_map_t::iterator         smap_iterator;
// typedef layer_map_t::iterator           lmap_iterator;

int main(int argc, char* argv[])
{
    namespace xcl = spurt::command_line;
    
    std::string filename;
    bool show_layers;
    bool verbose = false;
    bool world_coords = false;
    bool show_earth = false;
    
    xcl::option_traits 
        required_group("Required parameters", true, false),
        positional_group("Positional parameters", true, true),
        optional_group("Optional parameters", false, false);
        
    xcl::option_parser parser(argv[0], "Visualizing core-mantle points");
    
    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", filename, "Input file (string)", positional_group);
        parser.add_flag("verbose", verbose, "Use verbose output", optional_group);
        parser.add_flag("world", world_coords, "Display data in spherical coordinates", optional_group);
        parser.add_flag("earth", show_earth, "Show earth crust", optional_group);
        
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
                  << "Command line options enteredso far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    // ----------------------
    //
    // Central data structure
    //
    // ----------------------
    
    // data storage
    std::vector<vertex_t>   vertices;
    // std::vector<section_t>  section_data;
    std::vector<layer_t>    layer_data;
    
    // // data indexing
    // section_map_t           section_map;
    // layer_map_t             layer_map;
    vertex_set_t vertex_container;
    bbox_t bounds = core_mantle::read_text(vertices, filename, false);
    
    std::cout << "imported " << vertices.size() << " vertices across region "
        << bounds << '\n';
    
    for (int i=0 ; i<vertices.size() ; ++i) {
        vertex_container.insert(vertices[i]);
    }
    
    // test my understanding of boost::multi_index_container's API
    if (true) {
        // test 1. List all distinct longitude values 
        longitude_map& lo_map = get<longitude>(vertex_container);
        typedef longitude_map::iterator lo_iterator;
        std::cout << "Available longitudes:\n";
        for (lo_iterator it=lo_map.begin() ; it!=lo_map.end() ;) {
            scalar_t lo = it->longitude();
            std::cout << lo;
            std::pair<lo_iterator, lo_iterator> 
                cur_range = lo_map.equal_range(lo);
            const lo_iterator _begin = cur_range.first;
            const lo_iterator _end = cur_range.second;
            std::cout << ": " 
                      << std::distance(_begin, _end)
                      << " vertices sharing this longitude\n";
            it = cur_range.second;
            
            std::set<scalar_t> unique_latitudes;
            std::set<index_t> unique_sections;
            for (lo_iterator jt=_begin ; jt!=_end ; ++jt) {
                unique_latitudes.insert(jt->latitude());
                unique_sections.insert(jt->section_index());
            }
            std::cout << "\t" << unique_latitudes.size()
                      << " different latitudes at that longitude ranging from " 
                      << *unique_latitudes.begin() << " to "
                      << *unique_latitudes.rbegin() << '\n';
            std::cout << "\t" << unique_sections.size()
                      << " different sections at that longitude ranging from "
                      << *unique_sections.begin() << " to "
                      << *unique_sections.rbegin() << '\n';
        }
    }
    
#if 0
    
    // ----------
    //
    // Algorithm:
    //
    // ----------
    // 1. sort vertices by longitude angle
    // 2. for each longitude angle
    //          determine curves connecting vertices of same layer_id
    // 3. connect curves across longitudes
    

    // Identify all the available sections, differentiating by longitude,
    // section indices, and type. Then determine which section each vertex 
    // belongs to and add corresponding reference to section_data
    
    for (int i=0 ; i<vertices.size() ; ++i)
    {
        section_id_t id(vertices[i]);
        // location of a tentative new section in central repository
        int nb_sections = section_data.size(); 
        std::pair<smap_iterator, bool> where =
            section_map.insert(std::make_pair(id, nb_sections));
        if (!where.second) { // insertion failed: Section exists already
            // add current vertex id to this section
            section_data[where.first->second].push_back(i);
        }
        else {
            // create new section entry in central repository
            section_data.push_back(section_t());
            section_data.back().push_back(i);
        }
    }
    
    if (verbose) {
        std::cout << section_data.size() << " different sections available: \n";
        for (auto it=section_map.begin(); it!=section_map.end() ; ++it) {
            std::cout << it->first << "\n";
            const section_t& s = section_data[it->second];
            if (section_id_t(vertices[s[0]]) == section_id_t(165,2,2)) {
                std::cout << "Current section contains " 
                    << s.size() << " positions and those are \n";
                std::copy(s.begin(), s.end(), 
                          std::ostream_iterator<index_t>(std::cout, "\n"));
                std::cout << std::endl;
            }
        }
    }

    int nb_used_verts = 0;
    
    // Determine individual layers in each Section
    for (smap_iterator it=section_map.begin() ; it!=section_map.end() ; ++it) {
        section_t& sec = section_data[it->second];
        
        if (!core_mantle::sanity_check(vertices, sec)) {
            std::cout << "\n\nERROR encountered in section " 
                << section_id_t(vertices[sec.front()]) << std::endl;
        }
        
        int nlayers=0;
        
        // Assign vertices of current section to correct layer
        for (int i=0 ; i<sec.size() ; ++i) {
            index_t v_id = sec[i];
            vertex_t& v = vertices[v_id];
            layer_id_t l_id(v);
            // location of a tentative new layer in central repository
            int nb_layers = layer_data.size(); 
            // check if a new layer needs to be created
            std::pair<lmap_iterator, bool> where =
                layer_map.insert(std::make_pair(l_id, nb_layers));
            if (!where.second) { // layer already encountered
                // add current vertex id to existing layer
                layer_data[where.first->second].push_back(v_id);
            }
            else {
                // create new layer entry in central repository
                layer_data.push_back(layer_t());
                layer_data.back().push_back(v_id);
                // add reference to newly created layer in current section
                sec.layers.push_back(nb_layers);
                ++nlayers;
            }
        }
    }
    
    // Determine layer correspondence between consecutive Sections
    typedef std::map<layer_id_t, int> layer_matches;
    std::vector<interface_t> interfaces;
    
    // iterate over all sections
    for (smap_iterator it=section_map.begin(); it!=section_map.end(); ++it) {
        const section_t sec = section_data[it->second];
        section_id_t sec_id(vertices[sec[0]]);
        
        // for each layer, determine potential match with previously
        // processed layers / section
        for (int i=0 ; i<sec.layers.size() ; ++i) {
            const layer_t& cur_layer = layer_data[sec.layers[i]];
            layer_id_t lay_id(vertices[cur_layer[0]]);
        }
    }
     
    
    // Display individual layers in their respective section
    std::vector<pos_t> positions;
    positions.reserve(vertices.size());
    std::vector<float> values;
    values.reserve(vertices.size());
    std::for_each(vertices.begin(), vertices.end(), [&](const vertex_t& v)
    {
        if (!world_coords) {
            positions.push_back(v.position);
            // rescaling to improve aspect ratio
            positions.back()[1] *= 2;
            positions.back()[2] /= 3;
        }
        else {
            scalar_t r = core_radius + v.position[2];
            scalar_t la = deg2rad(v.latitude());
            scalar_t lo = deg2rad(v.longitude());
            pos_t p(r*sin(la)*cos(lo),
                    r*sin(la)*sin(lo),
                    r*cos(la));
            positions.push_back(p);
        }
        
        values.push_back(v.layer_id);
    });
    
    typedef nvis::fvec3 color_t;
    
    // create rendering engine and render window
    VTK_CREATE(vtkRenderer, renderer);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(800, 800);
    
    vtkSmartPointer<vtkPolyData> pd = vtk_utils::make_points(positions);
    if (verbose) {
        pd->ComputeBounds();
        double _bounds[6];
        pd->GetBounds(_bounds);
        std::cout << "bounds of created polydata are: \n"
                  << _bounds[0] << " -> " << _bounds[1] << "\n"
                  << _bounds[2] << " -> " << _bounds[3] << "\n"
                  << _bounds[4] << " -> " << _bounds[5] << '\n';
    }
    vtk_utils::add_scalars(pd, values);
    
    VTK_CREATE(vtkCellArray, curves);
    for (int i=0 ; i<section_data.size() ; ++i) {
        const section_t& sec = section_data[i];
        section_id_t sid(vertices[sec.front()]);
        if (sid.t == -1) continue; // skip type -1 sections
#if 0
        VTK_CREATE(vtkPlaneSource, psource);
        psource->SetCenter(0, sec.longitude(), 0);
        psource->SetNormal(0, 1, 0);
        psource->SetPoint1(100, sec.longitude(), 0);
        psource->SetPoint2(0, sec.longitude(), 100);
        psource->Update();
        vtkSmartPointer<vtkPolyData> plane(psource->GetOutput());
        VTK_MAKE_ACTOR(plane_actor, plane);
        plane_actor->GetProperty()->SetColor(0, 0, 0.5);
        renderer->AddActor(plane_actor);
#endif
        for (int l=0 ; l<sec.layers.size() ; ++l) {
            const layer_t& layer = layer_data[sec.layers[l]];
            
            scalar_t ref_long = vertices[layer.front()].longitude();
            
            curves->InsertNextCell(layer.size());
            std::for_each(layer.begin(), layer.end(), [&](const index_t& id)
            {
                curves->InsertCellPoint(id);
                if (vertices[id].longitude() != ref_long) {
                    std::cout << "ERROR in section " << sid 
                              << ": ref longitude is: " << ref_long
                              << " but position " << id 
                              << " within same layer has longitude "
                              << vertices[id].longitude() << '\n';
                    std::cout << "For reference: first position is " 
                              << vertices[layer.front()]
                              << " and current position is " 
                              << vertices[id] << '\n'
                              << "This section contains " 
                              << sec.layers.size() 
                              << " layers and current layer has size "
                              << layer.size() << std::endl;
                }
            });
        }
    }
    pd->SetLines(curves);
    
    VTK_CREATE(vtkTubeFilter, tubef);
    VTK_CONNECT(tubef, pd);
    tubef->SetRadius(1);
    tubef->SetNumberOfSides(8);
    tubef->Update();
    VTK_MAKE_ACTOR(tube_actor, tubef->GetOutput());
    
    VTK_CREATE(vtkColorTransferFunction, ctf);
    ctf->AddRGBPoint(0, 0.9, 0.9, 0.9);
    ctf->AddRGBPoint(1, 0, 0, 1);
    ctf->AddRGBPoint(2, 0, 0.5, 1);
    ctf->AddRGBPoint(3, 0, 1, 1);
    ctf->AddRGBPoint(4, 0, 1, 0.5);
    ctf->AddRGBPoint(5, 0, 1, 0);
    ctf->AddRGBPoint(6, 0.5, 1, 0);
    ctf->AddRGBPoint(7, 1, 1, 0);
    ctf->AddRGBPoint(8, 1, 0.5, 0);
    ctf->AddRGBPoint(9, 1, 0, 0);
    ctf->AddRGBPoint(10, 1, 0, 0.5);
    
    // tube_actor->GetProperty()->SetColor(1,0,0);
    tube_actor->GetMapper()->SetLookupTable(ctf);
    tube_actor->GetMapper()->ScalarVisibilityOn();
    renderer->AddActor(tube_actor);
    
#if 0
    // for each layer
    for (int i=1 ; i<2; ++i) {
        // create a polydata corresponding to the points included in the layer
        vtkPolyData* pd = vtk_utils::make_points(layer_pos[i]);
        // scale up x and y coordinates
        VTK_CREATE(vtkTransform, transform);
        transform->Scale(10, 10, 1);
        VTK_CREATE(vtkTransformFilter, filter);
        vtk_connect(filter, pd);
        filter->SetTransform(transform);
        // triangulate
        VTK_CREATE(vtkDelaunay2D, del2D);
        VTK_PLUG(del2D, filter);
        del2D->Update();
        // visualize
        VTK_MAKE_ACTOR(actor,del2D->GetOutput());
        actor->GetMapper()->ScalarVisibilityOff();
        color_t c = spurt::rainbow[1+((2*i)%15)];
        actor->GetProperty()->SetColor(c[0], c[1], c[2]);
        renderer->AddActor(actor);
    }
#endif
    
    if (show_earth) {
        std::cout << "showing earth crust\n";
        VTK_CREATE(vtkSphereSource, earth);
        earth->SetRadius(earth_radius);
        earth->SetCenter(0,0,0);
        earth->SetThetaResolution(50);
        earth->SetPhiResolution(50);
        earth->Update();
        VTK_MAKE_ACTOR(earth_actor, earth->GetOutput());
        earth_actor->GetProperty()->SetColor(0,0,0);
        earth_actor->GetProperty()->SetOpacity(0.1);
        renderer->AddActor(earth_actor);
    }
    
    // configure rendering window
    renderer->SetBackground(1,1,1);
    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    
    // enter interaction loop
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    window->Render();
    interactor->Initialize();
    interactor->Start();
#endif 
       
    return 0;
}