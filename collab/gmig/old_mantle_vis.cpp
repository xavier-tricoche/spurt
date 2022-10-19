#include <string>
#include <iostream>
#include <map>

#include "core_mantle.hpp"
#include "core_mantle_io.hpp"
#include "utils.hpp"

// xavier's utilities
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

using namespace xavier::gmig;

typedef float  scalar_t;
typedef int    integer_t;

const scalar_t _core_radius_  = 3486; // in km
const scalar_t _earth_radius_ = 6371; // in km

const scalar_t _infinity_ = std::numeric_limits<scalar_t>::max();

const integer_t invalid_index = static_cast<integer_t>(-1);

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

namespace boost {
typedef property<vertex_color_t, bool> vertex_prop_t;
typedef property<edge_weight_t, double> edge_prop_t;
typedef adjacency_list<vecS, vecS, bidirectionalS, vertex_prop_t, edge_prop_t> digraph_t;
typedef adjacency_list<vecS, vecS, undirectedS, vertex_prop_t, edge_prop_t> graph_t;
typedef graph_traits<graph_t> graph_traits_t;
typedef graph_traits_t::vertex_descriptor graph_vertex_t;
typedef graph_traits_t::edge_descriptor graph_edge_t;
typedef graph_traits<digraph_t> digraph_traits_t;
typedef digraph_traits_t::vertex_descriptor digraph_vertex_t;
typedef digraph_traits_t::edge_descriptor digraph_edge_t;
typedef graph_t::vertex_iterator graph_vertex_iterator;
typedef graph_t::edge_iterator graph_edge_iterator;
typedef digraph_t::vertex_iterator digraph_vertex_iterator;
typedef digraph_t::edge_iterator digraph_edge_iterator;
}

template<typename Key_>
class Counter : private std::map<Key_, size_t> 
{
    typedef typename std::map<Key_, size_t> base_t;
public:
    typedef Key_ key_t;
    typedef typename base_t::value_type value_t;
    typedef typename base_t::iterator iterator;
    typedef typename base_t::const_iterator const_iterator;
private:
    typedef std::pair<iterator, bool> answer_t;
public:
    Counter() : base_t(), __total_size(0), 
    __this(static_cast<base_t*>(this)),
    __cthis(static_cast<const base_t*>(this)) {}
    
    size_t increment(const key_t& key) {
        ++__total_size;
        answer_t what = __this->insert(value_t(key, 1));
        if (!what.second) ++(what.first->second);
        return what.first->second;
    }
    
    size_t decrement(const key_t& key) {
        iterator it = __this->find(key);
        if (it == __this->end() || 
            it->second==0) return 0;
        --__total_size;
        return --(it->second);
    }
    
    size_t operator[](const key_t& key) const {
        iterator it = __cthis->find(key);
        if (it != __cthis->end()) 
            return it->second;
        else return 0;
    }
    
    std::pair<size_t, size_t> size() const {
        return std::make_pair(__cthis->size(), __total_size);
    }
    
    iterator begin() { return __this->begin(); }
    iterator end() { return __this->end(); }
    const_iterator begin() const { return __cthis->begin(); }
    const_iterator end() const { return __cthis->end(); }
    
private:
    size_t __total_size;
    base_t* __this;
    const base_t* __cthis;
};

typedef Counter<integer_t> counter_t;

// Program parameters and options
std::string filename;
bool show_layers;
bool verbose = false;
bool world_coords = false;
bool show_earth = false;
bool show_axes = false;
bool show_lines = true;
std::string color_map_name;

enum color_mapping_t {
    color_by_sector,
    color_by_interface,
    color_by_curve,
    color_by_latitude,
    color_by_longitude,
    color_by_height,
    color_by_type
};
int color_map = color_mapping_t::color_by_curve;

typedef nvis::fixed_vector<scalar_t, 3> color_t;

struct VertexDistance {
    scalar_t operator()(integer_t vid1, integer_t vid2) const {
        const pos_t& p1 = vertices_[vid1].position;
        const pos_t& p2 = vertices_[vid2].position;
        // ignore longitudinal distance
        pos_t subv1(p1[0], 0, p1[2]);
        pos_t subv2(p2[0], 0, p2[2]);
        return nvis::norm(subv1-subv2); 
    }
};

bool init(int argc, char** argv) {
    namespace xcl = xavier::command_line;
    
    xcl::option_traits 
        required_group("Required parameters", true, false),
        positional_group("Positional parameters", true, true),
        optional_group("Optional parameters", false, false);
        
    xcl::option_parser parser(argv[0], "Visualizing core-mantle points");
    
    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", filename, "Input file (string)", 
                         positional_group);
        parser.add_flag("verbose", verbose, "Use verbose output", 
                        optional_group);
        parser.add_flag("world", world_coords, 
                        "Display data in spherical coordinates", 
                        optional_group);
        parser.add_flag("earth", show_earth, "Show earth crust", 
                        optional_group);
        parser.add_flag("axes", show_axes, "Show coordinates axes", 
                        optional_group);
        parser.add_value("colormap", color_map_name, "Color mapping input",
                         optional_group);
        parser.add_value("lines", show_lines, "Show layers as lines",
                         optional_group);
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

int main(int argc, char* argv[])
{
    if (!init(argc, argv)) {
        std::cout << "Invalid command line options. Exiting.\n";
        exit(1);
    }
    
    bbox_t bounds = core_mantle::read_text(vertices_, filename, false);
    std::cout << vertices_.size() << " vertices_ imported across region "
              << bounds << '\n';
    
    srand48(time(NULL));
            
    // ----------
    //
    // Algorithm:
    //
    // ----------
    // 1. sort vertices_ by longitude angle
    // 2. for each longitude angle
    //          determine curves connecting vertices_ of same layer_id
    // 3. connect curves across longitudes
    

    // Identify all the available sections, differentiating by longitude,
    // section indices, and type. Then determine which section each vertex 
    // belongs to and add corresponding reference to sections_
    for (int i=0 ; i<vertices_.size() ; ++i)
    {
        // update section data
        sections_.insert(section_id_t(vertices_[i]), i);
    }
    if (verbose) {
        std::cout << sections_.size() << " sections identified\n";
    }
    
    // Determine the section contents of each sector
    for (auto it=sections_.begin() ; it!=sections_.end() ; ++it) {
        const vertex_t& v = vertices_[it->second.front()];
        
        // add current section to its sector
        sectors_.insert(sector_id_t(v), it->first);
    }
    if (verbose) {
        std::cout << sectors_.size() << " sectors identified\n";
    }
    
    if (verbose) {
        std::cout << sections_.size() << " different sections available: \n";
        // for (auto it=sections_.begin(); it!=sections_.end() ; ++it) {
        //     std::cout << it->first << "\n";
        //     const section_t& s = it->second;
        //     if (section_id_t(vertices_[s[0]]) == section_id_t(165,2,2)) {
        //         std::cout << "Current section contains "
        //             << s.size() << " positions and those are \n";
        //         std::copy(s.begin(), s.end(),
        //                   std::ostream_iterator<integer_t>(std::cout, "\n"));
        //         std::cout << std::endl;
        //     }
        // }
    }

    int nb_used_verts = 0;
    
    // Determine individual layers in each Section
    for (auto it=sections_.begin() ; it!=sections_.end() ; ++it) {
        section_t& section = it->second;
        
        if (verbose) {
            std::cout << "Processing section #" 
                << std::distance(sections_.begin(), it) << "/"
                << sections_.size() << ", " 
                << layers_.size() << " layers so far\n";
        }
        
        if (!core_mantle::sanity_check(vertices_, section, verbose)) {
            std::cout << "\n\nERROR encountered in section " 
                << section_id_t(vertices_[section.front()]) << std::endl;
        }
                
        // Assign vertices of current section to correct layer
        for (int i=0 ; i<section.size() ; ++i) {
            integer_t v_id = section[i];
            vertex_t& v = vertices_[v_id];
            
            if (layers_.insert(layer_id_t(v), v_id)) {
                // add position of newly created layer to current section
                section.layers.push_back(layer_id_t(v));
            }
        }
    }
    
    std::for_each(layers_.begin(), layers_.end(), 
    [&](const layer_container::base_type::value_type& v) {
       if (v.second.size() < 2) {
           if (verbose) {
               std::cout << "layer " << v.first << " contains only "
                   << v.second.size() << " vertices\n";
           }
       } 
    });
    
    /*
    
        Determine layer correspondence between consecutive Sections
    
    */
    
    // Start by matching layers within same sector only
    for (auto it=sectors_.begin() ; it!=sectors_.end() ; ++it) {
        const sector_t& sector = it->second;
        
        if (verbose) {
            std::cout << "Processing sector #"
                << std::distance(sectors_.begin(), it) 
                << " of " << sectors_.size() << "\n";
        }
        
        // indices of all the interfaces that were found in the last section
        // (initially empty)
        std::vector<integer_t> active_fronts; 
        
        // for debugging purposes
        std::vector<std::pair<layer_id_t, layer_id_t> > error1_layers;
        std::vector<std::pair<layer_id_t, layer_id_t> > error2_layers;
        
        // iterate over sections contained in current sector in
        // increasing order of longitude
        if (verbose) {
            std::cout << "sector contains " << sector.size() << " sections\n";
        }
        for (int i=0 ; i<sector.size() ; ++i) {
            const section_t& section = sections_[sector[i]];
            
            // if (verbose) {
            //     std::cout << "\n\nThere are currently "
            //         << interfaces_.size() << " interfaces\n";
            //     for (int n=0 ; n<interfaces_.size() ; ++n) {
            //         std::cout << "-interface #" << n << " contains "
            //             << interfaces_[n].size()
            //             << " sections and its front contains "
            //             << layers_[interfaces_[n].back()].size() << " points\n";
            //     }
            // }
            
            if (active_fronts.empty()) {
                if (verbose) {
                    std::cout << "no active fronts in section #" << i << '\n';
                    std::cout << "creating " << section.layers.size() 
                        << " new interfaces\n";
                }
                // no active interface. each layer is treated
                // as front of a new interface
                for (int j=0 ; j<section.layers.size() ; ++j) {
                    interfaces_.push_back(interface_t());
                    interfaces_.back().push_back(section.layers[j]);
                    // the front just created is now active
                    active_fronts.push_back(interfaces_.size()-1);
                }
                continue; // we are done processing this section
            }
            
            if (verbose) {
                std::cout << "\n\tProcessing section #" << i
                          << "/" << sector.size()
                          << ", " << active_fronts.size()
                          << " active fronts\n";
                for (int n=0 ; n<active_fronts.size() ; ++n) {
                    std::cout << "front #" << n << " is last layer of active interface #"
                        << active_fronts[n] << " and contains "
                            << layers_[interfaces_[active_fronts[n]].back()].size()
                                << " points\n";
                    std::cout << "interface #" << active_fronts[n]
                        << " itself contains " << interfaces_[active_fronts[n]].size()
                            << " layers\n";
                }
            }
            
            // We construct a bipartite graph fronts - layers 
            // - each layer points to closest front
            // - each front points to closest layer
            size_t nb_fronts = active_fronts.size();
            size_t nb_layers = section.layers.size();
            if (verbose) {
                std::cout << "\t\tnb_fronts=" << nb_fronts
                    << ", nb layers=" << section.layers.size() << '\n';
            }
            // closest layers from each front
            std::vector<edge_t> front_to_layer
                (nb_fronts, edge_t(invalid_index, _infinity_));
            // closest fronts from each layer
            std::vector<edge_t> layer_to_front
                (nb_layers, edge_t(invalid_index, _infinity_));
            
            // create graph data structure 
            boost::digraph_t digraph(nb_fronts + nb_layers);
            boost::property_map<boost::digraph_t, boost::vertex_color_t>::type vcolor;
            boost::property_map<boost::digraph_t, boost::edge_weight_t>::type eweight;
            vcolor = get(boost::vertex_color, digraph);
            eweight = get(boost::edge_weight, digraph);
            
            for (int n=0 ; n<nb_layers ; ++n) {
                vcolor[n] = true;
            }
            for (int n=0 ; n<nb_fronts ; ++n) {
                vcolor[n+nb_layers] = false;
            }
            
            for (int j=0 ; j<section.layers.size() ; ++j) {
                const layer_t& layer = layers_[section.layers[j]];
                
                if (verbose) {
                    std::cout << "\t\tProcessing layer #" << j
                              << "/" << section.layers.size() << '\n';
                }
                
                // Find closest match for current layer among active fronts
                for (int k=0 ; k<nb_fronts ; ++k) {
                    const interface_t& interface = 
                        interfaces_[active_fronts[k]];
                    const layer_t& front = layers_[interface.back()];
                    
                    if (verbose) {
                        std::cout << "\t\t\tProcessing front #" << k
                            << "/" << nb_fronts << '\n';
                    }
                    
                    // check type compatibility between front and layer
                    if (interface.back().type != section.layers[j].type) {
                        if (verbose) {
                            std::cout << "incompatible types. skip.\n";
                            continue;
                        }
                    } 
                    
                    scalar_t _mean, _std_dev;
                    scalar_t dist_jk = 
                        hausdorff_distance<integer_t, VertexDistance, scalar_t>
                            (layer, front, _mean, _std_dev);
                    if (verbose) {
                        std:cout << "\t\t\t\tdist=" << dist_jk << ", mean=" 
                            << _mean << ", std dev=" << _std_dev << '\n';
                    }
                    if (dist_jk < front_to_layer[k].second) {
                        front_to_layer[k] = edge_t(j, dist_jk);
                    }
                    if (dist_jk < layer_to_front[j].second) {
                        layer_to_front[j] = edge_t(k, dist_jk);
                    }
                }
            }
            
            for (int j=0 ; j<nb_layers ; ++j) {
                if (layer_to_front[j].first != -1) {
                    boost::add_edge(j, layer_to_front[j].first+nb_layers, 
                                    layer_to_front[j].second, digraph);
                }
            }
            for (int k=0 ; k<nb_fronts ; ++k) {
                if (front_to_layer[k].first != -1) {
                    boost::add_edge(k+nb_layers, front_to_layer[k].first,
                                    front_to_layer[k].second, digraph);
                }
            }
            
            if (verbose) {
                std::cout << "Results for layer-to-front connections:\n";
                for (int _j=0 ; _j<nb_layers ; ++_j) {
                    edge_t edge = layer_to_front[_j];
                    if (edge.first == -1) {
                        if (verbose) {
                            std::cout << "for some reason, layer #" 
                                << _j << " has no valid front neighbor\n";
                        }
                        continue;
                    }                    
                    integer_t interf_id = active_fronts[edge.first];
                    const interface_t& interface = interfaces_[interf_id];
                    const layer_t& front = layers_[interface.back()];
                    const layer_id_t& layer_id = section.layers[_j];
                    const layer_t& _layer = layers_[layer_id];

                    integer_t _k = edge.first;
                    scalar_t dist = edge.second;
                    scalar_t mean, std_dev;
                    scalar_t dummy = 
                        hausdorff_distance<integer_t, VertexDistance, scalar_t>
                            (_layer, front, mean, std_dev);
                    std::cout << "front #" << _k 
                        << " (" << active_fronts[_k] << ", "
                        << front.size() << " points)"
                        << " is closest (" << dist 
                        << ") to layer #" << _j << " (" << section.layers[_j] 
                        << ", " << _layer.size() << " points)"
                        << " - mean distance = " << mean << ", std dev = "
                        << std_dev << '\n';
                }
                std::cout << "Results for front-to-layer connections:\n";
                for (int _k=0 ; _k<nb_fronts ; ++_k) {
                    edge_t edge = front_to_layer[_k];
                    if (edge.first == -1) {
                        if (verbose) {
                            std::cout << "for some other reason, front #"
                                << _k << " has no valid layer neighbor\n";
                        }
                        continue;
                    }
                    integer_t _j = edge.first;
                    const layer_id_t& layer_id = section.layers[_j];
                    const layer_t& _layer = layers_[layer_id];
                    scalar_t dist = edge.second;
                    integer_t interf_id = active_fronts[_k];
                    const interface_t& interface = interfaces_[interf_id];
                    const layer_t& front = layers_[interface.back()];
                    std::cout << "layer #" << _j
                        << " (" << section.layers[_j] << ", "
                        << _layer.size() << " points)"
                        << " is closest (" << dist 
                        << ") to front #" << _k << " (" << active_fronts[_k]
                        << ", " << front.size() << " points)\n";
                }
            }
            
            
            // Determine pairing between fronts and layers
            std::vector<integer_t> new_active_fronts;
            
            // Check incident edges of each layer / front: 0 incident edge
            // means no proximity to another layer
            boost::digraph_t::vertex_iterator vit, vend;
            boost::tie(vit, vend) = boost::vertices(digraph);
            for (; vit!=vend ; ++vit) {
                std::cout << "Vertex " << *vit << " ("
                    << (vcolor[*vit] ? "front" : "layer")
                    << ") has in-degree = "
                    << boost::in_degree(*vit, digraph)
                    << '\n';
                boost::graph_traits<boost::digraph_t>::in_edge_iterator eit, eend;
                boost::tie(eit, eend) = boost::in_edges(*vit, digraph);
                for (; eit!=eend ; ++eit) {
                    std::cout << "incident edge: " << boost::source(*eit, digraph)
                        << " --> " << boost::target(*eit, digraph) 
                        << " (w=" << eweight[*eit] << ")"<< '\n';
                }
            }
            
            typedef std::set<integer_t> set_t;
            std::vector<set_t> to_fronts(active_fronts.size());
            std::vector<set_t> to_layers(section.layers.size());
            
            for (int j=0 ; j<section.layers.size() ; ++j) {
                if (layer_to_front[j].first == -1) {
                    // if (verbose) {
                    //     std::cout << "layer_to_front[" << j << "] is invalid\n";
                    // }
                    continue;
                }
                to_fronts[layer_to_front[j].first].insert(j);
                add_edge(j, layer_to_front[j].first+nb_layers, 
                         layer_to_front[j].second, digraph);
            }
            
            if (verbose) {
                std::cout << "\t\tto_fronts created\n";
            }
            
            for (int k=0 ; k<active_fronts.size() ; ++k) {
                // if (verbose) {
                //     std::cout << "\t\t\tfront_to_layer.size()="
                //         << front_to_layer.size() << '\n';
                //     std::cout << "\t\t\tfront_to_layer[" << k << "].first="
                //               << front_to_layer[k].first
                //               << "/" << to_layers.size() << '\n';
                // }
                if (front_to_layer[k].first == -1) {
                    // if (verbose) {
                    //     std::cout << "front_to_layer[" << k << "] is invalid\n";
                    // }
                    continue;
                }
                to_layers[front_to_layer[k].first].insert(k);
                add_edge(nb_layers+k, front_to_layer[k].first,
                         front_to_layer[k].second, digraph);
            }
                        
            if (verbose) {
                std::cout << "\t\tto_layers created\n";
            }
            
            for (int j=0 ; j<section.layers.size() ; ++j) {
                
                if (verbose) {
                    std::cout << "\t\t\tProcessing layer #"
                        << j << "/" << section.layers.size()
                            << '\n';
                }
                
                int layer_degree = to_layers[j].size();
                
                // if (verbose) {
                //     std::cout << "\t\t\t\tlayer degree="
                //         << layer_degree << '\n';
                // }
                
                if (layer_degree == 0) {
                    // Case 5: new front
                    interfaces_.push_back(interface_t());
                    interfaces_.back().push_back(section.layers[j]);
                    new_active_fronts.push_back(interfaces_.size()-1);
                    
                    if (verbose) {
                        std::cout << "this layer was not selected by any "
                            << " front -> new interface\n";
                    }
                }
                else if (layer_degree == 1) {
                    integer_t which_front = *to_layers[j].begin();
                    int front_degree = to_fronts[which_front].size();
                    if (front_degree == 1) {
                        integer_t which_layer =
                            *to_fronts[which_front].begin();
                        if (which_layer == j) {
                            // Case 1
                            integer_t front_i = active_fronts[which_front];
                            layer_id_t layer_id = section.layers[j];
                            interfaces_[front_i].push_back(layer_id);
                            // current front remains active
                            new_active_fronts.push_back(front_i);
                            if (verbose) {
                                std::cout << "expanding interface #" << front_i
                                    << '\n';
                            }
                        }
                        else {
                            // There is a closer layer for this layer's closest 
                            // front. Case 5: new front
                            interfaces_.push_back(interface_t());
                            interfaces_.back().push_back(section.layers[j]);
                            new_active_fronts.push_back(
                                interfaces_.size()-1);
                            if (verbose) {
                                std::cout << "creating a new interface\n";
                            }
                        }
                    }
                    else { 
                        // the unique front for which the current layer
                        // is closest is itself the closest neighbor of one or 
                        // more other layers.
                        std::cout << "Case ? not handled currently\n";
                        error1_layers.push_back(
                            std::make_pair(section.layers[j], 
                                           interfaces_[
                                               active_fronts[which_front]].
                                                   back())
                            );
                        std::cout << "error case. doing nothing\n";
                    }
                }
                else {
                    // current layer is closest neighbor of multiple fronts
                    // need to confirm that this layer's closest front is 
                    // among those to establish Case 2.
                    if (to_layers[j].find(layer_to_front[j].first) != 
                        to_layers[j].end()) {
                        // merge fronts
                        set_t::iterator it = to_layers[j].begin();
                        integer_t front1_i = active_fronts[*it];
                        interface_t& itf1 = interfaces_[front1_i];
                        itf1.push_back(section.layers[j]);
                        for (++it ; it!=to_layers[j].end() ; ++it) {
                            integer_t front_i = active_fronts[*it];
                            interface_t& itf = interfaces_[front_i];
                            std::copy(itf.begin(), itf.end(),
                                      std::back_inserter(itf1));
                            std::cout << "merging interface #"
                                << front_i << " into interface #"
                                << front1_i << '\n';
                            std::cout << "\tinterface #" << front_i 
                                << " is now empty\n";
                            itf.clear();
                        }
                        new_active_fronts.push_back(front1_i);
                        if (verbose) {
                            std::cout << "merging multiple layers into "
                                << "interface #" << front1_i << '\n';
                        }
                    }
                    else {
                        // funky case...
                        std::cout << "Funky case not handled currently\n";
                        error2_layers.push_back(
                            std::make_pair(section.layers[j],
                                           interfaces_[
                                                  active_fronts[
                                                      layer_to_front[j].first]].
                                                          back())
                            );
                        std::cout << "error case. doing nothing\n";
                    }
                }
            }
        }
    }
    
    counter_t counter;
    std::for_each(interfaces_.begin(), interfaces_.end(), 
                  [&](const interface_t& interface)
        {
            counter.increment(interface.size());
        });
        
    counter_t counter_layers;
    std::for_each(layers_.begin(), layers_.end(), 
                  [&](const std::pair<layer_id_t,layer_t>& layer)
        {
            counter_layers.increment(layer.second.size());
        });
        
    counter_t counter_sections;
    std::for_each(sections_.begin(), sections_.end(), 
                  [&](const std::pair<section_id_t, section_t>& section)
        {
            counter_sections.increment(section.second.size());
        });
    
    // Display individual interfaces
    std::vector<pos_t> positions;
    positions.reserve(vertices_.size());
    std::for_each(vertices_.begin(), vertices_.end(), 
                  [&](const vertex_t& v)
    {
        if (!world_coords) {
            positions.push_back(v.position);
            // rescaling to improve aspect ratio
            positions.back()[1] *= 2;
            positions.back()[2] /= 3;
        }
        else {
            scalar_t r = _core_radius_ + v.position[2];
            scalar_t la = deg2rad(v.latitude());
            scalar_t lo = deg2rad(v.longitude());
            pos_t p(r*sin(la)*cos(lo),
                    r*sin(la)*sin(lo),
                    r*cos(la));
            positions.push_back(p);
        }
    });

    std::vector<float> values(vertices_.size(), 0);
    
    std::cout << "Summary:\n"
              << '\t' << vertices_.size() << " vertices\n"
              << '\t' << sectors_.size() << " sectors\n"
              << '\t' << sections_.size() << " sections\n"
              << '\t' << layers_.size() << " layers\n"
              << '\t' << interfaces_.size() << " interfaces\n";

    std::cout << "\nInterface statistics:\n";
    std::for_each(counter.begin(), counter.end(), 
        [&](counter_t::value_t& v){
            std::cout << v.first << ": " << v.second 
                      << " (" << 100.*(float)v.second/interfaces_.size() 
                      << "\%)\n";
        });
    // for (int i=0 ; i<interfaces_.size() ; ++i) {
    //     std::cout << "interface #" << i
    //               << " contains " << interfaces_[i].size()
    //               << " layers\n";
    // }
    std::cout << "\nLayer statistics:\n";
    std::for_each(counter_layers.begin(), counter_layers.end(), 
        [&](counter_t::value_t& v){
            std::cout << v.first << ": " << v.second 
                      << " (" << 100.*(float)v.second/layers_.size() 
                      << "\%)\n";
        });
    std::cout << "\nSection statistics:\n";
    std::for_each(counter_sections.begin(), counter_sections.end(), 
        [&](counter_t::value_t& v){
            std::cout << v.first << ": " << v.second 
                      << " (" << 100.*(float)v.second/sections_.size() 
                      << "\%)\n";
        });
    
    if (color_map_name == "color_by_sector") {
        color_map = color_mapping_t::color_by_sector;
    }
    else if (color_map_name == "color_by_interface") {
        color_map = color_mapping_t::color_by_interface;
    }
    else if (color_map_name == "color_by_curve") {
        color_map = color_mapping_t::color_by_curve;
    }
    else if (color_map_name == "color_by_latitude") {
        color_map = color_mapping_t::color_by_latitude;
    }
    else if (color_map_name == "color_by_longitude") {
        color_map = color_mapping_t::color_by_longitude;
    }
    else if (color_map_name == "color_by_height") {
        color_map = color_mapping_t::color_by_height;
    }
    else if (color_map_name == "color_by_type") {
        color_map = color_mapping_t::color_by_type;
    }
    
    std::vector<bool> vertex_mask(vertices_.size(), true);
    std::map<int, int> interface_sizes;
    
    switch(color_map) {
        case color_mapping_t::color_by_sector: {
            // assign sector index as vertex value
            int i=0;
            for (auto it=sectors_.begin() ; it!=sectors_.end() ; ++it, ++i) {
                const sector_t& sector = it->second;
                for (int j=0 ; j<sector.size() ; ++j) {
                    const section_t& section = 
                        sections_[sector[j]];
                    for (int k=0 ; k<section.size() ; ++k) {
                        values[section[k]] = i;
                    }
                }
            }
            break;
        }
        case color_mapping_t::color_by_interface: {
            scalar_t min_longitude, max_longitude;
            min_longitude = std::numeric_limits<scalar_t>::max();
            max_longitude = std::numeric_limits<scalar_t>::min();
            // assign interface index as vertex value
            std::fill(vertex_mask.begin(), vertex_mask.end(), false);
            for (int i=0 ; i<interfaces_.size() ; ++i) {
                // if (/*i && */i!=2629) continue;
                const interface_t& itf = interfaces_[i];
                if (itf.empty()) continue;
                
                if (verbose) {
                    std::cout << "Interface #" << i 
                              << " contains "
                              << itf.size() << " layers\n";
                }
                
                scalar_t loc_min_long, loc_max_long;
                loc_min_long = itf[0].longitude;
                loc_max_long = itf.back().longitude;
                if (loc_min_long < min_longitude) 
                    min_longitude = loc_min_long;
                if (loc_max_long > max_longitude)
                    max_longitude = loc_max_long;
                
                scalar_t val = drand48();
                std::cout << "Assigning value " << val << " to interface " 
                    << i << '\n';
                for (int j=0 ; j<itf.size() ; ++j) {
                    const layer_t& layer = layers_[itf[j]];
                    for (int k=0 ; k<layer.size() ; ++k) {
                        vertex_mask[layer[k]] = true;
                        values[layer[k]] = val;
                    }
                }
            }
            /*
            for (int i=0 ; i<vertices_.size() ; ++i) {
                if (vertices_[i].longitude() <= max_longitude && 
                    vertices_[i].longitude() >= min_longitude &&
                    values[i] == 0) {
                    vertex_mask[i] = true;
                }
            }
            */
            break;
        }
        case color_mapping_t::color_by_curve: {
            // assign layer index as vertex value
            for (auto it=layers_.begin() ; it!=layers_.end() ; ++it) {
                const layer_t& layer = it->second;
                scalar_t val = drand48();
                for (int j=0 ; j<layer.size() ; ++j) {
                    values[layer[j]] = val;
                }
            }
            break;
        }
        case color_mapping_t::color_by_latitude: {
            // assign interface index as vertex value
            for (int i=0 ; i<vertices_.size() ; ++i) {
               values[i] = vertices_[i].latitude();
            }
            break;
        }
        case color_mapping_t::color_by_longitude: {
            // assign interface index as vertex value
            for (int i=0 ; i<vertices_.size() ; ++i) {
               values[i] = vertices_[i].longitude();
            }
            break;
        }
        case color_mapping_t::color_by_height: {
            // assign interface index as vertex value
            for (int i=0 ; i<vertices_.size() ; ++i) {
               values[i] = vertices_[i].height();
            }
            break;
        }
        case color_mapping_t::color_by_type: {
            if (verbose) {
                std::cout << "coloring by type\n";
            }
            // assign interface index as vertex value
            for (int i=0 ; i<vertices_.size() ; ++i) {
               values[i] = vertices_[i].type;
            }
            break;
        }
        default: {
            // uniform color
        }

    }
    
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
    
    if (show_lines) {
        VTK_CREATE(vtkCellArray, curves);
        for (auto it=sections_.begin() ; it!=sections_.end() ; ++it) {
            section_id_t sid = it->first;
            const section_t& sec = it->second;
            // if (sid.type == -1) continue; // skip type -1 sections
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
                const layer_t& layer = layers_[sec.layers[l]];        
                if (layer.empty() || !vertex_mask[layer.front()]) continue;
                curves->InsertNextCell(layer.size());
                std::for_each(layer.begin(), layer.end(), [&](const integer_t& id)
                {
                    curves->InsertCellPoint(id);
                });
            }
        }
        pd->SetLines(curves);
    }
    
    // create color map
    VTK_CREATE(vtkColorTransferFunction, ctf);
    
    if (color_map == color_mapping_t::color_by_type) {
        ctf->AddRGBPoint(-1, 1, 1, 0);
        ctf->AddRGBPoint(1, 0, 0, 1);
    }
    else {
        // adaptive color map with spiral color scale
        std::vector<color_t> colors(20);
        xavier::spiral_scale(colors, 20, 0.2);
        scalar_t minval = *std::min_element(values.begin(), values.end());
        scalar_t maxval = *std::max_element(values.begin(), values.end());
        scalar_t dval = (maxval-minval)/19;
        
        for (int i=0 ; i<20 ; ++i) {
            scalar_t value = minval + i*dval;
            color_t color = colors[i];
            ctf->AddRGBPoint(value, color[0], color[1], color[2]);
        }
    }
    
    if (show_lines) {
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
    }
    else {
        vtkSmartPointer<vtkPolyData> spheres = 
            vtk_utils::make_spheres(pd, 1);
        VTK_MAKE_ACTOR(sphere_actor, spheres);
        sphere_actor->GetMapper()->SetLookupTable(ctf);
        sphere_actor->GetMapper()->ScalarVisibilityOn();
        renderer->AddActor(sphere_actor);
    }
    
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
        color_t c = xavier::rainbow[1+((2*i)%15)];
        actor->GetProperty()->SetColor(c[0], c[1], c[2]);
        renderer->AddActor(actor);
    }
#endif
    
    if (show_earth) {
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
    renderer->SetBackground(1,1,1);
    
    vtkSmartPointer<vtkAxesActor> axes = 
        vtkSmartPointer<vtkAxesActor>::New();
    
    double __cyl_radius = axes->GetCylinderRadius();
    double __cone_radius = axes->GetConeRadius();
    int __cyl_res = axes->GetCylinderResolution();
    std::cout << "cylinder radius = " << __cyl_radius << '\n';
    std::cout << "cone radius = " << __cone_radius << '\n';
    std::cout << "cylinder resolution = " << __cyl_res << '\n';
    
    axes->SetCylinderRadius(0.05);
    axes->SetShaftTypeToCylinder();
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,0);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,0);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,0);
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetBackgroundColor(0.75,0.75,0.75);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetBackgroundColor(0.75,0.75,0.75);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetBackgroundColor(0.75,0.75,0.75);
 
    vtkSmartPointer<vtkOrientationMarkerWidget> widget = 
      vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor( 1, 1, 1 );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( interactor );
    widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
    if (show_axes) {
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