#ifndef __GRAPH_DEFINITIONS_HPP__
#define __GRAPH_DEFINITIONS_HPP__

// Boost graph library
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/utility.hpp>
#include <boost/graph/copy.hpp>

namespace spurt {

template <typename Scalar_, typename PropertyMap_>
struct property_filter {
    typedef Scalar_ scalar_t;
    typedef PropertyMap_ property_map_t; 
    
    property_filter(scalar_t threshold=0) {}
    property_filter(property_map_t pmap, scalar_t threshold=0)
        : m_pmap(pmap), m_threshold(threshold) {}
                  
    template <typename GraphElement_>
    bool operator()(const GraphElement_& e) const {
      return get(m_pmap, e)<m_threshold;
    }
      
    property_map_t m_pmap;
    scalar_t m_threshold;
};

template <typename Graph_>
using graph_traits=boost::graph_traits<Graph_>;

template <typename VertexProp_, typename EdgeProp_>
using graph=boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProp_, EdgeProp_>;

template <typename VertexProp_, typename EdgeProp_>
using digraph=boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProp_, EdgeProp_>;

template <typename Graph_>
using vertex=typename boost::graph_traits<Graph_>::vertex_descriptor;

template <typename Graph_>
using edge=typename boost::graph_traits<Graph_>::edge_descriptor;

template <typename Graph_>
using vertex_iterator=typename boost::graph_traits<Graph_>::vertex_iterator;

template <typename Graph_>
using edge_iterator=typename boost::graph_traits<Graph_>::edge_iterator;

template <typename Graph_>
using out_edge_iterator=typename boost::graph_traits<Graph_>::out_edge_iterator;

template <typename Graph_, typename Scalar_, typename PropertyMap_>
using filtered_graph=boost::filtered_graph<Graph_, property_filter<Scalar_, PropertyMap_> >;

template <typename Graph_>
using edge_weight_map=typename boost::property_map<Graph_, boost::edge_weight_t>::type;
    
template <typename Graph_>
using weight_filtered_graph=boost::filtered_graph<Graph_, edge_weight_map<Graph_> >;
    
template <typename Graph_, typename Tag_>
using property_map=typename boost::property_map<Graph_, Tag_>::type;

template <typename GraphIn_, typename GraphOut_>
using vertex_prop_copier=boost::detail::vertex_copier<GraphIn_, GraphOut_>;

template <typename GraphIn_, typename GraphOut_>
using edge_prop_copier=boost::detail::edge_copier<GraphIn_, GraphOut_>;
    
} // namespace spurt

#endif // __GRAPH_DEFINITIONS_HPP__