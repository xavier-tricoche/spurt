#ifndef __CONNECTIVITY_HPP__
#define __CONNECTIVITY_HPP__

#include <vector>
#include <graph/definitions.hpp>
#include <data/locator.hpp>
#include <misc/progress.hpp>

namespace {
template<typename Pos_>
struct __norm {};

template<typename T, size_t N>
struct __norm< nvis::fixed_vector<T, N> > {
    typedef nvis::fixed_vector<T, N> vector_t;
    T operator()(const vector_t& v) const {
        return nvis::norm(v);
    }
};

template<typename Pos_, typename Point_>
struct closer_than {
    closer_than() : from(0,0,0) {}
    closer_than(const Pos_& _from) : from(_from) {}
    closer_than(const Point_& _point) : from(_point.coordinate()) {}

    bool operator()(const Point_& a, const Point_& b) const {
        return length(a.coordinate()-from) < 
                          length(b.coordinate()-from);
    }

    Pos_ from;
    __norm<Pos_> length;
};
}

namespace spurt {
template <typename Graph_, typename Position_, size_t N>
void compute_connectivity(Graph_& graph, double radius, const std::vector<Position_>& vertices) {
    
    typedef Graph_ graph_t;
    typedef Position_ vertex_t;
    typedef typename Position_::value_type scalar_t;
    typedef spurt::point_locator<scalar_t, size_t, N> locator_t;
    typedef typename locator_t::point_type point_t;
    typedef typename std::list<point_t> neighborhood_t;
    typedef closer_than<vertex_t, point_t> less_dist;
    
    locator_t locator;
    for (size_t i=0; i<vertices.size(); ++i) {
        const vertex_t& v=vertices[i];
        locator.insert(point_t(v, i));
    }
    
    graph.clear();
    
    spurt::ProgressDisplay progress;
    progress.start(locator.size());
    size_t n=0;
    __norm<vertex_t> length;
    std::for_each(locator.begin(), locator.end(), 
                  [&] (const point_t& p) {
        neighborhood_t neighbors;
        size_t i=p.data();
        const vertex_t& x=p.coordinate();
        locator.find_within_range(neighbors, x, radius);
        neighbors.sort(less_dist(x));
        typename neighborhood_t::iterator it=neighbors.begin();
        for (++it; it!=neighbors.end(); ++it) {
            scalar_t dist=length(x-it->coordinate());
            if (dist>radius) break;
            boost::add_edge(i, it->data(), dist, graph);
        }
        progress.update(n++);
    });
    progress.end();
}

namespace {
class no_predicate {
public:
    template<typename T>
    bool operator()(const T& t) const { return true; }
};
}

template <typename GraphIn_, typename GraphOut_, 
          typename VertexPredicate_=no_predicate, 
          typename EdgePredicate_=no_predicate>
void extract_subgraph(GraphOut_& subg, GraphIn_& g, 
                      const VertexPredicate_& vp=no_predicate(), 
                      const EdgePredicate_& ep=no_predicate()) {
    typedef GraphIn_ graph_input_t;
    typedef GraphOut_ graph_output_t;
    typedef VertexPredicate_ vertex_pred_t;
    typedef EdgePredicate_ edge_pred_t;
    typedef typename boost::graph_traits<graph_input_t>::vertex_descriptor vertex_t;
    typedef typename boost::graph_traits<graph_input_t>::vertex_iterator vertex_iter_t;
    typedef typename boost::graph_traits<graph_input_t>::edge_iterator edge_iter_t;
    typedef typename boost::graph_traits<graph_output_t>::edge_descriptor edge_t;
    typedef spurt::vertex_prop_copier<graph_input_t, graph_output_t> vertex_prop_copy_t;
    typedef spurt::edge_prop_copier<graph_input_t, graph_output_t> edge_prop_copy_t;
        
    subg.clear();
    vertex_iter_t vit, vend;
    edge_iter_t eit, eend;
    vertex_prop_copy_t vertex_prop_copier(g, subg);
    edge_prop_copy_t edge_prop_copier(g, subg);
    
    // compute re-indexing
    std::map<int, int> old2new;
    auto index_map=boost::get(boost::vertex_index, subg);
    {
        int n=0;
        for (boost::tie(vit, vend)=boost::vertices(g); vit!=vend; ++vit) {
            if (!vp(*vit)) continue;
            old2new[index_map[*vit]]=n++;
        }
    }
    
    for (boost::tie(eit, eend)=boost::edges(g) ; eit!=eend; ++eit) {
        vertex_t a, b, a_copy, b_copy;
        a = boost::source(*eit, g);
        b = boost::target(*eit, g);
        auto aid=index_map[a];
        auto bid=index_map[b];
        if (!vp(a) || !vp(b) || !ep(*eit)) continue; // skip that edge
        auto edge_copy = boost::add_edge(old2new[a], old2new[b], subg);
        edge_prop_copier(*eit, edge_copy.first);
        a_copy = boost::source(edge_copy.first, subg);
        b_copy = boost::target(edge_copy.first, subg);
        vertex_prop_copier(a, a_copy);
        vertex_prop_copier(b, b_copy);
    }
}

template<typename Graph_, typename Int_>
size_t connected_components(std::vector<Int_>& cc_ids,Graph_& graph) {
    cc_ids.resize(boost::num_vertices(graph));
    return boost::connected_components(graph, &cc_ids[0]);
}

} // namespace spurt

#endif

