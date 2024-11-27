#ifndef __CONNECTIVITY_HPP__
#define __CONNECTIVITY_HPP__

#include <vector>
#include <graph/definitions.hpp>
#include <data/locator.hpp>
#include <misc/progress.hpp>
#include <misc/meta_utils.hpp>

template<typename Pos_, typename Point_>
struct closer_than {
    typedef spurt::data_traits<Pos_> traits_type;

    closer_than() : from(0,0,0) {}
    closer_than(const Pos_& _from) : from(_from) {}
    closer_than(const Point_& _point) : from(_point.position()) {}

    bool operator()(const Point_& a, const Point_& b) const {
        return traits_type::norm(a.position()-from) < 
                          traits_type::norm(b.position()-from);
    }

    Pos_ from;
};

namespace spurt {
    
template <typename Graph_, typename Position_, size_t N>
void compute_knn_connectivity(Graph_& graph, size_t knn, const std::vector<Position_>& vertices) {
    typedef Graph_ graph_t;
    typedef Position_ vertex_t;
    typedef data_traits<Position_> traits_t;
    typedef typename traits_t::value_type scalar_t;
    typedef spurt::point_locator<vertex_t, size_t> locator_t;
    typedef typename locator_t::point_type point_t;
    typedef typename std::list<point_t> neighborhood_t;
    typedef closer_than<vertex_t, point_t> less_dist;
    
    size_t npts = vertices.size();
    
    std::vector<point_t> data_pts(npts);
    for (size_t i=0; i<npts; ++i) {
        data_pts[i] = point_t(vertices[i], i);
    }
    locator_t locator(data_pts.begin(), data_pts.end());

    graph.clear();
    
    spurt::ProgressDisplay progress;
    progress.start(npts);
    size_t n=0;
    typedef std::vector<point_t> neighbors_t;
    std::for_each(data_pts.begin(), data_pts.end(), 
                  [&] (const point_t& p) {
        neighbors_t nns;
        size_t i = p.data();
        locator.find_n_nearest_points(nns, p.position(), knn);
        std::sort(nns.begin(), nns.end(), less_dist(p));
        typename neighbors_t::iterator it=nns.begin();
        for (++it; it!=nns.end(); ++it) {
            scalar_t dist=traits_t::norm(p.position()-it->position());
            boost::add_edge(i, it->data(), dist, graph);
        }
        progress.update(n++);
    });
    progress.end();
}


template <typename Graph_, typename Position_, size_t N>
void compute_connectivity(Graph_& graph, double radius, const std::vector<Position_>& vertices) {
    
    typedef Graph_ graph_t;
    typedef Position_ vertex_t;
    typedef data_traits<Position_> traits_t;
    typedef typename traits_t::value_type scalar_t;
    typedef spurt::point_locator<vertex_t, size_t> locator_t;
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
    std::for_each(locator.begin(), locator.end(), 
                  [&] (const point_t& p) {
        neighborhood_t neighbors;
        size_t i=p.data();
        const vertex_t& x=p.position();
        locator.find_within_range(neighbors, x, radius);
        neighbors.sort(less_dist(x));
        typename neighborhood_t::iterator it=neighbors.begin();
        for (++it; it!=neighbors.end(); ++it) {
            scalar_t dist=traits_t::norm(x-it->position());
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

