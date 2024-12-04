#include <list>
#include <utility>
#include <vector>

namespace spurt {

template<typename VAttribute_, typename EAttribute_>
class Graph : public std::vector<std::pair<VAttribute_,
                                 std::list<std::pair<EAttribute_,int> > > > {
public:
    typedef VAttribute_                              v_attribute_t;
    typedef EAttribute_                              e_attribute_t;    
    typedef std::pair<e_attribute_t, int>            edge_t;
    typedef std::list<edge_t>                        neighborhood_t;
    typedef std::pair<v_attribute_t, neighborhood_t> vertex_t;
    
    Graph() : adjacency(), directed(false) {}
    Graph(size_t nvertices) : adjacency(nvertices), directed(false) {}
    
    void add_vertex(v_attribute_t att) {
        adjacency.push_back(vertex_t(att, neighboorhood_t()));
    }
    
    void add_edge(int from, int to, e_attribute_t att) {
        assert(from < num_verts() && to < num_verts());
        adjacency[from].push_back(edge_t(att, to));
        if (!directed) {
            adjacency[to].push_back(edge_t(att, from));
        }
    }
    
    vertex_t& vertex(int v) { return adjacency[v]; }
    const vertex_t& vertex(int v) const { return adjacency[v]; }
    
    neighborhood_t& neighborhood(int v) { return adjacency[v].second; } 
    const neighborhood_t& neighborhood(int v) const { return adjacency[v].second; }

    size_t num_verts() const { return adjacency.size(); }
    size_t num_edges() const { 
        size_t sz=0;
        std::for_each(adjacency.begin(), adjacency.end(), [&](const vertex_t& v)
        {
            sz += v.second.size();
        });
        return sz / (directed ? 1 : 2); 
    }
    
    // iterators
private:
    std::vector<vertex_t> adjacency;
    bool directed;
};

class DFS {
    std::vector<bool> _marked; 
    std::list<int> _visited;
    
    void dfs(int s) {
        _marked[s] = true;
        _visited.push_back(s);
        for (auto it=_graph.neighbors[s].begin() ; 
            it!=_graph.neighbors[s].end() ; ++it) {
            if (!_marked[it->first]) {
                dfs(it->first);
            }
        } 
    }
public:
    DFS(const Graph& graph) : _graph(graph) {}
    
    void search(int s) {
        _marked.resize(graph.size());
        std::fill(_marked.begin(), _marked.end(), false);
        _visited.clear();
        dfs(_graph, s);
    }
    
    const std::list<int>& visited() const { return _visited; }
    bool visited(int v) const { return _marked[v]; }
}; 

} // namespace spurt