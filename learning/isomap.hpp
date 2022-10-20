#ifndef __ISOMAP_HPP__
#define __ISOMAP_HPP__

#include <stdexcept>

// nvis
#include <util/timer.hpp>

// Boost graph library
#include <graph/definitions.hpp>
#ifndef _DONT_USE_SPARSE_SOLUTION_TO_APSP_PROBLEM_
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#endif

// Eigen
#include <Eigen/Dense>

// own
#include <learning/multidimensional_scaling.hpp>
#include <misc/progress.hpp>

namespace {
// Wrapper for Eigen matrices (required by BGL's implementation of
// Johnson's APSP algorithm)

template<typename Row_>
struct row_wrapper {
    typedef Row_ row_t;
    typedef typename Row_::Base::Scalar value_t;
    
    row_wrapper(row_t _r) : row(_r) {}
    value_t& operator[](int c) {
        return row(0, c);
    }
    const value_t& operator[](int c) const {
        return row(0, c);
    }
    
    row_t row;
};

template<typename Matrix_>
struct matrix_wrapper {
    typedef Matrix_ matrix_t;
    typedef typename matrix_t::RowXpr matrix_row_t;
    typedef row_wrapper<matrix_row_t> row_t;
    
    matrix_wrapper(Matrix_& _m) : m(_m) {}
    row_t operator[](int r) {
        return row_t(m.row(r));
    }
    const row_t operator[](int r) const {
        return row_t(m.row(r));
    }
        
    matrix_t& m;
};
}


namespace spurt {

template<typename Scalar_, typename Graph_>
class Isomap {
public:
    typedef Scalar_ scalar_t; 
    typedef Graph_ graph_t;
    typedef MultidimensionalScaling<Scalar_> mds_t;
    typedef typename mds_t::vector_t vector_t;
    typedef typename mds_t::row_vector_t row_vector_t;
    typedef typename mds_t::matrix_t matrix_t;

private:
    typedef boost::graph_traits<graph_t> traits_t;
    typedef boost::property_map<graph_t, boost::edge_weight_t> weight_map;
    typedef typename traits_t::vertex_descriptor vertex_t;
    typedef typename traits_t::edge_descriptor edge_t;
    typedef typename traits_t::out_edge_iterator edge_iterator_t;
    typedef typename weight_map::type weight_map_t;
    typedef typename weight_map::value_type weight_t;
    
    void all_pairs_distances(matrix_t& dist) {
        typedef size_t index_t;
        
        const size_t N=boost::num_vertices(M_graph);
        
        // initialize distance matrix
        dist.resize(N, N);
        
#ifndef _DONT_USE_SPARSE_SOLUTION_TO_APSP_PROBLEM_
        matrix_wrapper<matrix_t> D(dist);
        boost::johnson_all_pairs_shortest_paths(M_graph, D);
#else
        const scalar_t INFTY=std::numeric_limits<scalar_t>::max()/2-1;
        weight_map_t weight_map=boost::get(boost::edge_weight, M_graph);
        dist.setConstant(INFTY);
        edge_iterator_t eit, eend;    
        for (size_t i=0 ; i<N; ++i) {
            dist(i,i)=0;
            for (boost::tie(eit, eend)=boost::out_edges(i, M_graph); 
                 eit!=eend; ++eit) {            
                index_t j=boost::target(*eit, M_graph);
                dist(i,j)=weight_map[*eit];
            }
        }
        // std::cout << "After initialization, dist matrix =\n" << dist << '\n';
        
        // Floyd-Warshall method to compute shortest path between all pairs
        ProgressDisplay progress;
        std::cout << "Floyd-Warshall computation\n";
        progress.start(N*N*N);
        for (size_t k=0; k<N; ++k) {
            for (size_t i=0; i<N; ++i) {
                progress.update(N*(i+N*k));
                for (size_t j=0; j<N; ++j) {
                    dist(i,j)=std::min(dist(i,j), dist(i,k)+dist(k,j));
                }
            }
        }
        progress.end();
        // std::cout << "After Floyd-Warshall, dist matrix=\n" << dist << '\n';
#endif
    }

public:
    Isomap(graph_t& graph) : M_graph(graph) {
        std::vector<int> cc(boost::num_vertices(M_graph));
        int num=boost::connected_components(M_graph, &cc[0]);
        if (num>1) {
            std::cout << "ERROR: There are more than one (" << num 
                      << ") connected component!!\n\n";
            throw std::invalid_argument("Isomap constructor: Provided graph is not connected!");
        }
        else {
            std::cout << "Isomap: " 
                << boost::num_vertices(M_graph) << " vertices and "
                << boost::num_edges(M_graph) << " edges in input\n";
        }
    }
    
    template<typename I_=int>
    void embed() {
        nvis::timer t;
        all_pairs_distances(distances);
        std::cout << "all pair shortest path computation took " 
                  << t.elapsed() << " seconds.\n";
        t.restart();
        mds.embed(distances);
        std::cout << "MDS embedding took " << t.elapsed() << " seconds.\n";
    }
    
    vector_t eigenvalues(int ndim) const {
        return mds.eigenvalues().tail(ndim);
    } 
    
    matrix_t coordinates(int ndim) const {
        return mds.coordinates(ndim);
    }
    
    scalar_t l2_error(int ndim) const {
        return mds.l2_error(ndim);
    }
    
private:
    graph_t M_graph;
    mds_t mds;
    matrix_t distances;
};

} // namespace spurt

#endif