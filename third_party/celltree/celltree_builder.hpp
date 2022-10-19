/*************************************************************************
* celltree: boundary interval hierarchy tree-based solution for point 
*           location and interpolation in unstructured meshes. This
*           technique was first proposed in:
*           "Fast, Memory-Efficient Cell location in Unstructured Grids
*            for Visualization"
*                by Christop Garth and Kenneth I Joy, VisWeek 2011.
* Author: Christoph Garth
* Copyright (C) 2011 Christoph Garth
* All rights reserved
**************************************************************************/

#ifndef __celltree_builder_hpp
#define __celltree_builder_hpp

#include <limits>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "celltree.hpp"
#include "mesh_traits.hpp"
#include "progress.hpp"

namespace CellTree {

class celltree_builder
{
    typedef CellTree::celltree::node node;
private:
    
    struct bucket
    {
        double        min;
        double        max;
        unsigned int cnt;
    
        bucket()
        {
            cnt = 0;
            min =  std::numeric_limits<double>::max();
            max = -std::numeric_limits<double>::max();
        }
    
        void add( const double _min, const double _max )
        {
            ++cnt;
        
            if( _min < min )
                min = _min;
            
            if( _max > max )
                max = _max;
        }
    };

    struct box
    {
        double min[3];
        double max[3];
        
        box()
        {
            for( unsigned int d=0; d<3; ++d )
            {
                min[d] = std::numeric_limits<double>::max();
                max[d] = std::numeric_limits<double>::max();
            }
        }
        
        void add( const box& other )
        {
            for( unsigned int d=0; d<3; ++d )
            {
                if( other.min[d] < min[d] )  min[d] = other.min[d];
                if( other.max[d] < max[d] )  max[d] = other.max[d];
            }
        }
    };

    struct per_cell 
    {
        uint32_t  ind;
        double     min[3];
        double     max[3];
    };

    struct center_order
    {
        unsigned int d;
    
        center_order( unsigned int _d ) : 
            d(_d)
        {
        }

        bool operator()( const per_cell& pc0, const per_cell& pc1 )
        {
            return pc0.min[d] + pc0.max[d] < pc1.min[d] + pc1.max[d];
        }
    };

    struct left_predicate
    {
        unsigned int       d;
        double              p;
    
        left_predicate( unsigned int _d, double _p ) : 
            d(_d), p(_p)
        {
        }
   
        bool operator()( const per_cell& pc )
        {
            return (pc.min[d] + pc.max[d])/2.0 < p;
        }
    };


    // -------------------------------------------------------------------------

    void find_min_max( const per_cell* begin, const per_cell* end,  
        double* min, double* max )
    {
        if( begin == end )
            return;
            
        for( unsigned int d=0; d<3; ++d )
        {
            min[d] = begin->min[d];
            max[d] = begin->max[d];
        }
        
        while( ++begin != end )
        {
            for( unsigned int d=0; d<3; ++d )
            {
                if( begin->min[d] < min[d] )    min[d] = begin->min[d];
                if( begin->max[d] > max[d] )    max[d] = begin->max[d];
            }
        }
    }

    // -------------------------------------------------------------------------
    
    void find_min_d( const per_cell* begin, const per_cell* end,  
                     unsigned int d, double& min )
    {
        min = begin->min[d];
        
        while( ++begin != end )
            if( begin->min[d] < min )    
                min = begin->min[d];
    }

    void find_max_d( const per_cell* begin, const per_cell* end,  
                     unsigned int d, double& max )
    {
        max = begin->max[d];
        
        while( ++begin != end )
            if( begin->max[d] > max )    
                max = begin->max[d];
    }

    // -------------------------------------------------------------------------

    void recursive_split( unsigned int index, per_cell* begin, per_cell* end )
    {
        size_t size = end - begin;
    
        if( size < m_leafsize )
        {
            m_prog += size;
            return;
        }
        
        double cost = std::numeric_limits<double>::max();
        double plane;
        int   dim;
        double clip[2];

        double min[3], max[3], ext[3];
        find_min_max( begin, end, min, max );

        for( unsigned int d=0; d<3; ++d )
            ext[d] = max[d]-min[d];



        per_cell* mid = begin;

        const unsigned int nbuckets = m_buckets;

        bucket* b[3];
        b[0] = new bucket[nbuckets];
        b[1] = new bucket[nbuckets];
        b[2] = new bucket[nbuckets];
        // bucket b[3][nbuckets];

        for( const per_cell* pc=begin; pc!=end; ++pc )
        {
            for( unsigned int d=0; d<3; ++d )
            {
                if( ext[d] == 0 )
                    continue;

                double cen = 0.5*(pc->min[d] + pc->max[d]);

                int ind = (int)ceil( nbuckets * (cen-min[d])/ext[d] ) - 1;
                b[d][ind].add( pc->min[d], pc->max[d] );
            }
        }
        
        for( unsigned int d=0; d<3; ++d )
        {
            unsigned int sum = 0;

            for( unsigned int n=0; n<nbuckets-1; ++n )
            {
                // printf( "dim %u, bucket plane %f %f\n", d, min[d] + (n+1)*ext[d]/nbuckets, (min[d]+max[d])/2.0 );

                double lmax = -std::numeric_limits<double>::max();
                double rmin =  std::numeric_limits<double>::max();

                for( unsigned int m=0; m<=n; ++m )
                    if( b[d][m].max > lmax )
                        lmax = b[d][m].max;

                for( unsigned int m=n+1; m<nbuckets; ++m )
                    if( b[d][m].min < rmin )
                        rmin = b[d][m].min;

                sum += b[d][n].cnt;

                double lvol = (lmax-min[d])/ext[d];
                double rvol = (max[d]-rmin)/ext[d];

                double c = lvol*sum + rvol*(size-sum);

                if( sum > 0 && sum < size && c < cost )
                {
                    cost    = c;
                    dim     = d;
                    plane   = min[d] + (n+1)*ext[d]/nbuckets;
                    clip[0] = lmax;
                    clip[1] = rmin;
                }
            }
        }

        // fallback
        if( cost != std::numeric_limits<double>::max() )
            mid = std::partition( begin, end, left_predicate( dim, plane ) );

        if( mid == begin || mid == end )
        {
            dim = std::min_element( ext, ext+3 ) - ext;
            
            mid = begin + (end-begin)/2;
            std::nth_element( begin, mid, end, center_order( dim ) );

        }

        find_max_d( begin, mid, dim, clip[0] );
        find_min_d( mid,   end, dim, clip[1] );

        node child[2];
        child[0].make_leaf( begin - m_pc, mid-begin );
        child[1].make_leaf( mid   - m_pc, end-mid );
        
        //#pragma omp critical
        {
            m_nodes[index].make_node( m_nodes.size(), dim, clip );
            m_nodes.insert( m_nodes.end(), child, child+2 );
        }
        
        //#pragma omp task
        recursive_split( m_nodes[index].left(), begin, mid );
    
        //#pragma omp task
        recursive_split( m_nodes[index].right(), mid, end );
    }
     
public:

    celltree_builder()
    {
        m_buckets =  5;
        m_leafsize = 4;
    }
    
    template<typename M>
    void build( celltree& ct, const M& mesh )
    {   
        typedef mesh_traits<M> tr;
        
        const int size = tr::cells_size( mesh );
        
        m_pc = new per_cell[size];
        
        #pragma omp parallel for
        for( int i=0; i<size; ++i )
        {
            m_pc[i].ind = i;
            tr::minmax( mesh, i, m_pc[i].min, m_pc[i].max );
        }
                
        node root;
        root.make_leaf( 0, size );
        m_nodes.push_back( root );

        m_prog.restart( "building cell tree", size );

        #pragma omp parallel
        {
            #pragma omp single nowait
            recursive_split( 0, m_pc, m_pc+size );
        }
        
        m_prog.finished();
        
        ct.nodes.resize( m_nodes.size() );
        ct.nodes[0] = m_nodes[0];

        std::vector<node>::iterator ni = ct.nodes.begin();
        std::vector<node>::iterator nn = ct.nodes.begin()+1;

        for( ; ni!=ct.nodes.end(); ++ni )
        {
            if( ni->is_leaf() )
                continue;
            
            *(nn++) = m_nodes[ni->left()];
            *(nn++) = m_nodes[ni->right()];

            ni->set_children( nn-ct.nodes.begin()-2 );
        }
        
        ct.leaves.resize( size );

        for( int i=0; i<size; ++i )
            ct.leaves[i] = m_pc[i].ind;
        
        delete[] m_pc;
    }

public:

    unsigned int        m_buckets;
    unsigned int        m_leafsize;
    per_cell*           m_pc;
    progress            m_prog;
    std::vector<node>   m_nodes;
};

} // namespace CellTree


#endif // __celltree_builder_hpp
