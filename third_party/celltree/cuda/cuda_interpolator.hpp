#include <cassert>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

// #include <cuda/mem.hpp>
// #include <cuda/copy.hpp>
// #include <cuda/interval_timer.hpp>
#include "invert.hpp"

#include <mesh.hpp>
#include <variable.hpp>
#include <celltree.hpp>
#include <cutil_math.h>

// -------------------------------------------------------------------------

struct cuda_interpolator
{
    thrust::device_vector<unsigned int>   m_cells;
    thrust::device_vector<celltree::node> m_nodes;
    thrust::device_vector<float3>         m_pts;
    thrust::device_vector<float3>         m_var;

    cuda_interpolator( const mesh* m, 
                       const variable* v, 
                       const celltree* ct )
    {
        using namespace thrust;
        
        // reorder celltree nodes 
        std::vector<celltree::node> nodes( ct->nodes.size() );
        std::vector<celltree::node>::iterator ni, nn;

        nodes[0] = ct->nodes[0];

        for( ni=nodes.begin(), nn=nodes.begin()+1; ni!=nodes.end(); ++ni )
        {
            if( ni->is_leaf() )
                continue;
            
            assert( nn != nodes.end() );

            *(nn++) = ct->nodes[ni->left()];
            *(nn++) = ct->nodes[ni->right()];

            ni->set_children( nn-nodes.begin()-2 );
        }

        // inline cell types with indices
        std::vector<unsigned int> cells( m->nindices );
        std::vector<unsigned int>::iterator ci = cells.begin();

        for( ni=nodes.begin(); ni!=nodes.end(); ++ni )
        {
            if( ni->is_node() )
                continue;
        
            const unsigned int* begin = &(ct->leaves.front()) + ni->start();
            const unsigned int* end   = begin + ni->size();
        
            ni->st = ci - cells.begin();
                
            for( ; begin!=end; ++begin )
            {
                const mesh::cell& cell = m->cells[*begin];
                unsigned char k = cell.kind;
            
                const unsigned int* iind = m->indices + cell.start;
            
                for( unsigned int i=0; i<cell_size(cell.kind); ++i )
                {
                    *(ci++) = *(iind++) << 2 | (cell.kind & 3);
                    k >>= 2;
                }
            }
        }

        assert( ci == cells.end() );

        // output some stats
        fprintf( stdout, "celltree GPU mem: %lu bytes points\n"
                         "                  %lu bytes cells\n"
                         "                  %lu bytes nodes\n"
                         "                  %lu bytes variable\n"
                         "           total: %lu bytes\n\n",
                 (unsigned long)(m->npoints   * sizeof(float3)),
                 (unsigned long)(cells.size() * sizeof(unsigned int)),
                 (unsigned long)(nodes.size() * sizeof(celltree::node)),
                 (unsigned long)(m->npoints   * sizeof(float3)),
                 (unsigned long)(m->npoints   * sizeof(float3)) +
                 (unsigned long)(cells.size() * sizeof(unsigned int)) +
                 (unsigned long)(nodes.size() * sizeof(celltree::node)) +
                 (unsigned long)(m->npoints   * sizeof(float3)) );

        // copy to device (not optimal, should use pinned memory, but ok for now)
        const float3* pts = (const float3*)m->points;
        const float3* var = (const float3*)v->data;

        m_pts = device_vector<float3>( pts, pts + m->npoints );
        m_var = device_vector<float3>( var, var + m->npoints );

        m_nodes = nodes;
        m_cells = cells;
    }

    struct interpolator
    {
        float3*         m_pts;
        float3*         m_var;
        celltree::node* m_nodes;
        unsigned int*   m_cells;
        
        __device__ float3 operator()( const float3& pos ) const
        {
            unsigned int stack[32];
            
            unsigned int si = 0;
            unsigned int ni = 0;
            
            float3 c, pts[8];
            celltree::node n;
            unsigned int ind[8], kind;
            
            while( true )
            {
                n = m_nodes[ni];
    
                if( (n.index & 3) == 3 )
                {
                    for( unsigned int i=0, ii=n.st; i<n.sz; ++i )
                    {
                        ind[0] = m_cells[ii++];
                        kind   = ind[0] & 3;
                        ind[0] = ind[0] >> 2;
    
                        switch( kind )
                        {
                        case TETRAHEDRON:
    
                            for( unsigned int j=1; j<4; ++j )
                                ind[j] = m_cells[ii++] >> 2;
                            for( unsigned int j=0; j<4; ++j )
                                pts[j] = m_pts[ind[j]];
    
                            if( invert_tet( pts, c, pos ) )
                            {
                                for( unsigned int j=0; j<4; ++j )
                                    pts[j] = m_var[ind[j]];
                                
                                return interpolate_tet( pts, c );
                            }

                            break;
                        
                        case HEXAHEDRON:
                        
                            for( unsigned int j=1; j<8; ++j )
                                ind[j] = m_cells[ii++] >> 2;
                            for( unsigned int j=0; j<8; ++j )
                                pts[j] = m_pts[ind[j]];
                        
                            if( invert_hex( pts, c, pos ) )
                            {
                                for( unsigned int j=0; j<8; ++j )
                                    pts[j] = m_var[ind[j]];

                                return interpolate_hex( pts, c );
                            }

                            break;

                        case PRISM:

                            for( unsigned int j=1; j<6; ++j )
                                ind[j] = m_cells[ii++] >> 2;
                            for( unsigned int j=0; j<6; ++j )
                                pts[j] = m_pts[ind[j]];

                            if( invert_prism( pts, c, pos ) )
                            {
                                for( unsigned int j=0; j<6; ++j )
                                    pts[j] = m_var[ind[j]];

                                return interpolate_prism( pts, c );
                            }

                            break;

                        case PYRAMID:

                            for( unsigned int j=1; j<5; ++j )
                                ind[j] = m_cells[ii++] >> 2;
                            for( unsigned int j=0; j<5; ++j )
                                pts[j] = m_pts[ind[j]];

                            if( invert_pyramid( pts, c, pos ) )
                            {
                                for( unsigned int j=0; j<5; ++j )
                                    pts[j] = m_var[ind[j]];
                                    
                                return interpolate_pyramid( pts, c );
                            }

                            break;
                        }
                    }

                    if( si==0 )
                        break;

                    ni = stack[--si];
                }
                else
                {
                    unsigned int dim  = n.index & 3;
                    unsigned int left = n.index >> 2;

                    float p = dim == 0 ? pos.x : (dim == 1 ? pos.y : pos.z);

                    bool l = p <= n.lm;
                    bool r = p >  n.rm;
        
                    if( l && r )
                    {    
                        if( n.lm-p > p-n.rm )
                        {
                            ni = left;
                            stack[si++] = left + 1;
                        }
                        else
                        {
                            ni = left + 1;
                            stack[si++] = left;
                        }
                    }
                    else if( l )
                        ni = left;
                    else if( r )
                        ni = left + 1;
                    else if( si > 0 )
                        ni = stack[--si];
                    else
                        break;
                }
            }
            
            return make_float3( 0, 0, 0 );
        }
    };
    
    interpolator get_interpolator()
    {
        interpolator intp;
        
        intp.m_pts = thrust::raw_pointer_cast( &m_pts[0] );
        intp.m_var = thrust::raw_pointer_cast( &m_var[0] );
        intp.m_cells = thrust::raw_pointer_cast( &m_cells[0] );
        intp.m_nodes = thrust::raw_pointer_cast( &m_nodes[0] );

        return intp;
    }
};
