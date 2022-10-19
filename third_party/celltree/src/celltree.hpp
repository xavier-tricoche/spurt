#ifndef __celltree_hpp
#define __celltree_hpp

#include <stdint.h>
#include <float.h>
#include <vector>
#include <cmath>
#include <cstdio>

#ifndef _HD
#define _HD
#endif 

namespace celltree {

struct celltree
{
    struct node
    {   
        uint32_t index;
    
        union {
            struct {
                float lm;
                float rm;
            };
        
            struct {
                uint32_t sz;
                uint32_t st;
            };
        };

        void make_node( uint32_t left, unsigned int d, float b[2] )
        {
            index = (d & 3) | (left << 2);
            lm = b[0];
            rm = b[1];
        }

        void set_children( uint32_t left )
        {
            index = dim() | (left << 2);
        }

        
        bool is_node() const
        {
            return (index & 3) != 3;
        }

        
        unsigned int left() const
        {
            return (index >> 2);
        }

         
        unsigned int right() const
        {
            return (index >> 2) + 1;
        }

        
        unsigned int dim() const
        {
            return index & 3;
        }
    
        
        const float& lmax() const
        {
            return lm;
        }
        
        
        const float& rmin() const
        {
            return rm;
        }

        // ---
    
        void make_leaf( uint32_t start, uint32_t size )
        {
            index = 3;
            sz = size;
            st = start;
        }
    
        
        bool is_leaf() const
        {
            return index == 3;
        }
    
        
        unsigned int start() const
        {
            return st;
        }
    
        
        unsigned int size() const
        {
            return sz;
        }
    };

    std::vector<node>       nodes;
    std::vector<uint32_t>   leaves;
    
    size_t memsize() const
    {
        return nodes.size() * sizeof(node) + leaves.size() * sizeof(uint32_t);
    }
    
    template<typename V>
    bool traverse( const float* pos, const V& visitor ) const
    {
        unsigned int stack[64] = { 0 };
        unsigned int* sp = stack + 1;
        
        while( true )
        {
            const node& n = nodes[*(--sp)];
    
            if( n.is_leaf() )
            {
                const uint32_t* begin = &leaves[n.start()];
                const uint32_t* end   = begin + n.size();
                
                for( ; begin != end; ++begin )
                    if( visitor( *begin ) )
                        return true;
            }
            else
            {
                const float p = pos[n.dim()];
                const uint32_t left = n.left();

                bool l = p <= n.lmax();
                bool r = p > n.rmin();

                if( l && r )
                {    
                    if( n.lmax()-p < p-n.rmin() )
                    {
                        *(sp++) = left;
                        *(sp++) = left+1;
                    }
                    else
                    {
                        *(sp++) = left+1;
                        *(sp++) = left;
                    }
                }
                else if( l )
                    *(sp++) = left;
                else if( r )
                    *(sp++) = left+1;
            }
            
            if( sp == stack )
                return false;
        }
    }

    struct point_traversal
    {
        const celltree& m_ct;
        unsigned int    m_stack[32];
        unsigned int*   m_sp;
        const float*    m_pos;
    
        unsigned int*   m_nodecnt;

        point_traversal( const celltree& ct, const float* pos, unsigned int* nodecnt ) : 
            m_ct(ct), m_pos(pos), m_nodecnt(nodecnt)
        {
            m_stack[0] = 0;
            m_sp = m_stack + 1;
        }
    
        const celltree::node* next()
        {
            while( true )
            {
                if( m_sp == m_stack )
                    return 0;

                const celltree::node* n = &m_ct.nodes.front() + *(--m_sp);
    
                ++(*m_nodecnt);

                if( n->is_leaf() )
                    return n;
            
                const float p = m_pos[n->dim()];
                const uint32_t left = n->left();

                bool l = p <= n->lmax();
                bool r = p > n->rmin();

                if( l && r )
                {    
                    if( n->lmax()-p < p-n->rmin() )
                    {
                        *(m_sp++) = left;
                        *(m_sp++) = left+1;
                    }
                    else
                    {
                        *(m_sp++) = left+1;
                        *(m_sp++) = left;
                    }
                }
                else if( l )
                    *(m_sp++) = left;
                else if( r )
                    *(m_sp++) = left+1;
            }
        }
    };

#if 0
    struct ray_traversal
    {
        const celltree& m_ct;
        unsigned int    m_stack[64];
        unsigned int*   m_sp;
        const float*    m_roff;
        const float*    m_rdir;

        ray_traversal( const celltree& ct, const float* roff, const float* rdir ) :
            m_ct(ct), m_roff(roff), m_rdir(rdir)
        {
            m_stack[0] = 0;
            m_sp = m_stack + 1;
        }

        template<typename F>
        void traverse( F& f, float t0, float t1, unsigned int ni=0 ) 
        {
            const celltree::node n = m_ct.nodes[ni];
    
	    // printf( "traverse %u, %f %f: ", ni, t0, t1 );

            if( n.is_leaf() )
            {
		// printf( "leaf\n" );

                f( m_roff, t0, t1, m_rdir, &n );
                return;
            }

            const float d = m_rdir[n.dim()];
            const float o = m_roff[n.dim()];

            const unsigned int left = n.left();

	    // printf( "direction %u, ", n.dim() );

            if( d == 0 )
            {
		// printf( "straight, o = %f, lmax = %f, rmin = %f\n", o, n.lmax(), n.rmin() );

                if( o <= n.lmax() )
                    traverse( f, t0, t1, left );
                if( o >= n.rmin() )
                    traverse( f, t0, t1, left+1 );
            }
            else
            {
                float lt = (n.lmax() - o)/d;
                float rt = (n.rmin() - o)/d;

		// printf( "across, lt = %f, rt = %f\n", lt, rt );

                if( d>0 )
                {
                    if( lt > t0 && lt < t1 )
                        traverse( f, t0, lt, left );

                    if( rt > t0 && rt < t1 )
                        traverse( f, rt, t1, left+1 );
                }
                else
                {
                    if( rt > t0 && rt < t1 )
                        traverse( f, t0, rt, left+1 );

                    if( lt > t0 && lt < t1 )
                        traverse( f, lt, t1, left );
                }
            }
        }
    };
#endif
};

// -------------------------------------------------------------------------

#if 0
inline static void nccheck( int result )
{
    if( result != NC_NOERR )
        throw std::runtime_error( nc_strerror(result) );
}

void save_tree( const char* filename, const celltree& tree )
{
    int ncid;
    nccheck( nc_create( filename, NC_NETCDF4, &ncid ) );

    int ndim[2], nvar;
    nccheck( nc_def_dim( ncid, "nnodes", tree.nodes.size(), &ndim[0] ) );
    nccheck( nc_def_dim( ncid, "snodes", 3,                 &ndim[1] ) );
    nccheck( nc_def_var( ncid, "nodes", NC_UINT, 2, ndim, &nvar ) );

    int ldim, lvar;
    nccheck( nc_def_dim( ncid, "nleaves", tree.leaves.size(), &ldim ) );
    nccheck( nc_def_var( ncid, "leaves", NC_UINT, 1, &ldim, &lvar ) );

    nccheck( nc_put_var_uint( ncid, nvar, (unsigned int*)&tree.nodes.front() ) );
    nccheck( nc_put_var_uint( ncid, lvar, (unsigned int*)&tree.leaves.front() ) );

    nc_close( ncid );
}
#endif 

} // namespace celltree

#endif // __celltree_hpp
