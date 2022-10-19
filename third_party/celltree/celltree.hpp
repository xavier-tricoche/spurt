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

namespace CellTree {

struct celltree
{
    struct node
    {
        uint32_t index;

        union {
            struct {
                double lm;
                double rm;
            };

            struct {
                uint32_t sz;
                uint32_t st;
            };
        };

        void make_node( uint32_t left, unsigned int d, double b[2] )
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


        const double& lmax() const
        {
            return lm;
        }


        const double& rmin() const
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

    bool write(const std::string& filename) const {
        std::fstream output(filename, std::ios::out | std::ios::binary);
        if (!output) {
            std::cerr << "unable to open " << filename << " to export celltree\n";
            return false;
        }
        size_t nnodes=nodes.size();
        size_t nleaves=leaves.size();
        output.write(reinterpret_cast<const char*>(&nnodes), sizeof(size_t));
        output.write(reinterpret_cast<const char*>(&nleaves), sizeof(size_t));
        output.write(reinterpret_cast<const char*>(&nodes[0]), nnodes*sizeof(node));
        output.write(reinterpret_cast<const char*>(&leaves[0]), nleaves*sizeof(uint32_t));
        output.close();
        return true;
    }

    bool read(const std::string& filename) {
        std::fstream input(filename, std::ios::in | std::ios::binary);
        if (!input) {
            std::cerr << "unable to open " << filename << " to import celltree\n";
            return false;
        }
        size_t n;
        input.read(reinterpret_cast<char*>(&n), sizeof(size_t));
        nodes.resize(n);
        input.read(reinterpret_cast<char*>(&n), sizeof(size_t));
        leaves.resize(n);
        input.read(reinterpret_cast<char*>(&nodes[0]), nodes.size()*sizeof(node));
        input.read(reinterpret_cast<char*>(&leaves[0]), leaves.size()*sizeof(uint32_t));
        return true;
    }

    template<typename V>
    bool traverse( const double* pos, const V& visitor ) const
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
                const double p = pos[n.dim()];
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
        const double*   m_pos;

        unsigned int*   m_nodecnt;

        point_traversal( const celltree& ct, const double* pos, unsigned int* nodecnt ) :
            m_ct(ct), m_pos(pos), m_nodecnt(nodecnt)
        {
            m_stack[0] = 0;
            m_sp = m_stack + 1;
        }

        const node* next()
        {
            while( true )
            {
                if( m_sp == m_stack )
                    return 0;

                const node* n = &m_ct.nodes.front() + *(--m_sp);

                ++(*m_nodecnt);

                if( n->is_leaf() )
                    return n;

                const double p = m_pos[n->dim()];
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
        const double*    m_roff;
        const double*    m_rdir;

        ray_traversal( const celltree& ct, const double* roff, const double* rdir ) :
            m_ct(ct), m_roff(roff), m_rdir(rdir)
        {
            m_stack[0] = 0;
            m_sp = m_stack + 1;
        }

        template<typename F>
        void traverse( F& f, double t0, double t1, unsigned int ni=0 )
        {
            const CellTree::node n = m_ct.nodes[ni];

        // printf( "traverse %u, %f %f: ", ni, t0, t1 );

            if( n.is_leaf() )
            {
        // printf( "leaf\n" );

                f( m_roff, t0, t1, m_rdir, &n );
                return;
            }

            const double d = m_rdir[n.dim()];
            const double o = m_roff[n.dim()];

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
                double lt = (n.lmax() - o)/d;
                double rt = (n.rmin() - o)/d;

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

} // namespace CellTree

#endif // __celltree_hpp
