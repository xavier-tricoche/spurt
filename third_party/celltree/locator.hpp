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

#ifndef __locator_hpp
#define __locator_hpp

#include <iostream>
#include <map>

#include "celltree.hpp"
#include "mesh_traits.hpp"
#include "mesh.hpp"

namespace CellTree {

class locator
{
    typedef mesh_traits<mesh> mtraits;

public:

    locator( const mesh* m, const celltree& ct ) :
        m_mesh(m), m_ct(ct), m_nodes(0), m_cells(0), m_evals(0), m_missed(0)
    {
    }

    bool find_cell( const double* pos,
                    mtraits::size_type& cell,
                    mtraits::coord_type* coord ) const
    {
        using namespace CellTree;
        ++m_evals;

        celltree::point_traversal pt( m_ct, pos, &m_nodes );

        unsigned int cnt = 0;

        while( const celltree::node* n = pt.next() )
        {
            const uint32_t* begin = &m_ct.leaves[n->start()];
            const uint32_t* end   = begin + n->size();

            for( ; begin!=end; ++begin )
            {
                ++m_cells;

                double min[3], max[3];
                mtraits::minmax( *m_mesh, *begin, min, max );

                if( pos[0] < min[0] || pos[0] > max[0] ||
                    pos[1] < min[1] || pos[1] > max[1] ||
                    pos[2] < min[2] || pos[2] > max[2] )
                    continue;

                ++cnt;

                if( mtraits::invert( *m_mesh, *begin, pos, coord ) )
                {
                    cell = *begin;

#ifdef HISTOGRAM
                    m_stats[cnt]++;
#endif

                    return true;
                }
            }
        }

        #ifdef HISTOGRAM
        m_stats[cnt]++;
        #endif

        ++m_missed;
        return false;
    }

    void get_stats( unsigned int& evals, unsigned int& nodes, unsigned int& cells, unsigned int& missed )
    {
        evals = m_evals;
        cells = m_cells;
        nodes = m_nodes;
        missed = m_missed;
    }

    void write_stats( std::ostream& out )
    {
#ifdef HISTOGRAM
        if( m_stats.size() )
        {
            unsigned int max = (--m_stats.end())->first;
            unsigned int sum = 0, csum = 0;

            std::vector<double> frac( max );

            for( std::map<int,int>::iterator mi=m_stats.begin(); mi!=m_stats.end(); ++mi )
            {
                sum += mi->second;
                csum += mi->first * mi->second;
            }

            out << '\n';

            for( unsigned int i=0; i<=max; ++i )
            {
                out << i << ", " << m_stats[i] << "\n";
                sum += m_stats[i];
            }

            out << '\n';

            printf( "sum = %u\n", sum );
        }
#endif
    }


    void reset_stats()
    {
        m_cells = 0;
        m_nodes = 0;
        m_missed = 0;
        m_evals = 0;
#ifdef HISTOGRAM
        m_stats.clear();
#endif
    }

protected:

    const mesh*             m_mesh;
    const celltree&         m_ct;

    mutable unsigned int      m_cells;
    mutable unsigned int      m_nodes;
    mutable unsigned int      m_evals;
    mutable unsigned int      m_missed;

#ifdef HISTOGRAM
    mutable std::map<int,int> m_stats;
#endif
};

} // namespace CellTree


#endif // __locator_hpp
