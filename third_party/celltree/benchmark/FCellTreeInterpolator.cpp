#include "FCellTreeInterpolator.hpp"
#include "fantom_mesh_traits.hpp"

#include <celltree.hpp>
#include <celltree_builder.hpp>

// -------------------------------------------------------------------------

FCellTreeInterpolator::FCellTreeInterpolator( const shared_ptr<FTensorField> tf ) 
    : m_tf(tf)
{
    m_ct = shared_ptr<celltree>( new celltree );
    
    celltree_builder builder;
    builder.build( *m_ct, tf );
}

// -------------------------------------------------------------------------

bool 
FCellTreeInterpolator::operator()( const double& time, 
                                   const double* pos, 
                                   double* result ) const
{
    FTensor r;
    FPosition p( pos[0], pos[1], pos[2] );

    unsigned int tmp;

    const float _x[3] = { pos[0], pos[1], pos[2] };
    celltree::point_traversal pt( *m_ct, _x, &tmp );

    while( const celltree::node* n = pt.next() )
    {
        const unsigned int* begin = &(m_ct->leaves[n->start()]);
        const unsigned int* end   = begin + n->size();

        for( ; begin!=end; ++begin )
        {
            shared_ptr<FCell> cell = m_tf->getCell( *begin );

            FBoundingBox box = cell->getBoundingBox();

            if( /*box.isInside( p ) &&*/ cell->isInside( p ) )
            {
                // cell->setTensors( m_tf->getTensorSet().get() );
                cell->interpolate( r, p );

                result[0] = r(0);
                result[1] = r(1);
                result[2] = r(2);
                return true;
            }
        }
    }    
    
    return false;
}
