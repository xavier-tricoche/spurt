#include "FDefaultInterpolator.hpp"

#include <FGrid3DArbitrary.hh>

// -------------------------------------------------------------------------

FDefaultInterpolator::FDefaultInterpolator( const shared_ptr<FTensorField> tf )
    : m_tf(tf)
{
    shared_ptr<FGrid3DArbitrary> g3da = 
        shared_dynamic_cast<FGrid3DArbitrary>( m_tf->getGrid() );

    if( g3da )
        g3da->buildKdTree();
}

// -------------------------------------------------------------------------
    
bool
FDefaultInterpolator::operator()( const double& time, 
                                  const double* pos, 
                                  double* result ) const
{
    FTensor r;
    FPosition p( pos[0], pos[1], pos[2] );

    if( m_tf->interpolate( r, p ) )
    {
        result[0] = r(0);
        result[1] = r(1);
        result[2] = r(2);
        return true;
    }
    
    return false;
}

