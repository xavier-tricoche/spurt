#ifndef __FDefaultInterpolator_hpp
#define __FDefaultInterpolator_hpp

#include <FTensorField.hh>

class FDefaultInterpolator
{
public:
    
    typedef double value_type;
    typedef double coord_type;

    FDefaultInterpolator( const shared_ptr<FTensorField> tf );
    
    bool operator()( const double& time, 
                     const double* pos, 
                     double* result ) const;

protected:    
    
    const shared_ptr<FTensorField> m_tf;
};


#endif // __FDefaultInterpolator_hpp