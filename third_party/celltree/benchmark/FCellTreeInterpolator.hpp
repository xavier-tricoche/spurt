#ifndef __FCellTreeInterpolator_hpp
#define __FCellTreeInterpolator_hpp

#include <FTensorField.hh>

class celltree;

class FCellTreeInterpolator
{
public:

    typedef double value_type;
    typedef double coord_type;

    FCellTreeInterpolator( const shared_ptr<FTensorField> tf );

    bool operator()( const double& time, 
                     const double* pos,
                     double* result ) const;
   
protected:    
    
    const shared_ptr<FTensorField> m_tf;
    shared_ptr<celltree>           m_ct;
};

#endif // __FCellTreeInterpolator_hpp
