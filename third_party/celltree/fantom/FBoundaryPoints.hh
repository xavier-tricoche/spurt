//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FBoundaryPoints.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/11 17:41:49 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FBoundaryPoints_hh
#define __FBoundaryPoints_hh

#include "FFeature.hh"
#include <vector>

struct FBoundaryPoints : public FFeature
{  
    FBoundaryPoints( vector< bool >& p_input, vector< bool >& c_input ) 
	: FFeature() 
	{ 
	    p_boundary.swap( p_input );
	    c_boundary.swap( c_input );
	};

    virtual std::string getType() const { return "boundarypoints"; };

    vector< bool > p_boundary;
    vector< bool > c_boundary;

    virtual ~FBoundaryPoints() {};
};

#endif // __FBoundaryPoints_hh
