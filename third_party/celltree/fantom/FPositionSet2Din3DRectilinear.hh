//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile$
// Language:  C++
// Date:      $Date$
// Author:    $Author$
// Version:   $Revision$
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet2Din3DRectilinear_hh
#define __FPositionSet2Din3DRectilinear_hh

#include "FPositionSet2DRectilinear.hh"
#include <vector>
#include "FMatrix.hh"

/** 
    the class has the functionality its name says
*/
class FPositionSet2Din3DRectilinear : public FPositionSet2DRectilinear
{
    //  MAKE_SERIALIZABLE( FPositionSet2Din3DRectilinear );

    //==== Constructors ========================================================
public:

    /** 
     *\par Description:
     *Constructor
     *\pre
     *Xcoords and YCoords save the x/y -coordinate 
     *values in increasing order
     *\post 
     *Xcoords and YCoords are empty vectors 
     *(because of swap operation)
     *\exception
     *none
     *\param
     *none
     *\param XCoords: vector with x coordinate values
     *\param YCoords: vector with y coordinate values
     */
    FPositionSet2Din3DRectilinear( vector<double>& XCoords, vector<double>& YCoords,
				   const FMatrix& _trans );

    //=== Member Functions ====================================================

    //=========================================================================
    //=== QUERY FUNCTIONS  ====================================================
    //=========================================================================

    /** 
     *\par Description:
     *Gets the class name.
     *\pre
     *none
     *\post
     *none
     *\exception
     *none
     */
    virtual const FString& getClassName() const;

    void get3DPosition( FArray& pos3D, const FIndex& id ) const;
    void get2DPosition( FArray& pos2D, const FArray& pos3D ) const;
    
protected:

    FMatrix trans;
};

#endif // __FPositionSet_hh
