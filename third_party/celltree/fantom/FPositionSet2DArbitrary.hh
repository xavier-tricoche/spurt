//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet2DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:01 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet2DArbitrary_hh
#define __FPositionSet2DArbitrary_hh

#include <boost/shared_ptr.hpp>
#include "FPositionSet.hh"

using namespace boost;
// forward-declaration
class FGrid2DArbitrary;

/** 
 *   the class has the functionality its name says
 */
class FPositionSet2DArbitrary: public FPositionSet
{
    //==== Constructors ========================================================
public:

    /** 
     *\par Description:
     *Constructor (float-based)
     *\pre
     *coords contains the x and y component for each position.
     *( pos0[x], pos0[y], pos1[x], pos1[y] ... )
     *\post 
     *coords has become an empty vector (because of swap operation)
     *\exception
     *none
     *\param
     *none
     *\param 
     */
    FPositionSet2DArbitrary( vector<double>& coords );
    FPositionSet2DArbitrary( shared_ptr<FanyArray<double> > coords );


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
    const FString& getClassName() const;

    /** 
     *\par Description:
     *Gets a position.
     *\pre
     *finalize() has been invoked
     *0 <= posId < nbPositions
     *\post
     *none
     *\exception
     *FIndexOutOfBoundsException
     *\param
     *resultPos FPosition to store the searched position (or vector to receive
     *the double coordinates).
     *\param
     *pIndex: index of Position to get
     */           
    void getPosition(FPosition& resultPos, const FIndex& pIndex) const;


    /** 
     *\par Description:
     *Gets the components of a position
     *\pre
     *0 <= posId < nbPositions; resultPos must be at least size 2
     *\post
     *none
     *\exception
     *FIndexOutOfBoundsException
     *\param
     *resultPos vector of doubles
     *the double coordinates).
     *\param
     *pIndex: index of Position to get
     */           
    void getPosition(vector<double>& resultPos, const FIndex& pIndex) const;

    positive getNbPositions() const;

    positive memSize() const;

    friend class FGrid2DArbitrary;

    virtual FanyArray< double >::const_iterator 
    getCoordinatesPtr( positive id ) const
	{
	    return positions->begin()+(2*id);
	}

    shared_ptr<FanyArray<double> > getInternalData() const
    {
      return positions;
    }

private:

    positive nbPositions;

    shared_ptr<FanyArray<double> > positions;
};

#endif // __FPositionSet2DArbitrary_hh
