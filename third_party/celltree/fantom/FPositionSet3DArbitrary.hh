//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet3DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:01 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet3DArbitrary_hh
#define __FPositionSet3DArbitrary_hh

#include "FPositionSet.hh"
#include "FanyArray.hh"
#include <boost/shared_ptr.hpp>
using namespace boost;


// forward-declaration
class FGrid3DArbitrary;

/** 
    the class has the functionality its name says
*/
class FPositionSet3DArbitrary: public FPositionSet
{
    //==== Constructors ========================================================
public:

    /** 
     *\par Description:
     *Constructor
     *\pre
     *coords contains the x,y and z component for each position.
     *( pos0[x], pos0[y], pos0[z], pos1[x] ... )
     *\post 
     *none
     *\exception
     *none
     *\param
     *none
     *\param 
     */
  FPositionSet3DArbitrary( shared_ptr< FanyArray<double> >  coords, const FBoundingBox*bb=0);
  

  void setTree(shared_ptr<FkdTree> atree);

  FPositionSet3DArbitrary( vector<double> & coords );

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
    void getPosition( FPosition& resultPos, const FIndex& pIndex ) const;

    void getPosition( vector<double>& resultPos, const FIndex& pIndex ) const;

    positive getNbPositions() const;

    virtual positive memSize() const;

    virtual FanyArray< double >::const_iterator 
    getCoordinatesPtr( positive id ) const
        {
	  return positions.begin()+3*id;
	}

  //--------------------------------------------------------------------------- 
  
  void getDistribution(vector<positive>& sizes,
		       vector<string> & names) const;

  shared_ptr<FanyArray<double> > getInternalData() const
  {
    return sharedPositionPtr;
  }


private:
  
  shared_ptr< FanyArray<double> > sharedPositionPtr;
  FanyArray<double> & positions;

  //buffer vector holding one position.
  mutable vector<double> bufvector;
  
};

#endif // __FPositionSet3DArbitrary_hh
