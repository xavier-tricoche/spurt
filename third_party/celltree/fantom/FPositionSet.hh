//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet.hh,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:01 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet_hh
#define __FPositionSet_hh

#include "FObject.hh"
#include "FBoundingBox.hh"

#include "stdAliases.hh"
#include "FException.hh"
#include "FanyArray.hh"
#include <boost/shared_ptr.hpp>

#include <iostream>

class FkdTree;
class FNeighborhoodData;
class FIndex;

struct CoordSystem
{
    CoordSystem()
	: valid( false )
	{}
	
    bool valid;
    FArray center;
    FArray e0, e1, e2;

    FArray worldc( const FArray& p )
	{      
	    try{
	    if ( !valid ) return p;

	    if ( p.size() == 2 )
		return center + p(0)*e0 + p(1)*e1;
	    else
		return center + p(0)*e0 + p(1)*e1 + p(2)*e2;
	    }
	    catch( FException& e )
	    {
		std::cerr << "something went wrong in worldc" 
			  << std::endl
			  << "p=" << p << std::endl;
		throw;
	    }
	}

    void worldc( FArray& p )
	{  
	    try{
	    if ( !valid ) return;

	    if ( p.size() == 2 )
		p = center + p(0)*e0 + p(1)*e1;
	    else 
		p = center + p(0)*e0 + p(1)*e1 + p(2)*e2;
	    }
	    catch( FException& e )
	    {
		std::cerr << "something went wrong in worldc" 
			  << std::endl
			  << "p=" << p << std::endl;
		throw;
	    }
		
	}
};


/** 
 * FPositionSet is an abstract basis class for the management of 2D or 3D
 * position sets. It provides access to each contained position as well
 * as the number of positions and the dimension of the euclidean space 
 * they are contained in.
 * Specialized FPositionSet classes must be designed to optimize the memory
 * and access costs for typical or particular position sets, like the 
 * rectilinear ones, for instance.
 *
 *   \ingroup DataSet
 */
class FPositionSet: public FObject
{
    //==== Constructors ========================================================
public:

    /** 
     *\par Description:
     *Constructor: provides an empty PositionSet.
     *\pre
     *none
     *\post
     *none
     *\exception
     *none
     *\param
     *none
     */
    FPositionSet( positive dim );

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
    virtual const FString& getClassName() const = 0;

    /**
     *\par Description:
     *gets a reference (const) to the bounding box
     *\pre
     *the bounding box has been computed already
     *\post
     *none
     *\exception
     *none
     *\param
     *none
     */
    const FBoundingBox& getBoundingBox() const;

    /** 
     *\par Description:
     *Returns the number of positions.
     *\pre
     *none
     *\post
     *none
     *\exception
     *none
     */
    virtual positive getNbPositions() const = 0;

    /** 
     *\par Description:
     *Returns the dimension of the positions belonging to this FPositionSet.
     *\pre
     *none
     *\post
     *none
     *\exception
     *none
     *\param
     *none
     */
    positive getDimension() const;

    /** 
     *\par Description:
     *Gets a position.
     *\pre
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
    virtual void getPosition( FPosition& resultPos, 
			      const FIndex& pIndex ) const = 0;

    virtual void getPosition( vector<double>& resultPos, 
			      const FIndex& pIndex ) const = 0;

    /*
     * \return
     * approximate size in bytes
     */
    virtual positive memSize() const = 0;

  


    /*
     * \return
     * a kdtree built from this positionset
     */
   boost::shared_ptr<FkdTree> getKdTree(const FNeighborhoodData*nbdat=0) const;

  /**
   *\par if posset is distributed,
   * get distribution of the positions 
   */
  virtual void getDistribution(vector<positive>& sizes,
			       vector<string> & names) const;

    
  
    /*
     * This is a hack!!!!
     */
    virtual FanyArray< double >::const_iterator 
    getCoordinatesPtr( positive /* unused: id*/ ) const
	{
	    throw FNotImplementedException( "FPositionSet::getCoordinatesPtr" );
	}
    

    mutable CoordSystem coords;
    
protected:

    friend class FCellLocator;
    friend class FGrid;

    // Dimension of the space
    positive dimension;

    // bounding box
    FBoundingBox bBox;
  
    mutable boost::shared_ptr<FkdTree> tree;

};

#endif // __FPositionSet_hh
