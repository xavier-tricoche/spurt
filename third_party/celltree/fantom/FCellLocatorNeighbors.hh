//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellLocatorNeighbors.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:04 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellLocatorNeighbors_hh
#define __FCellLocatorNeighbors_hh

#include <set>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "FCellLocator.hh"
#include "FArray.hh"
#include "FIndex.hh"
#include "FCell.hh"

class FkdTree;
class FGrid;
class FNeighborhoodData;

//using namespace boost;
using boost::shared_ptr;
//using namespace std;
using std::ostream;
using std::istream;
using std::vector;

/**
 * FCellLocatorNeighbors is a class that reaches a cell by
 * going from a start cell center to the cell which contains the
 * destination position through all trhe cells which lie in between
 * step by step, using the neighborhood info
 */

class FCellLocatorNeighbors : public FCellLocator
{
public:
   struct AccelData {
    AccelData()
    {
      lastGonePositions=0;lastGoneCells=0;
    }
    ~AccelData()
    {
     if(lastGonePositions)
       delete lastGonePositions;

     if(lastGoneCells)
       delete lastGoneCells;
    }
    mutable vector<FArray>* lastGonePositions;
    mutable vector<FIndex>* lastGoneCells;

    //const FGrid* const grid;
    //const FNeighborhoodData* neighborData;
    
    mutable FIndex lastCell;

    /// parameters for last cell

    // last cell
    mutable shared_ptr<FCell> cell; 

    mutable FPosition cellCenter;
    mutable double radiusBSphereSq;//radius of bounding sphere,squared

    void computeCellBall();
  };

   void* allocateData() const{return new AccelData;}
   void freeData(void * data) const{ delete static_cast<AccelData*>( data );}
 
  /**
   *\par Description:
   *Constructor: provides an new FCellLocatorNeighbors 
   *\pre
   * the grid can already execute the getCell method
   * and also the getNeighborData function
   *\param grid
   *grid for which celllocator is needed
   */
  FCellLocatorNeighbors( const FGrid* grid );

  /**
   *\par Description:
   *Destructor, has to be overloaded by subclasses
   */
  virtual ~FCellLocatorNeighbors();

  
  //#####################################################################
  //  USER INTERFACE : REQUESTING CELL INDICES FOR POSITIONS
  //#####################################################################
  
  /**
   *\par Description:
   *Returns the index of the (a) cell containing a given position.
   *Note that this might be impossible in some cases.
   *\pre
   *none
   *\post
   *aIndex contains the id of the found cell if any
   *\exception
   *none
   *\param
   *aPosition: the input position
   *\param
   *aIndex: the found cell index
   *\return
   *\c true if a cell was found
   *\c false if no cell was found containing the position
   */
  bool searchCell( FIndex& aIndex, const FPosition& aPosition ) const;
  bool searchCell( FIndex& aIndex, const FPosition& aPosition, AccelData*data ) const;

  /**
   *\par Description:
   *
   *\pre
   *
   *\post
   *
   *\exception
   *
   *\param tmp
   * stream to write
   *
   */
  virtual void info( ostream& tmp ) const;

  /**
   *  for visAlgo cellLocatorNeighborsTest/showkdtree
   */

  shared_ptr<const FkdTree> getKdTree() const;

#if 0
  void setFillLastGone(bool b=true) const;
  const vector<FArray> & getLastGonePositions() const;
  const vector<FIndex> & getLastGoneCells() const;
#endif
  
  virtual positive memSize() const;
  
private:

  shared_ptr<const FkdTree> tree;

protected:
  const FGrid* const grid;
  const FNeighborhoodData* neighborData;
/*
  mutable vector<FArray>* lastGonePositions;
  mutable vector<FIndex>* lastGoneCells;

  
  mutable FIndex lastCell;

  /// parameters for last cell

  // last cell
  mutable shared_ptr<FCell> cell; 

  mutable FPosition cellCenter;
  mutable double radiusBSphereSq;//radius of bounding sphere,squared
*/
  mutable AccelData mydata;
  
  //compute cellCenter+radius
  //precondition: cell is valid
  void computeCellBall() const;

public:

  /**
   *\par Description:
   *asks all cells neighboring to start cell or start pos
   * if they contain searchPos, cell layer by cell layer
   *
   *\param start
   *index of start cell or start position (dependant on startIsCell)
   *
   *\param searchPos
   *position to find cell which contains it
   *
   *\param maxlayers
   * maximum number of layers to search in ( complexity of this algorithm
   *is approximately quadratic with number of layers )
   *
   *\return
   *either index of cell which contains searchPos
   *or invalid (if maxlayers layers were reached )
   */
  FIndex isInsideNeighbors(FIndex start,const FPosition&searchPos,
			   int maxlayers,bool startIsCell=true) const;
    
};

#endif // __FCellLocatorNeighbors_hh
