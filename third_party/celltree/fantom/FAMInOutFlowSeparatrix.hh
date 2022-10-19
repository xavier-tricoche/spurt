//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMInOutFlowSeparatrix.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.6 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMInOutFlowSeparatrix_hh
#define __FAMInOutFlowSeparatrix_hh
 
#include "FAMElement.hh"
#include "FPosition.hh"
#include "FIndex.hh"
#include "FVector.hh"

#include <vector>

/** 
 * Contains separatrices for the inflow/outflow field, defined on the boundary
 * of a cellcluster or the whole grid. Therefore the relation to cells may be
 * helpful and can be stored additionally.
 */
class FAMInOutFlowSeparatrix : public FAMElement
{
public:

  //---------------------------------------------------------------
  //            Class specific type definitions
  //---------------------------------------------------------------
  /**
   *  This holds a single line segment of the separatrix,
   *  which is no more than a sequence of points defining a
   *  linestrip. The separatrix itself consists of the end points
   *  of the directed separatrix-pieces, the start point being
   *  the end point of the last piece.
   *  Viewed from the outside of the grid, to the right of the line
   *  segmentwe have inflow and to the left outflow.
   */
  typedef struct SeparatrixPieceStruct
  {
    SeparatrixPieceStruct()
	: point(3)
    { };

    /**
     * \param a
     *  The endpoint of a segment.
     * \param tan
     *  Boolean representing wether the flow is tangential from the
     *  outside of the grid or not (which means tangential from the
     *  inside).
     */
    SeparatrixPieceStruct(const FVector& a,
				      bool tan)
	:point(a), isOutTangential(tan)
    { };
    
    FVector point;
    bool isOutTangential;
  } SeparatrixPiece;
  
  /**
   * \par Description:
   *  Standard constructor. Empty object.
   */
  FAMInOutFlowSeparatrix();

  /**
   * \par Description:
   *  Destructor.
   */
  ~FAMInOutFlowSeparatrix();
  
  /**
   * \par Description:
   *  Relict from the FObject - inheritance.
   */
  virtual const FString& getClassName() const;

  /**
   * \par Description:
   *  Inherited from FAMElement().
   */
  virtual void save( const FString& fileNameBase );

  /**
   * \par Description:
   *  Inherited from FAMElement().
   */
  virtual void load( const FString& fileNameBase );

  /**
   * \par Description:
   *  Lets you retrieve the stored inflow / outflow separatrices, the
   *  std::vector of separatrix pieces being adjacent pieces of the line strip
   *  bordering one region.
   *  Each separatrix piece also contains information wether the flow touches
   *  the boundary from the inside or the outside. (consult the documentation
   *  of the FAMInOutSeparatrix class for further details)
   * \pre
   *  The separatrices have been stored.
   * \post
   *  outSeparatrices contains the stored separatrices.
   * \param outSeparatrices
   *  Vector containing the strips.
   */
  void getInOutFlowSeparatrices(std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &inSeparatrices ) const;

  /**
   * \par Description:
   *  If available, this method lets you retrieve the index of the 2Din3D cell
   *  where a separatrix piece is located as well as the separatrices
   *  themselves.
   *
   *  It is an enhanced version of one-argument "getInOutFlowSeparatrices".
   * \pre
   *  The separatrices and separatrix to cell relations have been stored.
   * \post
   *  outSeparatrices contains the stored separatrices, outSepToCellRelations
   *  contains the corresponding cell indices.
   * \param outSeparatrices
   *  Vector containing the strips.
   * \param outSepToCellRelations
   *  Vector containing cellindices corresponding to the separatrix pieces.
   *  The index of the cell containing a particular separatrix piece
   *  ( a single line segment ) can be found in outSepToCellRelations at the
   *  same adress in the 2-dimensional array.
   */
  void getInOutFlowSeparatrices(std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &inSeparatrices,
                                std::vector< std::vector< FIndex > > &inSepToCellRelations) const;

  /**
   * \par Description:
   *  Adds a single strip of separatrix pieces to the stored ones. Know
   *  that you must not add separatrices with AND without cell relations
   *  at the same time since that would leave the stored information in
   *  an undefined state.
   * \pre
   *  none.
   * \post
   *  One more separatrix strip has been added to the stored collection.
   *  If there have been separatrix to cell relations stored before, the
   *  state is undefined.
   * \param inSeparatrix
   *  Vector containing the strip.
   */
  void addInOutFlowSeparatrix(const std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> &inSeparatrix);

  /**
   * \par Description:
   *  Adds a single strip of separatrix pieces and the corresponding relations
   *  to the including cells to the stored ones. Know that you must not add
   *  separatrices with AND without cell relations
   *  at the same time since that would leave the stored information in
   *  an undefined state.
   * \pre
   *  none.
   * \post
   *  One more separatrix strip has been added to the stored collection.
   *  The cell relations of the strip have been added to the stored collection.
   *  If there have been separatrices without cell relations stored before, the
   *  state is undefined.
   * \param inSeparatrix
   *  Vector containing the strip.
   * \param inSepToCellRelations
   *  The corresponding relations to cells. See the two-argument
   *  getInOutFlowSeparatrices().
   */
  void addInOutFlowSeparatrix(const std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> &inSeparatrix,
                              const std::vector< FIndex > &inSepToCellRelation);

  /**
   * \par Description:
   *  Stores the separatrices. See "getInOutFlowSeparatrices()".
   * \pre
   *  none.
   * \post
   *  The separatrices have been stored.
   * \param inSeparatrices
   *  Vector containing the strips.
   */
  void setInOutFlowSeparatrices(const std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &inSeparatrices);

  /**
   * \par Description:
   *  Stores the separatrices and their relations to the containing cells.
   *  See "getInOutFlowSeparatrices()" with two arguments.
   * \pre
   *  none.
   * \post
   *  The separatrices have been stored.
   * \param inSeparatrices
   *  Vector containing the strips.
   * \param inSepToCellRelations
   *  Vector containing cellindices corresponding to the separatrix pieces.
   *  The index of the cell containing a particular separatrix piece
   *  ( a single line segment ) can be found in outSepToCellRelations at the
   *  same adress in the 2-dimensional array.
   */
  void setInOutFlowSeparatrices(const std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &inSeparatrices,
                                const std::vector< std::vector< FIndex > > &inSepToCellRelations);

  /**
   * \par Description:
   *  Tells you if the information is stored on which cell each separatrix
   *  piece lies in.
   * \pre
   *  none.
   * \post
   *  none
   * \return
   *  True if the information is available.
   */
  bool hasRelations( void ) const;
  
private:
  std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > separatrices;  //!  all our separatrices...
  
  std::vector< std::vector<FIndex> > sepToCellRelations;
  //! to be interpreted as this:
  //! separatrices[][i] to separatrices[][i+1] iss located on the cell
  //! sepToCellRelation[][i]. usually this array is one less in size
};

//===========================================================================

#endif // __FAMInOutFlowSeparatrix_hh
