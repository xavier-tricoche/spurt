//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FDataSet.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:04 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FDataSet_hh
#define __FDataSet_hh

#include <vector>
#include <boost/shared_ptr.hpp>

#include "FIndex.hh"
#include "FArray.hh"
#include "FTensorSet.hh"
#include "FPositionSet.hh"
#include "FTimeDependentTensorField.hh"

#include "FObservable.hh"

class FCell;
class FOctree;
class FGrid;
class FTimeDependentTensorField;
class FTensorField;
class FNeighborhoodData;
class FCellDefinitions;

using namespace boost;

//===========================================================================

/**
 *The FDataSet class provides a global environment for the management of
 *geometric data (grid and cell structures) associated with numerical
 *data (namely tensor values of order 0, 1, 2). Several tensor fields as 
 *well as several times-dependent tensor fields are supported in several 
 *topological resolutions. 
 *
 * \ingroup DataSet
 */
class FDataSet : public FObservable
{
public:

  /** 
   *\par Description:
   *Constructor: provides an empty FDataSet.
   */
  FDataSet();

  /** 
   *\par Description:
   *Destructor.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  ~FDataSet();

  /** 
   *\par Description:
   *Returns the current number of tensor fields
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   * \return number of tensor fields
   */
  positive getNbTensorFields() const;

  /**
   *\par Description:
   *Returns the name of the tensorfield with a given Index.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
	 * \param tensorSetID: index of an existing tensor field
   * \return name of tensorfield
	*/
  std::string getTensorFieldName(FIndex myField) const;

  //-------------------------------------------------------------------------
  // 
  // FDataSet extension: 
  //
  //-------------------------------------------------------------------------

  FIndex addTensorField( shared_ptr<FTensorField> tensorField );

  shared_ptr<FTensorField> getTensorField( const FIndex& tensorFieldId ) const;

  void removeAllTensorFields();

  void removeTensorField( shared_ptr<const FTensorField> field );

  void setCurrentTensorfieldId( const FIndex& tensorfieldId );

  // --- begin temporary hack

  FIndex addTimeDependentTensorField( shared_ptr<FTimeDependentTensorField> newTimeDependentTensorField );
  shared_ptr<FTimeDependentTensorField> getTimeDependentTensorField( const FIndex& timeDependentTensorFieldId ) const;
  shared_ptr<FTimeDependentTensorField> getTimeDependentTensorField() const;

  void removeAllTimeDependentTensorFields();
  void removeTimeDependentTensorField( shared_ptr<const FTimeDependentTensorField> timeDependentTensorField );
  void setCurrentTimeDependentTensorFieldId( const FIndex& timeDependentTensorFieldId );

  positive getNbTimeDependentTensorFields() const;
  std::string getTimeDependentTensorFieldName(const FIndex& myDependentTensorFieldId) const;

  // --- end temporary hack

  positive memSize() const;

private:

  std::string name;

  // data structures in FDataSet
  vector< shared_ptr<FTensorField> > tensorfield;
  vector< shared_ptr<FTimeDependentTensorField>   > timeDependentTensorFields;

  FIndex currentTensorfieldId;
  FIndex currentTimeDependentTensorFieldId;
};


// the global data set instance
extern FDataSet *theDataSet;

#endif // __FDataSet_hh


/*! 
  \defgroup DataSet DataSet (FAnToM's Data Structure for Data Management)

  \brief
  This module implements FAnToM's second-generation data management structures.

  \par WARNING
  Be aware that this module consists of many more classes that are not listed here because noone has added them yet.
*/
