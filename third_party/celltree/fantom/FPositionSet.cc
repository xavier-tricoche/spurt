//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:08 $
// Author:    $Author: garth $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#include "FPositionSet.hh"
#include "FkdTree.hh"
//---------------------------------------------------------------------------

boost::shared_ptr<FkdTree>
FPositionSet::getKdTree(const FNeighborhoodData*nbdat) const
{
  if(!tree)
    tree.reset(new FkdTree( this, nbdat));

  return tree;
}


FPositionSet::FPositionSet( positive dim )
  : dimension(dim)
{
}

//---------------------------------------------------------------------------

positive FPositionSet::getDimension() const
{
    return dimension;
}

//---------------------------------------------------------------------------

const FBoundingBox& FPositionSet::getBoundingBox() const
{
    return bBox;
}
//--------------------------------------------------------------------------- 

void FPositionSet::
getDistribution(vector<positive>& /*sizes*/,
		vector<string> & /*names*/) const
{
  THROW_EXCEPTION(FNotImplementedException," this class is not distributed! ");
}
