
//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:36:59 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FGrid.hh"
#include "FPositionSet.hh"
#include "FCellLocator.hh"
#include "FCellDefinitions.hh"

// ugh, but we need this for constructSuitableGrid(), 
// which is a kludge by itself
#include "FCellDefinitions2DStructured.hh"
#include "FCellDefinitions3DStructured.hh"
#include "FCellDefinitions2DUnstructured.hh"
#include "FCellDefinitions3DUnstructured.hh"
#include "FCellDefinitions2Din3DUnstructured.hh"
#include "FCellDefinitions3DTriangulation.hh"

#include "FPositionSet2DRectilinear.hh"
#include "FPositionSet2DCurvilinear.hh"
#include "FPositionSet2DArbitrary.hh"
#include "FPositionSet3DRectilinear.hh"
#include "FPositionSet3DCurvilinear.hh"
#include "FPositionSet3DArbitrary.hh"

#include "FGrid2DRectilinear.hh"
#include "FGrid2DCurvilinear.hh"
#include "FGrid2DArbitrary.hh"
#include "FGrid3DRectilinear.hh"
#include "FGrid3DCurvilinear.hh"
#include "FGrid3DArbitrary.hh"
#include "FGrid2Din3DArbitrary.hh"
#include "FCell2D.hh"

// 1D
#include "FGrid1Din3DArbitrary.hh"
#include "FGrid1Din2DArbitrary.hh"
#include "FCellDefinitions2DLines.hh"
#include "FCellDefinitions3DLines.hh"

// 0D
#include "FGrid0Din3DArbitrary.hh"
#include "FGrid0Din2DArbitrary.hh"
#include "FCellDefinitions2DPoints.hh"
#include "FCellDefinitions3DPoints.hh"

#include "eassert.hh"
#include <iostream>

//---------------------------------------------------------------------------

FGrid::FGrid( shared_ptr<FPositionSet> posSetPtr,
	      shared_ptr<FCellDefinitions> cellDefPtr,
	      const std::string& aName )
{
    name = aName;

    eassert( posSetPtr );
    eassert( cellDefPtr );

    posSet = posSetPtr;
    cellDef = cellDefPtr;
}

//---------------------------------------------------------------------------

FGrid::~FGrid()
{
    // locator is implicitly deleted (shared ptr)
}

//---------------------------------------------------------------------------

const FString& FGrid::getClassName() const
{
  static FString name( "FGrid" );
  return name;
}

//---------------------------------------------------------------------------

const std::string& FGrid::getName() const
{
    return name;
}

//---------------------------------------------------------------------------

void FGrid::setName(const std::string& newname)
{
    name = newname;
}

//---------------------------------------------------------------------------

shared_ptr<const FCellLocator> FGrid::getCellLocator() const
{
    return locator;
}
//---------------------------------------------------------------------------

shared_ptr<FPositionSet> FGrid::getPositionSet() const
{
    return posSet;
}

//---------------------------------------------------------------------------

shared_ptr<FCellDefinitions> FGrid::getCellDefinitions() const
{
    return cellDef;
}

//---------------------------------------------------------------------------

shared_ptr<FCell> FGrid::getCell( const FIndex& cellId ) const
{
    eassert( cellId < cellDef->getNbCells() );
    
    shared_ptr<FCell> cell = cellDef->getCellTorso(cellId);
    cell->setPositions( posSet.get() );

    return cell;
}

//---------------------------------------------------------------------------

FIndex FGrid::searchCellIndex( const FPosition& aPosition ) const
{
    eassert( locator );

    FIndex aIndex;
    locator->searchCell( aIndex, aPosition );

    return aIndex;
}

//---------------------------------------------------------------------------

shared_ptr<FCell> FGrid::searchCell( const FPosition& aPosition ) const
{
    return getCell( searchCellIndex( aPosition ) );
}

//---------------------------------------------------------------------------

positive FGrid::memSize() const
{
  return posSet->memSize() + cellDef->memSize();
}

//---------------------------------------------------------------------------

shared_ptr<FGrid> FGrid::constructSuitableGrid( const std::string& name,
						shared_ptr<FCellDefinitions> cellDef,
						shared_ptr<FPositionSet> posSet )
{
    // create grid type suitable to both cell definitions and
    // position set
   
  try{
    char posSetCase;
    char cellDefCase;

    FGrid *newgrid;

    if ( shared_dynamic_cast<FPositionSet2DRectilinear>(posSet) )
	posSetCase = 0;
    else if ( shared_dynamic_cast<FPositionSet2DCurvilinear>(posSet) )
	posSetCase = 1;
    else if ( shared_dynamic_cast<FPositionSet2DArbitrary>(posSet) )
	posSetCase = 2;
    else if ( shared_dynamic_cast<FPositionSet3DRectilinear>(posSet) )
	posSetCase = 3;
    else if ( shared_dynamic_cast<FPositionSet3DCurvilinear>(posSet) )
	posSetCase = 4;
    else if ( shared_dynamic_cast<FPositionSet3DArbitrary>(posSet) )
	posSetCase = 5;
    else {
	FException e("ERROR: unknown position set type");
	throw e;
    }

    if ( shared_dynamic_cast<FCellDefinitions2DStructured>(cellDef) )
	cellDefCase = 0;
    else if ( shared_dynamic_cast<FCellDefinitions2DUnstructured>(cellDef) )
	cellDefCase = 1;
    else if ( shared_dynamic_cast<FCellDefinitions3DStructured>(cellDef) )
      // this handles FCellDefinitions3DRectilinear, too.
	cellDefCase = 2;
    else if ( shared_dynamic_cast<FCellDefinitions3DUnstructured>(cellDef) )
	cellDefCase = 3;
    else if ( shared_dynamic_cast<FCellDefinitions2Din3DUnstructured>(cellDef) )
	cellDefCase = 4;
    else if ( shared_dynamic_cast<FCellDefinitions3DTriangulation>(cellDef) )
	cellDefCase = 5;
    else if ( shared_dynamic_cast<FCellDefinitions3DLineStrips>(cellDef ) )
      cellDefCase = 6;
    else if ( shared_dynamic_cast<FCellDefinitions3DLines>(cellDef ) )
      cellDefCase = 6;
    else if ( shared_dynamic_cast<FCellDefinitions2DLines>(cellDef ) )
      cellDefCase = 7;
    else if ( shared_dynamic_cast<FCellDefinitions3DPoints>(cellDef ) )
      cellDefCase = 8;
    else if ( shared_dynamic_cast<FCellDefinitions2DPoints>(cellDef ) )
      cellDefCase = 9;
    else {
	THROW_EXCEPTION( FException ,"ERROR: unknown cell definitions module type (-)");
    }
    
    switch ( posSetCase ) 
    {
    case 0:
	// 2D rectilinear position set
	if ( cellDefCase == 0 )
	    // 2D structured cell definitions
	    newgrid =  new FGrid2DRectilinear( posSet,
					    cellDef,
					    name );
	else if ( cellDefCase == 1 )
	    // 2D unstructured cell definitions
	    newgrid =  new FGrid2DArbitrary( posSet,
					  cellDef,
					  name );
	else 
	{
	    THROW_EXCEPTION( FException , ("ERROR: incompatible position set and cell definitions (0)"));
	}
	break;

    case 1: 
	// 2D curvilinear position set
	if ( cellDefCase == 0 )
	    // 2D structured cell definitions
	    newgrid =  new FGrid2DCurvilinear( posSet,
					    cellDef,
					    name );
	else if ( cellDefCase == 1 )
	    // 2D unstructured cell definitions
	    newgrid =  new FGrid2DArbitrary( posSet,
					  cellDef,
					  name );
	else {
	    THROW_EXCEPTION( FException ,("ERROR: incompatible position set and cell definitions (1)"));
	}
	break;

    case 2:
	// 2D arbitrary positionset
	if ( cellDefCase == 1 )
	    // 2D unstructured cell definitions
	    newgrid =  new FGrid2DArbitrary( posSet,
					  cellDef,
		              name );
    else if (cellDefCase == 7 )
      newgrid = new FGrid1Din2DArbitrary( posSet, cellDef, name );
    else if (cellDefCase == 9 )
      newgrid = new FGrid0Din2DArbitrary( posSet, cellDef, name );
	else {
	    THROW_EXCEPTION( FException ,("ERROR: incompatible position set and cell definitions (2)"));
	}
	break;

    case 3: 
	// 3D rectilinear position set
	if ( cellDefCase == 2 )
	    // 3D structured cell definitions
	    newgrid =  new FGrid3DRectilinear( posSet,
					    cellDef,
					    name );
	else if ( cellDefCase == 3 )
	    // 3D unstructured cell definitions
	    newgrid =  new FGrid3DArbitrary( posSet,
					  cellDef,
					  name );
	else if ( cellDefCase == 4 )
	    // 2Din3D unstructured cell definitions
	    newgrid =  new FGrid2Din3DArbitrary( posSet,
					      cellDef,
					      name );
    else if ( cellDefCase == 6 )
        newgrid = new FGrid1Din3DArbitrary(posSet,
                          cellDef,
                          name );
    else if ( cellDefCase == 8 )
        newgrid = new FGrid0Din3DArbitrary(posSet,
                          cellDef,
                          name );
	else {
	    THROW_EXCEPTION( FException ,"ERROR: incompatible position set and cell definitions (3)");
	}
	break;

    case 4:
	// 3D curvilinear position set
	if ( cellDefCase == 2 )
	    // 3D structured cell definitions
	    newgrid =  new FGrid3DCurvilinear( posSet,
					       cellDef,
					       name );
	else if ( cellDefCase == 3 )
	    // 3D unstructured cell definitions
	    newgrid =  new FGrid3DArbitrary( posSet,
					     cellDef,
					     name );
	else if ( cellDefCase == 4 )
	    // 3D unstructured cell definitions
	    newgrid =  new FGrid2Din3DArbitrary( posSet,
						 cellDef,
						 name );
	else {
	    THROW_EXCEPTION(FException ,"ERROR: incompatible position set and cell definitions (4)");
	}
	break;

    case 5:
    {
	// 3D arbitrary positionset
	if ( cellDefCase == 3 )
	    // 3D unstructured cell definitions
	    newgrid =  new FGrid3DArbitrary( posSet,
					     cellDef,
					     name );
    else if ( cellDefCase == 2 )
        // HELP: don't know if that really works,
        // but seems to be fine afai can see (Mario)
        newgrid =  new FGrid3DArbitrary( posSet,
                        cellDef,
                        name );
	else if ( cellDefCase == 4 || cellDefCase == 5 )
	    // 3D unstructured cell definitions
	    newgrid =  new FGrid2Din3DArbitrary( posSet,
						 cellDef,
						 name );
    else if (cellDefCase == 6 )
      newgrid = new FGrid1Din3DArbitrary( posSet, cellDef, name );
    else if (cellDefCase == 8 )
      newgrid = new FGrid0Din3DArbitrary( posSet, cellDef, name );
	else {
	    THROW_EXCEPTION(FException ,"ERROR: incompatible position set and cell definitions (5)");
	}
	break;
    }
    default: 
    {
	THROW_EXCEPTION( FException, "ERROR: cases not supported");
    }
    }
    
    return shared_ptr<FGrid>(newgrid);
  }
  CATCH_N_RETHROW( FException );
}

//---------------------------------------------------------------------------

bool FGrid::getNextCell(shared_ptr<FCell> & cell , FIndex & cellId ,
			FPosition & pos, const FPosition& newPos,
			vector<FIndex> * indicesGone,
			vector<FPosition> * positionsGone) const
{

  //variables for detection of repeated iteration
  static const positive bits=3;
  static const positive noInds=1<<bits;
  //maximum number of cells repeated at once in the last noInds cells
  static const positive maxrepeat=noInds;
  
  positive repeat=0;

  FPosition dir;

  double length;
  FIndex newInd,nbr;
  vector<FIndex> newIndV;

  positive lastIndices[noInds];
  for(positive i=0;i<noInds;i++)
    lastIndices[i]=FIndex::invalid;


  //variables for cell correction for num. errors
  FArray pos2;
  FArray dir2;
  positive itNo;

  double epsilon=0;
  

  const FNeighborhoodData* nd = getCellDefinitions()->getNeighborhoodData();


  for(itNo=0;itNo<10000;itNo++){

    //    cout<<"cellId:"<<cellId<<endl;


    positive i;
#if NULL
    for(i=0;i<noInds;i++)
      if(lastIndices[i]==cellId)
	break;
    
    if(i<noInds){
      repeat++;
      epsilon+=epsilon+1e-6;	       
      cout<<repeat<<endl;
    }      
    else{ 
      repeat=0;
      epsilon=0;
      lastIndices[itNo&(noInds-1)]=cellId;
      } // itNo&(noInds-1) ^= itNo % noInds
#endif
    
    //.......................................................................
    //get number of next face which is cut by ray
    
    dir = newPos-pos;

    if(dir.normSquare()<=1e-29*pos.normSquare())
      return true;

    if(indicesGone)
      indicesGone->push_back(cellId);

    if(positionsGone)
      positionsGone->push_back(pos);
    
    if(getPositionSet()->getDimension()==2){
      FCell2D*cell2d = (FCell2D*)cell.get();      
      nbr = cell2d->edgeToGo(length,pos,dir);
    }
    else
      nbr = cell->faceToGo(length,pos,dir);

#if NULL
    //check every 4th cell if ray is inside
    if(repeat<maxrepeat/4 && ( !nbr.isValid()||(itNo&3)==2 ) )
      if(!cell->isInside(pos+dir*(length*0.5)))
	{
	  //	  cout<<"2nd check!!"<<endl;
	  //check if pos is still in actual cell,otherwise change cell

	  //.................................................................
	  //get number of next face which is cut by ray from pos2 to pos
	  pos2=cell->getBarycenter();
	  dir2=pos-pos2; 

	  //if dir2=0
	  if(dir2.normSquare()<=1e-29*pos2.normSquare())
	    return true;
	  
	  if(getPositionSet()->getDimension()==2){
	    FCell2D*cell2d = (FCell2D*)cell.get();      
	    nbr = cell2d->edgeToGo(length2,pos2,dir2);
	  }
	  else
	    nbr = cell->faceToGo(length2,pos2,dir2);

	  if(!nbr.isValid()){
	    cout<<"get next cell: can not continue"<<endl;
	    return 1;
	  }
	  
	  //.................................................................
	  //if pos was outside, correct cell number
	  if(length2 < 1+epsilon){
	    
	    
	    if(getPositionSet()->getDimension()==2){
	      nd->getCellEdgeNeighbors(cellId,nbr,newIndV);
	      
	      if(newIndV.empty()) 
		return (length>1-cell->getPrecision());
	      
	      cellId = newIndV[0];
	    }
	    else{
	      nd->getCellFaceNeighbor(cellId,nbr,newInd);
	      
	      //if border is reached,and newPos is almost at border,
	      //return true, otherwise false
	      if(!newInd.isValid())
		return (length>1-cell->getPrecision());
	      
	      cellId = newInd;
	    }

	    cell=getCell(cellId);
	    continue;
	  }
	}                     
#endif
    //.......................................................................
    //check if pos+dir is inside cell

    //if way to newpos is shorter than way to cell border
    //(cell extended by epsilon)
    if(length>=1-epsilon) {
      if(length<1)
	pos+=length*dir;
      else
	pos=newPos;

      return true;
    }
    
    //.......................................................................
    //get next cell to iterate

    //pos is set to the point where
    //the last cell border was passed
    pos += length * dir;       
    
    if(!nbr.isValid()){
      cout<<"get next cell: can not continue"<<endl;
      return 0;
    }
    if(getPositionSet()->getDimension()==2){
      nd->getCellEdgeNeighbors(cellId,nbr,newIndV);
      
      if(newIndV.empty()) 
	return (length>1-cell->getPrecision());
      
      cellId = newIndV[0];
    }
    else{
      nd->getCellFaceNeighbor(cellId,nbr,newInd);
      
      //if border is reached,and newPos is almost at border,
      //return true, otherwise false
      if(!newInd.isValid())
	return (length>1-cell->getPrecision());
      
      cellId = newInd;
    }
    

    //.......................................................................

    cell = getCell(cellId);

    //.......................................................................
    //detection of repeated iteration of a cell
    
    for(i=0;i<noInds;i++)
      if(lastIndices[i]==cellId)
	break;
    
    if(i<noInds){
      repeat++;
      if(repeat>=maxrepeat){
	cout<<"cell number "<<cellId<<" iterated more than"
	    <<maxrepeat<<" times!!"<<endl
	    <<"assuming good cell found"<<endl;
	
	return true;
      }      
    }
    else
      { repeat=0;lastIndices[itNo&(noInds-1)]=cellId;} // itNo&(noInds-1) ^= itNo % noInds
    
  }

  
  cout<<"too many iterations!! (more than "<<itNo<<")"<<endl;
  cout<<"assuming outside of grid"<<endl;

  return false;
}
