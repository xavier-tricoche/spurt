//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensorField.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#include "FTensorField.hh"
#include "FException.hh"
#include "FCell.hh"

#include <cassert>
#include "eassert.hh"
#include <sstream>

using namespace std;

// --------------------------------------------------------------------------

FTensorField::FTensorField( const std::string& newname, 
			    shared_ptr<FTensorSet> newpTensorSet,
			    shared_ptr<FGrid> newpGrid,
			    bool isCellBased) 
    : hashtable(10007)
{
    assert( newpTensorSet );
    assert( newpGrid );

    if(    newpTensorSet->getNbTensors() == newpGrid->getPositionSet()->getNbPositions()
            || newpTensorSet->getNbTensors() == newpGrid->getCellDefinitions()->getNbCells())
      ;
    else{
      std::stringstream oss;
      oss << "Number of tensors in field \"" << newname << "\" unequal number of positions or number of cells.";
      THROW_EXCEPTION(FException, oss.str());
    }

    //initialization of variable
    this->cellBased=false;

    if(isCellBased)
      {
	if(newpTensorSet->getNbTensors() == newpGrid->getCellDefinitions()->getNbCells())
	  {
	    this->cellBased=isCellBased;
	  }
	else
	  {
	    THROW_EXCEPTION(FException,"Number of tensors unequal number of cells, but isCellBased is set. This is a contradiction");
	  }
      }
    else
      this->cellBased = 
	(newpTensorSet->getNbTensors() == newpGrid->getCellDefinitions()->getNbCells())
	&&(newpTensorSet->getNbTensors() != newpGrid->getCellDefinitions()->getNbPositions());
    
    name          = newname;
    pTensorSet    = newpTensorSet;
    pGrid         = newpGrid;

    // invalidate "last cell"
    lastCellID.setToInvalid();
}

// --------------------------------------------------------------------------

FTensorField::~FTensorField()
{ 
}

// --------------------------------------------------------------------------

bool FTensorField::interpolate( FTensor& outResult,
				const FPosition& position ) const
{
//#ifndef NODEBUG
//    eassert( !cellBased );
//#endif
    if( cellBased )
    {
      // on cell based to constant interpolation inside cell
      FIndex index = pGrid->searchCellIndex( position );
      if( !index.isValid() ) return false;
      pTensorSet->getTensor( outResult, index );
      return true;
    }
    else
    {
      // if point based, do default interpolation provided by
      // cell type
      shared_ptr<FCell> c = getCell( position );   

      if( c )
      {
        c->interpolate( outResult, position );
        return true;
      }  
      return false;
    }
}

// --------------------------------------------------------------------------

bool FTensorField::derivatives( FTensor& outResult,
				const FPosition& position ) const
{  
#ifndef NODEBUG
    eassert( !cellBased );
#endif
 
    shared_ptr<FCell> c = getCell( position );   
    
    if( c )
    {
	c->derivatives( outResult, position );
	return true;
    }  
    
    return false;
}  

// --------------------------------------------------------------------------

shared_ptr<FCell> FTensorField::getCell( const FIndex& cellID ) const
{
    // shortcut if same as last cell
    if( lastCellID.isValid() && cellID==lastCellID )
	return lastCell;

    // go through cell cache
    shared_ptr<FCell> result = hashtable[cellID];

    if( result==0 )
    {
        // cell not in cache
        result = pGrid->getCell( cellID );
        assert( result );

        if(!cellBased)
            result->setTensors( pTensorSet.get() );
	    
        //add cell to cache
        hashtable.add( result, cellID );
    }
	
    lastCell = result;
    lastCellID = cellID;

    return result;
}  

// --------------------------------------------------------------------------

void FTensorField::computeSingularities() const
{
    eassert (!cellBased);
    vector< FAMSingularPoint > sings;
	list< FAMSingularPoint >   tmpSings;
    list< FAMSingularPoint >::iterator iter; 

    // compute all singularities
    for( positive i=0 ; i<pGrid->getCellDefinitions()->getNbCells() ; i++ ) 
    {
	shared_ptr<FCell> tmpCell;

	tmpSings.clear();
	tmpCell = pGrid->getCell( i ); 
	tmpCell->setTensors( pTensorSet.get() );
	tmpCell->getZeros( tmpSings );

	for( iter = tmpSings.begin(); iter != tmpSings.end() ; iter++ ) 
	{
	    iter->setIncludingCellIndex( i );
	    sings.push_back(*iter);
	}
    }
    
    // complete AnalysisModuleData
    analysisModuleData.addAMObject (shared_ptr<FAMSingularPoints> (new FAMSingularPoints (sings)));
}

// --------------------------------------------------------------------------

shared_ptr<FCell> FTensorField::getCell( const FPosition& position ) const
{
    // shortcut if same as last cell
    if( lastCell!=0 && lastCell->isInside(position))
	return lastCell;

    // look for the cell index in the grid
    FIndex index = pGrid->searchCellIndex( position );

    // position not inside grid?
    if( !index.isValid () )
	return shared_ptr<FCell>();

    return getCell( index );
}

// --------------------------------------------------------------------------

const FString& FTensorField::getClassName() const
{
  static FString className("FTensorField");
  return className;
}

// --------------------------------------------------------------------------
const string& FTensorField::getName() const
{
    return name;
}

// --------------------------------------------------------------------------

void FTensorField::setName( const std::string& newname )
{
    name = newname;
}

// --------------------------------------------------------------------------

shared_ptr<FTensorSet> FTensorField::getTensorSet() const
{
    return pTensorSet;
}

// --------------------------------------------------------------------------

shared_ptr<FGrid> FTensorField::getGrid() const
{
    return pGrid;
}

// --------------------------------------------------------------------------

void FTensorField::getConvenienceInfo( FTensorFieldInfo& info ) const
{
    info.grid       = pGrid;
    info.tensorSet  = pTensorSet;

    // sanity check
    assert( info.grid );
    assert( info.tensorSet );

    info.cellDef = pGrid->getCellDefinitions();
    info.posSet  = pGrid->getPositionSet();

    // sanity check
    assert( info.cellDef );
    assert( info.posSet );

    info.nbPos       = info.posSet->getNbPositions();
    info.posDim      = info.posSet->getDimension();
    info.nbCells     = info.cellDef->getNbCells();
    info.nbTensors   = info.tensorSet->getNbTensors();
    info.tensorOrder = info.tensorSet->getOrder();
    info.tensorDim   = info.tensorSet->getDimension();

    info.tensorComp  = info.tensorOrder == 0 ? 1 : (int)pow( (double)info.tensorDim, 
							     (double)info.tensorOrder );

    assert( info.nbTensors == info.nbCells ||
	    info.nbTensors == info.nbPos );

    info.isCellBased = this->cellBased;

    info.analysis = getAnalysisModuleData();
}

// --------------------------------------------------------------------------

bool FTensorField::isCellBased() const
{
  return cellBased;
}

// --------------------------------------------------------------------------

FAnalysisModuleData* FTensorField::getAnalysisModuleData() const
{
    return &analysisModuleData;
}


positive FTensorField::memSize() const
{
    return pGrid->memSize() + pTensorSet->memSize();
}
