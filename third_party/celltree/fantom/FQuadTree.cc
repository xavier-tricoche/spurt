//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuadTree.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:11 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FQuadTree.hh"
#include <iostream>

FQuadTree::FQuadTree()
{
  maxDepth = 0;
  maxCells = 8;
  root = 0;
  
}

//---------------------------------------------------------------------------

FQuadTree::FQuadTree( const FBoundingBox& mainBBox,
		      const positive& newmaxDepth,
		      const positive& newmaxCells )
{
  dimension = 2;
  maxDepth = newmaxDepth;
  maxCells = newmaxCells;
  // root is the only node with parent == 0, depth == 0
  // dimension is 2 now
  root = new FNode( dimension,
		    maxCells,
		    0,
		    maxDepth,
		    mainBBox );
		    
}

//---------------------------------------------------------------------------

FQuadTree::FQuadTree( const FBoundingBox& mainBBox,
		      const positive& newmaxDepth,
		      const positive& newmaxCells,
		      const vector< pair< FIndex, FPosition > >& newcells )
{
  dimension = 2;
  maxDepth = newmaxDepth;
  maxCells = newmaxCells;
  // root is the only node with parent == 0, depth == 0
  // dimension is 2 now
  root = new FNode( dimension,
		    maxCells,
		    maxDepth,
		    newcells,
		    mainBBox );
}

//---------------------------------------------------------------------------

FQuadTree::FQuadTree( const FBoundingBox& mainBBox,
		      const positive& newmaxDepth,
		      const positive& newmaxCells,
		      const vector< FPosition >& newcells )
{
  dimension = 2;
  maxDepth = newmaxDepth;
  maxCells = newmaxCells;
  // root is the only node with parent == 0, depth == 0
  // dimension is 2 now
  pair< FIndex, FPosition > temp;

  vector< pair< FIndex, FPosition > > cells;

  for( unsigned int i = 0; i < newcells.size(); i++ )
    {
      temp.first = FIndex(i);
      temp.second = newcells[i] ;
      cells.push_back( temp );
    }
      
  try
    {
      root = new FNode( dimension,
			maxCells,
			maxDepth,
			cells,
			mainBBox);
    }
  catch(FException &e)
    {
      cout << e << endl;
    }
}


//---------------------------------------------------------------------------

FQuadTree::~FQuadTree()
{
  
  delete root;
}

//---------------------------------------------------------------------------

void FQuadTree::addCell( const pair< FIndex, FPosition >& cell  )
{
  try
    {
      root->addCell( cell );
    }
  catch( FException &e )
    {
      throw e;
    }
}

//---------------------------------------------------------------------------


vector<FIndex> FQuadTree::search( const FPosition& pos )
{
  vector< pair< FIndex, FPosition > > searchResult;
  searchResult = root->search( pos );
  vector<FIndex> result;

  for( unsigned int i = 0; i < searchResult.size(); i++ )
    {
      result.push_back( searchResult[i].first );
    }
  return result;
}

//---------------------------------------------------------------------------

vector< pair< FIndex, FPosition > > FQuadTree::
searchboxes( const FPosition& pos )
{
  return root->search( pos );
}

//---------------------------------------------------------------------------

void FQuadTree::print()
{
  root->print();
}

//---------------------------------------------------------------------------

unsigned int FQuadTree::maxBucketSize()
{
  unsigned int size = 0;
  root->maxBucketSize( size );
  return size;
}

//---------------------------------------------------------------------------

void FQuadTree::save( char* filename )
{
  root->save( filename );
}

//---------------------------------------------------------------------------

void FQuadTree::load( char* filename )
{
  root->load( filename );
}
