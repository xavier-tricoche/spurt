//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FOcTree.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:08 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FOcTree.hh"

FOcTree::FOcTree()
{
  maxDepth = 0;
  maxCells = 8;
  root = 0;
  
}

//---------------------------------------------------------------------------

FOcTree::FOcTree( const FBoundingBox& mainBBox,
		      const positive& newmaxDepth,
		      const positive& newmaxCells )
{
  maxDepth = newmaxDepth;
  maxCells = newmaxCells;
  // root is the only node with parent == 0, depth == 0
  // dimension is 3 now
  dimension = 3;
  root = new FNode( dimension,
		    maxCells,
		    0,
		    maxDepth,
		    mainBBox );
		    
}

//---------------------------------------------------------------------------

FOcTree::FOcTree( const FBoundingBox& mainBBox,
		  const positive& newmaxDepth,
		  const positive& newmaxCells,
		  const vector< pair< FIndex, FPosition > >& newcells )
{
  maxDepth = newmaxDepth;
  maxCells = newmaxCells;
  // root is the only node with parent == 0, depth == 0
  // dimension is 3 now
  dimension = 3;
  root = new FNode( dimension,
		    maxCells,
		    maxDepth,
		    newcells,
		    mainBBox);
}

//---------------------------------------------------------------------------

FOcTree::FOcTree( const FBoundingBox& mainBBox,
		  const positive& newmaxDepth,
		  const positive& newmaxCells,
		  const vector< FPosition >& newcells )
{
  maxDepth = newmaxDepth;
  maxCells = newmaxCells;
  // root is the only node with parent == 0, depth == 0
  // dimension is 3 now

  pair< FIndex, FPosition > temp;
  vector< pair< FIndex, FPosition > > cells;
  for( unsigned int i = 0; i < newcells.size(); i++ )
    {
      temp.first = FIndex(i);
      temp.second = newcells[i] ;
      cells.push_back( temp );
    }
      
  dimension = 3;
  root = new FNode( dimension,
		    maxCells,
		    maxDepth,
		    cells,
		    mainBBox);
}


//---------------------------------------------------------------------------

FOcTree::~FOcTree()
{
  
  delete root;
}

//---------------------------------------------------------------------------

void FOcTree::addCell( const pair< FIndex, FPosition >& cell  )
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


vector<FIndex> FOcTree::search( const FPosition& pos )
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

vector< pair< FIndex, FPosition > > FOcTree::
searchboxes( const FPosition& pos )
{
  return root->search( pos );
}

// ---------------------------------------------------------------------------

void FOcTree::print()
{
  root->print();
}

// ---------------------------------------------------------------------------

unsigned int FOcTree::maxBucketSize()
{
  unsigned int size = 0;
  root->maxBucketSize( size );
  return size;
}

//---------------------------------------------------------------------------

void FOcTree::save( char* filename )
{
  root->save( filename );
}

//---------------------------------------------------------------------------

void FOcTree::load( char* filename )
{
  root->load( filename );
}
