//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNode.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:08 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FNode.hh"
#include "FException.hh"
#include <cmath>
#include <fstream>
#include <iostream>

// **************************************************************************
//
// CONSTRUCTORS & DESTRUCTOR
//
// **************************************************************************


// should not be called, just for completeness inserted

FNode::FNode()
{
  dimension = 0;
  maxCells = 0;
  depth = 0;
  maxDepth = 0;
  cells.clear();
  parent = 0;
  childs.clear();
}

//---------------------------------------------------------------------------

FNode::FNode( const unsigned int& newdimension, 
	      const unsigned int& maxNumberCells,
	      const unsigned int& newdepth,
	      const unsigned int& newmaxdepth,
	      const FBoundingBox& nodeBorders )
{ 
  // checking if given dimension of the node is valid
  if( ( newdimension != 2 ) && 
      ( newdimension != 3 ) )
    {
      THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
    }
  else
    {
      dimension = newdimension;
    }

  maxCells = maxNumberCells;

  depth = newdepth;

  maxDepth = newmaxdepth;

  // building a root node
  parent = 0; 

  // is the dimension of the bounding box valid?
  if( nodeBorders.getDimension() != dimension  )
    {
      THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
    }  
  else
    {
      nodeBox = nodeBorders;
    }

  // now the node is a leaf
  cells.clear();
  childs.clear();
}

//---------------------------------------------------------------------------

FNode::FNode( const unsigned int& newdimension, 
	      const unsigned int& maxNumberCells,
	      const unsigned int& newmaxdepth,
	      const vector< pair< FIndex, FPosition > >& newcells,  
	      const FBoundingBox& nodeBorders )
{
  // checking if given dimension of the node is valid
  if( ( newdimension != 2 ) && 
      ( newdimension != 3 ) )
    {
      THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
    }
  else
    {
      dimension = newdimension;
    }
  
  maxCells = maxNumberCells;

  // first node is a root node
  depth = 0;

  maxDepth = newmaxdepth;

  // additional constrain for a root node
  parent = 0;
  
  // is the dimension of the main bounding box valid?
  if( nodeBorders.getDimension() != dimension ) 
    {
      THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
    }
  else
    {
      nodeBox = nodeBorders;
      
      cells.clear();
      childs.clear();

      for( unsigned int i = 0; i < newcells.size(); i++ )
	{
	  // the bounding box of all cells inside the grid need to have
	  // the same dimension as the node
	  if( (newcells[i].second).getDimension() != dimension )
	    {
      THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
	    }
	  else
	    {
	      try
		{
		  addCell( newcells[i] );
		}
          CATCH_N_RETHROW( FException )
	    }
	}
    }
}

//---------------------------------------------------------------------------

FNode::FNode( char* filename )
{
  load( filename );
}

//---------------------------------------------------------------------------

FNode::~FNode()
{
  // destroying the node itself
  cells.clear();
  delete parent;
  childs.clear();
}


// **************************************************************************
//
// PUBLIC METHODS
//
// **************************************************************************


void FNode::addCell( const pair< FIndex, FPosition >& newcell )
{
  if( (newcell.second).getDimension() != dimension )
    {
      THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
    }
  else
    {
      if( isLeaf() ) 
	{
	  // leaf => cell can be added if space left or maxDepth is reached
	  if( ( !maximumReached() ) || 
	      ( depth >= maxDepth ) )
	    {
	      // does the bounding box of the cell intersect the one
	      // of the bucket?
	      // this check is needed because addCell is a public method,
	      // the FNode methods do this check from outside
	      if( nodeBox.isInside( newcell.second ) )
		{
		  // dimension and location of cell are ok
		  cells.push_back( newcell );
		}
	      else
		{
		  // cell doesn't intersect nodeBox, addCell was called
		  // from outside
          THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
		}
	    }
	  else 
	    {  
	      // overflow, the node is splitted now what means it looses 
	      // the leaf status
	      split();
	      // don't forget to add the cell
	      addCell( newcell );
	    }
	}
      else // no leaf
	{
	  // check newcells bounding box against the nodeBoxes of all childs
	  for( unsigned int i = 0; i < childs.size(); i++ )
	    {
	      // add cell to all buckets that intersect the cells
	      // bounding box
	      if( childs[i]->isInsideBBox( newcell ) )
		{
		  childs[i]->addCell( newcell );
		}
	    }

	}
    }
}

//---------------------------------------------------------------------------

vector< pair< FIndex, FPosition > > FNode::search( const FPosition& pos ) const
{
  // the cells and their bounding boxes are saved inside the leaves
  if( !isLeaf() )
    {
      // doing the isInside test for all childs
      for( unsigned int i = 0; i < childs.size(); i++ )
	{
	  if( childs[i]->isInsideBBox( pos ) )
	    {
	      return childs[i]->search( pos );
	    }
	}
    }
  return cells;
}

//---------------------------------------------------------------------------

void FNode::print() const
{
  if( !isLeaf() )
    {
      for( unsigned int i = 0; i < childs.size(); i++ )
	{
	  childs[i]->print();
	}
    }
  else
    {
      double xmin, xmax, ymin, ymax, zmin, zmax;
      nodeBox.getRange( xmin, xmax, 
			ymin, ymax, 
			zmin, zmax );
      cout << "Range of bounding box: " << endl
	   << "xmin, xmax = " << xmin << "," << xmax << endl
	   << "ymin, ymax = " << ymin << "," << ymax << endl
	   << "zmin, zmax = " << zmin << "," << zmax << endl
	   << "isLeaf = " << isLeaf() << endl
	   << "depth = " << depth << endl
	   << "Number of cells = " << cells.size() << endl
	   << "Number of childs = " << childs.size() << endl
	   << endl << endl;
    }

}  

//---------------------------------------------------------------------------

void FNode::maxBucketSize( unsigned int& size )
{
  if( isLeaf() )
    {
      if( cells.size() > size )
	{
	  size = cells.size();
	}
    }
  else
    {
      for( unsigned int i = 0; i < childs.size(); i++ )
	{
	  childs[i]->maxBucketSize( size );
	}
    }
}

//---------------------------------------------------------------------------

void FNode::save( char* /*filename*/ )
{
//   ofstream out;
//   out.open( filename,
// 	    ios::out );

//   if( !out )
//     {
//       THROW_DEFAULT_EXCEPTION( FInvalidFileException );
//     }
//   else
//     {
//       out << dimension << endl
// 	  << maxDepth  << endl
// 	  << maxCells  << endl;
//       dumpTree( out );
//     }
//   out.close();
}

//---------------------------------------------------------------------------

void FNode::load( char* /*filename*/ )
{
//   // watch up: for now, execute this method just for an empty tree!!!!

//   ifstream in;
//   in.open( filename,
// 	   ios::in );

//   int actDepth, cellsSize, childsSize, dim, maxC, maxD;
//   bool leaf;
//   // bounding box variables for the node
//   double xmin, xmax, ymin, ymax, zmin, zmax;
//   // bounding box variables for cells
//   double cxmin, cxmax, cymin, cymax, czmin, czmax;

//   FIndex tempIndex;
//   FBoundingBox tempBox;
//   char* currentLine = new( char );

//   // is there a root node already?
//   bool rootnode = false;

//   FNode* actNode = this;

//   int lastdepth = 0;

//   pair< FIndex, FBoundingBox > loadedcell;
//   if( !in )
//     {
//       THROW_DEFAULT_EXCEPTION( FInvalidFileException );
//     }
//   else
//     {
//       // reading global environment settings of a tree

//       // dimension
//       in.getline( currentLine, 10, '\n' );
//       dim = atoi( currentLine );
//       cout << currentLine << "\n";

//       // maximum depth
//       in.getline( currentLine, 10, '\n' );
//       maxD = atoi( currentLine );
//       cout << currentLine << "\n";

//       // maximum number of cells per bucket
//       in.getline( currentLine, 10, '\n' );
//       maxC = atoi( currentLine );
//       cout << currentLine << "\n";

//       // one execution of while reads the information of exactly one node
//       while( 1 )
// 	{
// 	  // reading depth of the node
// 	  in.getline( currentLine, 10, ',');

// 	  // is this line valid or did we reach EOF?
// 	  if( in.eof() )
// 	    break;

// 	  actDepth = atoi( currentLine );
// 	  cout << currentLine << ",";
	  
// 	  // reading if the node is a leaf
// 	  in.getline( currentLine, 10, ',');
// 	  leaf = atoi( currentLine );
// 	  cout << currentLine << ",";

// 	  // reading number of cells
// 	  in.getline( currentLine, 10, ',');
// 	  cellsSize = atoi( currentLine );
// 	  cout << currentLine << ",";

// 	  // reading number of childs
// 	  in.getline( currentLine, 10, ',');
// 	  childsSize = atoi( currentLine );
// 	  cout << currentLine << ",";

// 	  // reading bounding box of the bucket
// 	  in.getline( currentLine, 10, ',');
// 	  xmin = atof( currentLine );
// 	  cout << currentLine << ",";

// 	  in.getline( currentLine, 10, ',');
// 	  xmax = atof( currentLine );
// 	  cout << currentLine << ",";

// 	  in.getline( currentLine, 10, ',');
// 	  ymin = atof( currentLine );
// 	  cout << currentLine << ",";

// 	  in.getline( currentLine, 10, ',');
// 	  ymax = atof( currentLine );
// 	  cout << currentLine << ",";

// 	  in.getline( currentLine, 10, ',');
// 	  zmin = atof( currentLine );
// 	  cout << currentLine << ",";

// 	  in.getline( currentLine, 10, '\n');
// 	  zmax = atof( currentLine );
// 	  cout << currentLine << "\n";

// 	  if( dim == 2 )
// 	    {
// 	      tempBox.setBoundingBox( xmin, ymin,
// 				      xmax, ymax );
// 	    }
// 	  if( dim == 3 )
// 	    {
// 	      tempBox.setBoundingBox( xmin, ymin, zmin,
// 				      xmax, ymax, zmax );
// 	    }

// 	  // building the node
// 	  if( !rootnode )
// 	    {
// 	      dimension = dim;
// 	      maxCells  = maxC;
// 	      maxDepth  = maxD;
// 	      depth     = actDepth;
// 	      nodeBox   = tempBox;
// 	      cells.clear();
// 	      parent = 0;
// 	      childs.clear();
// 	      lastdepth = 0;
// 	      rootnode  = true;
// 	      cout << "root node added" << endl;
// 	    }
// 	  else
// 	    {
// 	      FNode* newNode = new FNode( dim,
// 					  maxC,
// 					  actDepth,
// 					  maxD,
// 					  tempBox );
// 	      cout << "Node created" << endl;
// 	      // node is on the RIGHT side of the last node
// 	      if( lastdepth == actDepth )
// 		{
// 		  // same parents
// 		  newNode->setParent( actNode->getParent() );
// 		}

// 	      if( lastdepth < actDepth )
// 		{
// 		  newNode->setParent( actNode );
// 		}

// 	      if( lastdepth > actDepth )
// 		{
// 		  newNode->setParent( ( actNode->getParent() )->getParent() );
// 		}
// 	      actNode = newNode;
// 	      ( actNode->getParent() )->addChild( actNode );	      
// 	      lastdepth = actDepth;
// 	      cout << "child node added" << endl;
// 	    }

// 	  // if the node is a leaf, the next lines define the cells
// 	  if( leaf )
// 	    {
// 	      for( int i = 0; i < cellsSize; i++ )
// 		{
// 		  // reading cell index
// 		  in.getline( currentLine, 10, ',');
// 		  tempIndex = atoi( currentLine );
// 		  cout << currentLine << ",";

// 		  // reading bounding box
// 		  in.getline( currentLine, 10, ',');
// 		  cxmin = atof( currentLine );
// 		  cout << currentLine << ",";

// 		  in.getline( currentLine, 10, ',');
// 		  cxmax = atof( currentLine );
// 		  cout << currentLine << ",";

// 		  in.getline( currentLine, 10, ',');
// 		  cymin = atof( currentLine );
// 		  cout << currentLine << ",";

// 		  in.getline( currentLine, 10, ',');
// 		  cymax = atof( currentLine );
// 		  cout << currentLine << ",";

// 		  in.getline( currentLine, 10, ',');
// 		  czmin = atof( currentLine );
// 		  cout << currentLine << ",";

// 		  in.getline( currentLine, 10, '\n');
// 		  czmax = atof( currentLine );
// 		  cout << currentLine << "\n";

// 		  if( dim == 2 )
// 		    {
// 		      tempBox.setBoundingBox( cxmin, cymin,
// 					      cxmax, cymax );
// 		    }
// 		  if( dim == 3 )
// 		    {
// 		      tempBox.setBoundingBox( cxmin, cymin, czmin,
// 					      cxmax, cymax, czmax );
// 		    }

// 		  loadedcell.first = FIndex(i);
// 		  loadedcell.second = tempBox;
// 		  actNode->addCell( loadedcell );
// 		  cout << "cell added" << endl;
// 		}
// 	    }
// 	}
//     }
//   in.close();
}


// **************************************************************************
//
// PRIVATE METHODS
//
// **************************************************************************


bool FNode::isInsideBBox( const FPosition& pos ) const
{
  return nodeBox.isInside( pos );
}

//---------------------------------------------------------------------------

bool FNode::isInsideBBox( const pair< FIndex, FPosition >& cell ) const
{
  return nodeBox.isInside( cell.second );
}

//---------------------------------------------------------------------------

vector< pair< FIndex, FPosition > > FNode::returnCells() const
{
  return cells;
}

//---------------------------------------------------------------------------

bool FNode::isLeaf() const
{
  if(  childs.size() == 0 ) 
    {
      return true;
    }
  return false;
}

//---------------------------------------------------------------------------

bool FNode::maximumReached() const
{
  return (cells.size() > maxCells);
}

//---------------------------------------------------------------------------

FNode* FNode::getParent() const
{
  return parent;
}

//---------------------------------------------------------------------------

void FNode::setParent( FNode* newparent )
{
  if( newparent != 0 )
    {
      parent = newparent;
    }
  else
    {
      THROW_DEFAULT_EXCEPTION( FNullPointerAssignmentException );
    }
}

//---------------------------------------------------------------------------


void FNode::split()
{
  double xmin, xmax, ymin, ymax, zmin, zmax;

  // getting the values of the actual nodes bounding box so we can
  // calculate the split boxes
  nodeBox.getRange( xmin, xmax, 
		    ymin, ymax, 
		    zmin, zmax );
  FBoundingBox newBBox;

  // the numbers of the cases are the same as the ones in the FNode.hh
  // documentation 

  // 2D CASE
  if( dimension == 2 )
    {
      for( unsigned int i = 0; i < 4; i++ )
	{
	  switch(i)
	    {
	    case 0: 
	      {
		newBBox.setBoundingBox( xmin, (ymin+ymax)/2.0,
					(xmin+xmax)/2.0, ymax );
		break;
	      }
	    case 1:
	      {
		newBBox.setBoundingBox( xmin, ymin, 
					(xmin+xmax)/2.0, (ymin+ymax)/2.0 );
		break;
	      }
	    case 2:
	      {
		newBBox.setBoundingBox( (xmin+xmax)/2.0, ymin, 
					xmax, (ymin+ymax)/2.0 );
		break;
	      }
	    case 3:
	      {
		newBBox.setBoundingBox( (xmin+xmax)/2.0, (ymin+ymax)/2.0, 
					xmax, ymax );
		break;
	      }
	    };
	  // building the new node
	  FNode* newnode = new FNode( dimension,
				      maxCells,
				      depth+1,
				      maxDepth,
				      newBBox );


	  // who's your daddy?
	  newnode->setParent( this );
	  childs.push_back( newnode );
	}
    }
  else // 3D CASE
    {
      for( unsigned int i = 0; i < 8; i++ )
	{
	  switch(i)
	    {
	    case 0:
	      {
		newBBox.setBoundingBox( xmin, 
					ymin, 
					zmin,
					(xmin+xmax)/2.0, 
					(ymin+ymax)/2.0,
					(zmin+zmax)/2.0 );
		break;
	      }
	    case 1:
	      {
		newBBox.setBoundingBox( (xmin+xmax)/2.0, 
					ymin,
					zmin,
					xmax,
					(ymin+ymax)/2.0,
					(zmin+zmax)/2.0 );
		break;
	      }
	    case 2:
	      {
		newBBox.setBoundingBox( (xmin+xmax)/2.0,
					(ymin+ymax)/2.0,
					zmin,
					xmax,
					ymax,
					(zmin+zmax)/2.0 );
		break;
	      }
	    case 3:
	      {
		newBBox.setBoundingBox( xmin,
					(ymin+ymax)/2.0,
					zmin,
					(xmin+xmax)/2.0, 
					ymax,
					(zmin+zmax)/2.0 );
		break;
	      }
	    case 4:
	      {
		newBBox.setBoundingBox( xmin, 
					ymin, 
					(zmin+zmax)/2.0,
					(xmin+xmax)/2.0, 
					(ymin+ymax)/2.0,
					zmax );
		break;
	      }
	    case 5:
	      {
		newBBox.setBoundingBox( (xmin+xmax)/2.0, 
					ymin,
					(zmin+zmax)/2.0,
					xmax,
					(ymin+ymax)/2.0,
					zmax );
		break;
	      }
	    case 6:
	      {
		newBBox.setBoundingBox( (xmin+xmax)/2.0,
					(ymin+ymax)/2.0,
					(zmin+zmax)/2.0,
					xmax,
					ymax,
					zmax );
		break;
	      }
	    case 7:
	      {
		newBBox.setBoundingBox( xmin,
					(ymin+ymax)/2.0,
					(zmin+zmax)/2.0,
					(xmin+xmax)/2.0, 
					ymax,
					zmax );
		break;
	      }
	    };

	  // building the new node
	  FNode* newnode = new FNode( dimension,
				      maxCells,
				      depth+1,
				      maxDepth,
				      newBBox );
	  
	  // who's your daddy?
	  newnode->setParent( this );
	  childs.push_back( newnode );
	}
    }

  // now we have to distribute the cells again
  vector< pair< FIndex, FPosition > > thiscells = returnCells();
  cells.clear();

  for( unsigned int i = 0; i < thiscells.size(); i++ )
    {
      addCell( thiscells[i] );
    }
}

//---------------------------------------------------------------------------

void FNode::dumpTree( ofstream& out )
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  nodeBox.getRange( xmin, xmax, 
		    ymin, ymax, 
		    zmin, zmax );
  
  out << depth << "," 
      << isLeaf() << ","
      << cells.size() << ","
      << childs.size() << ","
      << xmin << ","
      << xmax << ","
      << ymin << ","
      << ymax << ","
      << zmin << ","
      << zmax << endl;

//   if( isLeaf() )
//     {
//       for( unsigned int i = 0; i < cells.size(); i++ )
// 	{
// 	  (cells[i].second).getRange( xmin, xmax, 
// 				   ymin, ymax, 
// 				   zmin, zmax );
// 	  out << cells[i].first << ","
// 	      << xmin << ","
// 	      << xmax << ","
// 	      << ymin << ","
// 	      << ymax << ","
// 	      << zmin << ","
// 	      << zmax << endl;
// 	}
//     }
//   else
//     {
//       for( unsigned int j = 0; j < childs.size(); j++ )
// 	{
// 	  childs[j]->dumpTree( out );
// 	}
//     }
}

//---------------------------------------------------------------------------

void FNode::addChild( FNode* newchild )
{
  childs.push_back( newchild );
}
