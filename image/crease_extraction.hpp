#ifndef __CREASEEXTRACTION_HPP__
#define __CREASEEXTRACTION_HPP__

#include <ext/nvis/math/fixed_vector.hpp>
#include <ext/xavier/data/raster.hpp>
#include <ext/xavier/image/convert.hpp>
#include <ext/xavier/image/crease.hpp>
#include <ext/xavier/image/probe.hpp>
#include <ext/xavier/data/edge.hpp>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <teem/nrrd.h>

namespace xavier {
namespace crease {
std::vector< double > crease_strength;
std::vector< unsigned int > isolated;

    void extract( const Nrrd* nrrd, 
          const Raster::Grid& grid,
          double threshold,
          bool ridge,
          std::vector< nvis::vec2 >& intersections,
          std::vector< std::list< unsigned int > >& creases );

    // helper functions to navigate the cell structure
    int next( const std::pair< int, int >& choices, unsigned int cur );
    int right( unsigned int i, unsigned int j, const Raster::Grid& grid );
    int top( unsigned int i, unsigned int j, const Raster::Grid& grid );
    int left( unsigned int i, unsigned int j, const Raster::Grid& grid );
    int down( unsigned int i, unsigned int j, const Raster::Grid& grid );
    unsigned int cellid( unsigned int i, unsigned int j, 
             const Raster::Grid& grid );
  };
};
};


void xavier::crease::
extract( const Nrrd* nrrd, 
     const Raster::Grid& grid,
     double threshold,
     bool ridge,
     std::vector< nvis::vec2 >& intersections,
     std::vector< std::list< unsigned int > >& creases )
{
    using namespace nvis;
    using namespace datastructure;
    using namespace gage_interface;
    
    crease_strength.clear();
    isolated.clear();
    
    double val0, val1;
    std::vector< vec2 > x(5);
    std::vector< vec2> p(5);
    vec2 inter;
    
  intersections.clear();
  creases.clear();

  // loop over all cells
  for ( unsigned int i=0 ; i<grid.nx-1 ; i++ )
    for ( unsigned int j=0 ; j<grid.ny-1 ; j++ )
      {
    // cyclic list of vertices
    p[0] = grid(i+1,j);
    p[1] = grid(i+1,j+1);
    p[2] = grid(i,j+1);
    p[3] = grid(i,j);
    p[4] = p[0];
    
    // corresponding cyclic list of edges as a [ cell <-- edge --> cell ] 
    // relationship
    k[0] = Edge( cellid(i,j,grid), right(i,j,grid) );
    k[1] = Edge( cellid(i,j,grid), top(i,j,grid) );
    k[2] = Edge( cellid(i,j,grid), left(i,j,grid) );
    k[3] = Edge( cellid(i,j,grid), down(i,j,grid) );
        
    // loop over 4 edges
    allpointsincell.clear();
    for ( unsigned int l=0 ; l<4 ; l++ )
      { 
        std::map< Edge, unsigned int >::iterator _it = edge2point.find( k[l] );
        if ( _it!=edge2point.end() )
          {
        allpointsincell.push_back( _it->second );
        std::cout << "reusing crease point #" << allpointsincell.back()
              << std::endl;
        continue;
          }
        else if ( processed.find( k[l] ) != processed.end() ) continue;

        processed.insert( k[l] );
        
        // Nrrd raster is scaled to [0,nx-1] x [0,ny-1]
        // that assumes that grid is invariably a unit square
        x[l][0] = ( grid.nx-1 )*p[l][0];
        x[l][1] = ( grid.ny-1 )*p[l][1];
        if ( !gH.value( x[l], val0 ) ) continue;
        
        x[l+1][0] = ( grid.nx-1 )*p[l+1][0];
        x[l+1][1] = ( grid.ny-1 )*p[l+1][1];
        if ( !gH.value( x[l+1], val1 ) ) continue;

        if ( ( !ridge && ( val0>=threshold || val1>=threshold ) ) ||
         ( ridge && ( val0<=threshold || val1<=threshold ) ) )
          continue;
        
        double strength;
        if ( find_intersection( x[l], x[l+1], inter,
                    gH, ridge, strength ) )
          {
        intersections.push_back( inter );
        crease_strength.push_back( strength );
        unsigned int interid = intersections.size()-1;
        edge2point[k[l]] = interid; 
        allpointsincell.push_back( interid );
          }
      }

    if ( allpointsincell.size() == 1 )
      {
        std::cout << "found only one crease point" << std::endl;
        isolated.push_back( allpointsincell[0] );
      }
    else if ( allpointsincell.size() > 1 )
      {
        unsigned int id0, id1;
        if ( allpointsincell.size() == 2 )
          {
        id0 = allpointsincell[0];
        id1 = allpointsincell[1];
          }
        else 
          {
        std::cout << "picking best crease points out of "
              << allpointsincell.size() << ": ";

        // connect points corresponding to largest / smallest values
        // of considered scalar field
        std::vector< double > vals( allpointsincell.size() );
        for ( unsigned int k=0 ; k<allpointsincell.size() ; k++ )
          {
            unsigned int pid = allpointsincell[k];
            nvis::vec2 p = intersections[pid];
            double val;
            gH.value( p, val );
            vals[k] = val;
          }
        double max1 = vals[0];
        unsigned int max1i = 0;
        double max2 = vals[1];
        unsigned int max2i = 1;
        for ( unsigned int k=2 ; k<allpointsincell.size() ; k++ )
          {
            if ( ( ridge && vals[k]>max1 ) ||
             ( !ridge && vals[k]<max1 ) )
              {
            max1 = vals[k];
            max1i = k;
              }
            else if ( ( ridge && vals[k]<max2 ) &&
                  ( !ridge && vals[k]<max2 ) )
              {
            max2 = vals[k];
            max2i = k;
              }
          }

        std::cout << max1i << " and " << max2i << std::endl;
        id0 = allpointsincell[max1i];
        id1 = allpointsincell[max2i];
          }

        cellpoints[cellid(i,j,grid)] = std::pair< int, int >( id0, id1 );
      }
      }

  // we have:
  // cellpoints: cell --> pair of crease points
  // we want inverse mapping:
  // vertexcells: crease point --> pair of cells
  std::map< unsigned int, std::pair< int, int > > vertexcells;
  for ( std::map< unsigned int, std::pair< int, int > >::iterator 
      it = cellpoints.begin() ; it!=cellpoints.end() ; it++ )
    {
      int c = it->first; // start cell

      // corresponding crease points
      int v1, v2;
      v1 = it->second.first;
      v2 = it->second.second;
      if ( v1>=0 )
    if ( vertexcells.find(v1) != vertexcells.end() )
      vertexcells[v1].second = c;
    else 
      vertexcells[v1] = std::pair< int, int >( c, -1 );
      if ( v2>=0 )
    if ( vertexcells.find(v2) != vertexcells.end() )
      vertexcells[v2].second = c;
    else 
      vertexcells[v2] = std::pair< int, int >( c, -1 );
    }

  std::vector< bool > checked( intersections.size(), false );
  for ( std::map< unsigned int, std::pair< int, int > >::iterator 
      it = cellpoints.begin() ; it!=cellpoints.end() ; it++ )
    {
      int ccur, v1, v2;
      ccur = it->first; // current cell
      v1 = it->second.first; // current "left" vertex
      v2 = it->second.second; // current "right" vertex

      std::list< unsigned int > newlist;
      std::map< unsigned int, std::pair< int, int > >::iterator vit, cit;

      // current cell vertices have been previously included in a valley
      if ( v1<0 || checked[v1] || v2<0 ) continue;

      // otherwise loop in both directions
      newlist.push_front( v1 );
      while ( !checked[v1] )
    {
      checked[v1] = true;

      vit = vertexcells.find( v1 );
      if ( vit == vertexcells.end() ) break;

      int cnext = next( vit->second, ccur );
      if ( cnext<0 ) break;

      cit = cellpoints.find( cnext );
      if ( cit == cellpoints.end() ) break;
      ccur = cnext;

      int vnext = next( cit->second, v1 );
      if ( vnext<0 ) break;
      
      newlist.push_front( vnext );
      v1 = vnext;
    }

      ccur = it->first;
      newlist.push_back( v2 );
      while ( !checked[v2] )
    {
      checked[v2] = true;

      vit = vertexcells.find( v2 );
      if ( vit == vertexcells.end() ) break;

      int cnext = next( vit->second, ccur );
      if ( cnext<0 ) break;

      cit = cellpoints.find( cnext );
      if ( cit == cellpoints.end() ) break;
      ccur = cnext;

      int vnext = next( cit->second, v2 );
      if ( vnext<0 ) break;
      
      newlist.push_back( vnext );
      v2 = vnext;
    }

      creases.push_back( newlist );
    }
}

inline int crease::
next( const std::pair< int, int >& choices, unsigned int cur )
{
    // choices.first < 0 => choices.second < 0
    if ( choices.first != ( int )cur ) {
        return choices.first;
    } else {
        return choices.second;
    }
}

inline int crease::
right( unsigned int i, unsigned int j, const xavier::Raster::Grid& grid )
{
    if ( i<grid.nx-2 ) {
        return i+1+j*(grid.nx-1);
    } else {
        return -1;
    }
}

inline int crease::
top( unsigned int i, unsigned int j, const xavier::Raster::Grid& grid )
{
    if ( j<grid.ny-2 ) {
        return i+(j+1)*(grid.nx-1);
    } else {
        return -1;
    }
}

inline int crease::
left( unsigned int i, unsigned int, const xavier::Raster::Grid& )
{
    if ( i>0 ) {
        return i-1;
    } else {
        return -1;
    }
}

inline int crease::
down( unsigned int i, unsigned int j, const xavier::Raster::Grid& grid )
{
    if ( j>0 ) {
        return i+(j-1)*(grid.nx-1);
    } else {
        return -1;
    }
}

inline unsigned int crease::
cellid( unsigned int i, unsigned int j, const xavier::Raster::Grid& grid )
{
    return i+j*( grid.nx-1 );
}

#endif
