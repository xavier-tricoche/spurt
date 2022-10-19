//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMGradients.cc,v $
// Language:  C++
// Date:      $Date: 2000/12/18 14:25:44 $
// Author:    $Author: hlawit $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 

#include "FAMGradients.hh"
#include <fstream>
#include <iostream>
#include <iterator>

using namespace std;

//--------------------------------------------------------------------------- 

FAMGradients::FAMGradients ()
{
  changed = true;
  empty = true;
}

//--------------------------------------------------------------------------- 

FAMGradients::~FAMGradients()
{}

//--------------------------------------------------------------------------- 

void FAMGradients::change (bool flag)
{
  changed = flag;
}

//--------------------------------------------------------------------------- 

bool FAMGradients::isEmpty ( void ) const
{
  return (gradients.size()==0);
}

//--------------------------------------------------------------------------- 

void FAMGradients::load( const FString& fileNameBase )
{
  if ( fileNameBase.find( "scheme" ) == fileNameBase.size()-6 )
  {
    double diffusionTime;
    ifstream in( fileNameBase.c_str() );
    if ( !in ) return;
    std::cout << "reading " << fileNameBase << std::endl;

    int nbGradients;

    in >> diffusionTime;
    in >> nbGradients;

    std::vector<FVector> vecs;
    for ( int i=0; i<nbGradients; ++i )
    {
      FVector v( 3 );
      in >> v[ 0 ] >> v[ 1 ] >> v[ 2 ];
      vecs.push_back( v );
    }
  }
  else if ( fileNameBase.find( "scheme1" ) == fileNameBase.size()-7 )
  {
    char dummy[ 256 ];
    ifstream in( fileNameBase.c_str() );
    if ( !in ) return;
    std::cout << "reading " << fileNameBase << std::endl;

    in.getline( dummy, 255 );
    dummy[ 256 ] = '\0'; // just to be sure
    std::cout << "first line of file: " << dummy << std::endl;
    
    std::vector<FVector> vecs;
    while ( in )
    {
      FVector v( 3 );
      double G;
      double DELTA;
      double delta;
      double TE;
      in >> v[ 0 ] >> v[ 1 ] >> v[ 2 ];
      in >> G >> DELTA >> delta >> TE;
      if ( in )
        vecs.push_back( v );
    }
  }
  else
  {
    ifstream in(fileNameBase.c_str());
    std::cout << "reading " << fileNameBase << std::endl;
    if(!in) return;
    double tmp;
    unsigned int i=0;
    vector<double> bv;
    vector<FVector> vecs;
    char dummy[255];
    while(!in.eof())
    {
      if(!in){ std::cout << "break0: maybe there are ',' signs in the file?" << std::endl; break;}

      tmp = 1.;
      // ignore comments at beginning of line
        char ch;
        in >> ch;
        if(in.eof()) {std::cout << "break1" << std::endl; break;}
        bool finished = false;
        while(ch == '#' && !finished)
        {
//          std::cout << "comment" << std::endl;
          in.getline(dummy, 255 );
          if(in.eof()) {finished=true; std::cout << "break2" << std::endl; break;}
          in >> ch;
        }
        if(finished || in.eof()) { std::cout << "break3" << std::endl; break;}
        in.putback( ch );
//        in >> tmp; if(in.eof()) break; // a dummy value for b values in some files
        FVector d(3);
        in >> d[0] >> d[1] >> d[2];
        bv.push_back(tmp);
        vecs.push_back(d);
        ++i;
    }
#if 0
    // this part is deprecated as b0 vector may be anywhere in the file

    
    // in my format, first vector has to be zero
    // because b0 field is stored there
    if(vecs[0] != FVector(0.,0.,0.))
    {
      int size = vecs.size();
      FVector last = FVector(0,0,0);
      FVector tmp(3);
      for(int i=0; i< size; ++i)
      {
        tmp=vecs[i];
        vecs[i] = last;
        last = tmp;
      }
      vecs.push_back( last );
    }
#endif
    setGradients(vecs);
    empty = false;
#ifndef NODEBUG
    std::cout << "read " << i << " gradients" << std::endl;
    std::copy( vecs.begin(), vecs.end(), std::ostream_iterator<FVector>(std::cout,"\n"));
#endif
  }
}

void FAMGradients::save( const FString& fileNameBase )
{
    ofstream out(fileNameBase.c_str());
    if(!out) return;
    vector<double>::iterator bvit = bvalues.begin();
    vector<FVector>::iterator veit= gradients.begin();
    
    for(; bvit != bvalues.end() && veit != gradients.end(); ++bvit, ++veit)
    {
        out << *bvit << " " << (*veit)[0] << " " << (*veit)[1] << " " << (*veit)[2] << endl;
    }
}
