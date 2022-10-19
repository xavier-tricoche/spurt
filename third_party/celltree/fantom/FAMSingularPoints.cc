//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPoints.hh,v $
// Language:  C++
// Date:      $Date:  $
// Author:    $Author: jaenicke $
// Version:   $Revision:  $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularPoints.hh"

#include <fstream>
#include <iostream>
using namespace std;

FAMSingularPoints::FAMSingularPoints (std::vector<FAMSingularPoint>& inPoints) : points (inPoints)
{
}

FAMSingularPoints::FAMSingularPoints ()
{
}


FAMSingularPoints::~FAMSingularPoints ()
{
}


const FString& FAMSingularPoints::getClassName () const
{
  static FString className("FAMSingularPoints");
  return className;
}


void FAMSingularPoints::getSingularPoints (std::vector<FAMSingularPoint>& outPoints) const
{
  outPoints =  points;
}


void FAMSingularPoints::setSingularPoints (std::vector<FAMSingularPoint>& inSingularPoints)
{
  points = inSingularPoints;
}
//---------------------------------------------------------------------------

void FAMSingularPoints::save(const FString& fileNameBase) {
    ofstream out(fileNameBase.c_str());
    if(!out) return;

    for( unsigned int k=0; k<points.size(); ++k )
    {
      FIndex c;
      FAMSingularPoint::singularityType type;
      FAMSingularPoint::singularityNature nature;
      
      points[k].getIncludingCellIndex( c );
      points[k].getType( type );
      points[k].getNature( nature );
      
      FPosition pos;
      points[k].getPosition( pos );
      
      vector<FVector> evec;
      vector< complex<double> >  eval;
      
      points[k].getEigenvalues( eval );
      points[k].getEigenvectors( evec );
      
      out.precision( 20 );
      
      out << c << '\t' << pos << " :" << nature << ", "
	  << type << ", ev=[";
      
      for( unsigned int l=0; l<3; ++l )
	out << eval[l] << (l==2 ? "], " : ", ");
      
      complex<double> det = eval[0] * eval[1] * eval[2];
      
      out << "sig=" << (det.real()<0 ? -1 : 1) <<' '<<points[k].getOrder() ;
      out << evec[0] << evec[1] << evec[2] <<"\n";
    }
    
    out.close();
}

//---------------------------------------------------------------------------

void FAMSingularPoints::load(const FString& fileNameBase) {
    ifstream in(fileNameBase.c_str());
    if(!in) return;

    unsigned int i=0, cellIndex, tmpOrder;
    int sig;
    char dummy[255];
    string tmpStr;
    
    

    FAMSingularPoint tmpSingPoint;
    vector<FAMSingularPoint> tmpPoints;
    tmpPoints.clear();
    
    while(!in.eof())
    {
      // ignore comments at beginning of line
        char ch;
	string chType;
	chType.resize(9,' ');
        in >> ch;
        if(in.eof()) break;
        bool finished = false;
        while(ch == '#' && !finished)
        {
//          std::cout << "comment" << std::endl;
          in.getline(dummy, 255 );
          if(in.eof()) {finished=true; break;}
          in >> ch;
        }
        if(finished || in.eof()) break;
        in.putback( ch );
//        in >> tmp; if(in.eof()) break; // a dummy value for b values in some files
        FPosition d(3);
	FVector evalVec(6);
	vector<FVector> evecVec;
	evecVec.resize(3,FVector(3));
        in >> cellIndex >> ch >> d[0] >> ch >> d[1] >> ch >> d[2];
	do
	{
	  in >> ch;
	  //cout <<"---"<< tmpStr<<endl;
	}
	while(ch!=',');
	for(positive j=0;j<9;j++)
	  in>>chType[j];
	do
	{
	  in >> ch;
	  //cout <<"---"<< tmpStr<<endl;
	}
	while(ch!='(');
	in>>evalVec[0]>>ch>>evalVec[1]
	  >>ch>>ch>>ch
	  >>evalVec[2]>>ch>>evalVec[3]
	  >>ch>>ch>>ch
	  >>evalVec[4]>>ch>>evalVec[5]
	  >>ch>>ch>>ch
	  >>ch>>ch>>ch>>ch
	  >>sig
	  >>tmpOrder
	  >>ch
	  >>evecVec[0][0]>>ch
	  >>evecVec[0][1]>>ch
	  >>evecVec[0][2]
	  >>ch>>ch
	  >>evecVec[1][0]>>ch
	  >>evecVec[1][1]>>ch
	  >>evecVec[1][2]
	  >>ch>>ch
	  >>evecVec[2][0]>>ch
	  >>evecVec[2][1]>>ch
	  >>evecVec[2][2];
	  ;







	vector< complex < double > > evals;
	evals.resize(3);
	evals[0]=complex<double>(evalVec[0],evalVec[1]);
	evals[1]=complex<double>(evalVec[2],evalVec[3]);
	evals[2]=complex<double>(evalVec[4],evalVec[5]);


	cout<<"======== "<<cellIndex<<"\t"<<d<<"\t"<<tmpStr<<"\t"<<evals[0]<<evals[1]<<evals[2]<<"\t"<<sig<<"\t"<<tmpOrder<<endl;

	tmpSingPoint.setIncludingCellIndex(cellIndex);
	tmpSingPoint.setPosition(d);
	FAMSingularPoint::singularityType type;
       
	if(chType=="SPIRALFAM")type=FAMSingularPoint::SPIRAL_SADDLE_3D;
	else if(chType=="REPELLING")type=FAMSingularPoint::REPELL_FOCUS_3D;
	else if(chType=="FAMSingul")type=FAMSingularPoint::SADDLE_3D;
	else if(chType=="ATTRACTIN")type=FAMSingularPoint::ATTRACT_FOCUS_3D;
	else FException("Loading of a FAMSingularPoint type is not supported (yet?).");


	cout<<type<<endl;
	tmpSingPoint.setType(type);
	//tmpSingPoint.setNature();
	tmpSingPoint.setOrder(tmpOrder);
	//tmpSingPoint.setLinearNature();
	//tmpSingPoint.setDegenerateType();
	tmpSingPoint.setEigenvectors(evecVec);
	tmpSingPoint.setEigenvalues(evals);
	tmpSingPoint.getType(type);
	cout<<type<<endl;
       
	in.getline(dummy, 255 );
	
        ++i;
	tmpPoints.push_back(tmpSingPoint);
    }
    points=tmpPoints;

#ifndef NODEBUG
    std::cout << "read " << i << " singularities" << std::endl;
//     std::copy( vecs.begin(), vecs.end(), std::ostream_iterator<FVector>(std::cout,"\n"));
#endif
}

//---------------------------------------------------------------------------
