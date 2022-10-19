//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FColorMap.cc,v $
// Language:  C++
// Date:      $Date: 2004/08/26 12:11:00 $
// Author:    $Author: wiebel $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FColorMap.hh"
#include "FException.hh"
#include "eassert.hh"
//---------------------------------------------------------------------------

FColor FFixedColorMap::map( double /*value*/ ) const
{
    return c;
}

//---------------------------------------------------------------------------

FColor FSegmentationColorMap::map( double value ) const
{
    value *= 5.0;

    if( value < 0.5 )
	// light gray 
	//return FColor( 0.9, 0.9, 0.9 );
	// white 
	return FColor( 1.0, 1.0, 1.0 );
    else if( value < 1.5 )
	// red
	return FColor( 1.0, 0.0 , 0.0 );
    else if( value < 2.5 )
	// orange 
	return FColor( 1.0, 0.75 , 0.0 );
    else if( value < 3.5 )
	//  green
	return FColor( 0.0, 1.0, 0.0 );
    else if( value < 4.5 )
	// light blue
	return FColor( 0.0, 1.0, 1.0 );
    else if( value < 5.5 )
	// blue 
	return FColor( 0.0, 0.0, 1.0 );
    else if( value < 6.5 )
	// yellow 
	return FColor( 1.0, 1.0, 0.0 );
    else 
	return FColor( 1.0, 1.0, 1.0 );
}
//---------------------------------------------------------------------------


FColor FDualColorMap::map( double value ) const
{
    value *= 4.0;

    if( value < 0. )
	return FColor( 0.f, 0.f, 0.f );
    else if( value < 1. )
	// blue -> light blue
	return FColor( 0.0f, value, 1.0f );
    else if( value < 2. )
	// light blue -> white
	return FColor( value-1, 1.0f, 1.0f );
    else if( value < 3. )
	// white -> yellow
	return FColor( 1.0f, 1.0f, 3.-value );
    else if( value <= 4. )
	// yellow -> red
	return FColor( 1.0f, 4.-value, 0.0f );
    else 
	return FColor( 0.f, 0.f, 0.f );
}



//---------------------------------------------------------------------------


FColor FDefaultColorMap::map( double value ) const
{
    value *= 5.0;

    if( value < 0. )
	return FColor( 0.f, 0.f, 0.f );
    else if( value < 1. )
	// blue -> light blue
	return FColor( 0.0f, value, 1.0f );
    else if( value < 2. )
	// light blue -> green
	return FColor( 0.0f, 1.0f, 2.-value );
    else if( value < 3. )
	// green -> yellow
	return FColor( value-2., 1.0f, 0.0f );
    else if( value < 4. )
	// yellow -> red
	return FColor( 1.0f, 4.-value, 0.0f );
    else if( value <= 5. )
	// red -> magenta
	return FColor( 1.0f, 0.0f, value-4. );
    else 
	return FColor( 1.0f, 0.f, 1.0f );
}


//---------------------------------------------------------------------------

FColor FHSVColorMap::map( double value ) const
{
    FColor c;
    c.setHSV( value, 1.0f, 1.0f );

    return c;
}

//---------------------------------------------------------------------------

FColor FGrayScaleMap::map( double value ) const
{
    return FColor( value, value, value );
}


//---------------------------------------------------------------------------

FColor FConstantIntensityMap::map( double value ) const
{
    static const FColor col[13] ={
        FColor( 79.f/255.f,   0.f/255.f,  99.f/255.f), // 79, 00, 99
        FColor( 93.f/255.f,  35.f/255.f, 179.f/255.f), // 93, 35,179
        FColor( 75.f/255.f,  75.f/255.f, 195.f/255.f), // 75, 75,195
        FColor( 47.f/255.f, 103.f/255.f, 163.f/255.f), // 47,103,163
        FColor(  0.f/255.f, 123.f/255.f, 123.f/255.f), //  0,123,123
        FColor(  0.f/255.f, 135.f/255.f,  39.f/255.f), //  0,135, 39
        FColor( 35.f/255.f, 139.f/255.f,   0.f/255.f), // 35,139,  0
        FColor(127.f/255.f, 139.f/255.f,   0.f/255.f), //127,139,  0
        FColor(175.f/255.f, 131.f/255.f,   0.f/255.f), //175,131,  0
        FColor(211.f/255.f, 107.f/255.f,   0.f/255.f), //211,107,  0
        FColor(247.f/255.f,  51.f/255.f,   0.f/255.f), //247, 51,  0
        FColor(255.f/255.f,   0.f/255.f,  71.f/255.f), //255,  0, 71
        FColor(255.f/255.f,   0.f/255.f, 159.f/255.f), //255,  0,159
    };

    value *= 12;
    int sel= (int)floorf(( float )value);
    if(sel>=12) return col[12];
    if(sel<0)  return col[0];
    value -= sel;
    return col[sel+1]*value + col[sel]*(1.-value);
}

//---------------------------------------------------------------------------

double FMonotonIntensityMap::omega( double value ) const
{
    return sqrt(3./2.)*value*(1.-value);
}

double FMonotonIntensityMap::z( double value ) const
{
    return sqrt(3.)*value;
}

FColor FMonotonIntensityMap::map( double value ) const
{
    FColor col;
    static const float onepsqrt3= 1.f+sqrt(3.0f);
    static const float onemsqrt3= 1.f-sqrt(3.0f);
    float wvalue = omega(value);
    float zvalue = z(value);
    col.setRed  ((1.f/2.f/sqrt(3.0f))*(onepsqrt3*wvalue*sin(frequency*value+basecolor)
            +    onemsqrt3*wvalue*cos(frequency*value+basecolor)
            +   2.*zvalue));
    col.setGreen((1.f/2.f/sqrt(3.0f))*(onemsqrt3*wvalue*sin(frequency*value+basecolor)
            +    onepsqrt3*wvalue*cos(frequency*value+basecolor)
            +   2.*zvalue));
    col.setBlue ((1.f/2.f/sqrt(3.0f))*(-2.f*wvalue*sin(frequency*value+basecolor)
                -2.f*wvalue*cos(frequency*value+basecolor)
            +   2.f*zvalue));
    return col;
}

//---------------------------------------------------------------------------

#if USE_UGLY_ADAPTIVE_COLORMAP
FAdaptiveColorMap::FAdaptiveColorMap()
{
    adaptive = false;
    split = false;
    setColorMaxMin = false;
    setMax=1;
    setMin=0;
}

//---------------------------------------------------------------------------

FColor FAdaptiveColorMap::map( double value ) const
{
      if (tvals[4]==tvals[0])
      	  return FColor(0.0, 0.0, 0.0);
      if (value < tvals[0])
          return FColor(0.0, 0.0, 1.0);
      if (value > tvals[4])
          return FColor(1.0, 0.0, 0.0);

      if ( !split )
      {
	  if( value<tvals[1] )
	      // blue -> light blue
	      return FColor( 0.0, (value-tvals[0])/(tvals[1]-tvals[0]), 1.0 );
	  else
	  {
	      if( value<tvals[2] )
		  // light blue -> green
		      return FColor( 0.0, 1.0,
				   (tvals[2]-value)/(tvals[2]-tvals[1] ));
	      else
		  if( value<tvals[3] )
		      // green -> yellow
		      return FColor( (value-tvals[2])/(tvals[3]-tvals[2]), 1.0, 0.0 );
		  else
		      // yellow -> red
		      return FColor( 1.0, (tvals[4]-value)/(tvals[4]-tvals[3]), 0.0 );
	  }    
      }
      else
      {
	  if ( value<tvals[1] )
	      // blue -> light blue
	      return FColor( 0.0, (value-tvals[0])/(tvals[1]-tvals[0]), 1.0 );
	  else
	  {
	      if( value<tvals[2] )
		  // light blue -> green
		      return FColor( 0.0, 1.0,
				   (tvals[2]-value)/(tvals[2]-tvals[1]) );
	      else
		  if( value<tvals[3] )
		      // yellow -> magenta
			  return FColor( 1.0, (tvals[3]-value)/(tvals[3]-tvals[2]),
				       (value-tvals[2])/(tvals[3]-tvals[2]) );
		  else
		      // magenta -> red
			  return FColor( 1.0, 0.0, (tvals[4]-value)/(tvals[4]-tvals[3]) );
	  }    
      }
}

//---------------------------------------------------------------------------
bool isLess (double p1, double p2)
{
  if (p1 < p2)
    return true;
  else
    return false;
}

void FAdaptiveColorMap::calculateMap( std::vector<double> values)
{
  /**
   * \todo FIXME:this sort uses a self defined isLess because
   * with our current compiler this is much faster
   */
  sort(values.begin(), values.end(),isLess);

    double min, max;
    if(setColorMaxMin){
      min=setMin;
      max=setMax;
    }
    else{
      min = values[0];
      max = values.back();
    }

    if(split && ( min>0 || max<0 ) )
      split=false;

    unsigned int zeroId=0;
    if ( split )
    {
      // locate zero crossing in list
      while ( values[zeroId] < 0. ) { zeroId++; }
    }
    
    if ( adaptive )
    {
      if(!split)
      {
	tvals[0] = values[0];
	tvals[1] = values[values.size()/4];
	tvals[2] = values[values.size()/2];
	tvals[3] = values[values.size()*3/4];
	tvals[4] = values.back();
      }
      else
      {
	tvals[0] = values[0];
	tvals[1] = values[zeroId/2];
	tvals[2] = 0;
	tvals[3] = values[ (values.size() - zeroId)/2 + zeroId ];
	tvals[4] = values.back();
      }
    }
    else 
    {
      if(!split)
      {
	tvals[0] = min;
	tvals[1] = .25*(3*min+max);
	tvals[2] = .5*(min+max);
	tvals[3] = .25*(min+3*max);
	tvals[4] = max;
      }
      else
      {
	if(-min>max)max=-min;
	tvals[0] = -max;
	tvals[1] = .5*-max;
	tvals[2] = 0;
	tvals[3] = .5*max;
	tvals[4] = max;
      }
//	cout << "min = " << min << ", max = " << max << endl;
      }
/*    
    cout<<"COLOR *** MINIMUM : "<<values[0]<<endl
	<<"COLOR *** MAXIMUM : "<<values.back()<<endl
	<<"COLOR *** MEDIUM  : "<<sum/(double)nb<<endl
	<<"COLOR *** MEDIAN  : "<<values[nb/2]<<endl;
*/
}

void FAdaptiveColorMap::setQuartiles( const std::vector<double> &quartiles)
{
  eassert(quartiles.size()==5);
  for(int i=0; i<5; ++i)
    tvals[i] = quartiles[i];
}
#endif
//---------------------------------------------------------------------------
