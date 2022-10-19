//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FColor.cc,v $
// Language:  C++
// Date:      $Date: 2004/09/02 14:21:25 $
// Author:    $Author: hlawit $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#include "FColor.hh"
#include "FException.hh"
#include <cstdio>
#include <iostream>
#include <sstream>
#include <math.h>

using namespace std;

//---------------------------------------------------------------------------

FColor::~FColor()
{
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream& os, const FColor& c)
{
  os <<"[ RGB: "<<c.getRed()<<" , "<<c.getGreen()<<" , ";
  os <<c.getBlue()<<" , A: "<<c.getAlpha()<<" ]";
  return os;
}  

//---------------------------------------------------------------------------

istream& operator>> (istream& is, FColor& c)
{
  char tmp[256];
  float rgba[4];

  is.getline(tmp, 255);
  sscanf (tmp, "%*s %*s %f %*s %f %*s %f %*s %*s %f",
	  &rgba[0],
	  &rgba[1],
	  &rgba[2],
	  &rgba[3]);

  c.setRed (rgba[0]);
  c.setGreen (rgba[1]);
  c.setBlue (rgba[2]);
  c.setAlpha (rgba[3]);

  return is;
}  

//---------------------------------------------------------------------------

FColor operator+(const FColor& c, const FColor& c2)
{
    FColor n(c);
    n.rgba[0]+=c2.rgba[0];
    n.rgba[1]+=c2.rgba[1];
    n.rgba[2]+=c2.rgba[2];
    n.rgba[3]+=c2.rgba[3];
    return n;
}
//---------------------------------------------------------------------------

FColor operator*(const FColor& c, double scale)
{
    FColor n(c);
    n.rgba[0]*=scale;
    n.rgba[1]*=scale;
    n.rgba[2]*=scale;
    n.rgba[3]*=scale;
    return n;
}

//---------------------------------------------------------------------------

FColor operator*(double scale, const FColor& c)
{
    FColor n(c);
    n.rgba[0]*=scale;
    n.rgba[1]*=scale;
    n.rgba[2]*=scale;
    n.rgba[3]*=scale;
    return n;
}

//---------------------------------------------------------------------------

const FString& FColor::getClassName() const
{
  static const FString className("FColor");

  return className;
}
//---------------------------------------------------------------------------

FColor::FColor()
{
  rgba[0] = 1.0;
  rgba[1] = 1.0;
  rgba[2] = 1.0;
  rgba[3] = 1.0;
}

//---------------------------------------------------------------------------

FColor::FColor( const FColor& color ) : FObject()
{
  rgba[0] = color.rgba[0];
  rgba[1] = color.rgba[1];
  rgba[2] = color.rgba[2];
  rgba[3] = color.rgba[3];
}

//---------------------------------------------------------------------------

FColor::FColor(float red, float green, float blue, float alpha)

{
#ifndef NODEBUG
  if( checkRange(red) || checkRange(green) 
      || checkRange(blue) || checkRange(alpha) )
    {
      cout <<"ERROR: Illegal color value."<<endl; 
      cout <<"FColor::FColor(" << red << ", " << green << ", " << blue << ", " << alpha << ")"<<endl;
    }
#endif

  rgba[0] = red;
  rgba[1] = green;
  rgba[2] = blue;
  rgba[3] = alpha;
}

//---------------------------------------------------------------------------

FColor FColor::operator= (FColor color)
{
  rgba[0] = color.rgba[0];
  rgba[1] = color.rgba[1];
  rgba[2] = color.rgba[2];
  rgba[3] = color.rgba[3];

  return *this;
}

//---------------------------------------------------------------------------

float FColor::getRed() const
{
  return (float)rgba[0];
}

//---------------------------------------------------------------------------

float FColor::getGreen() const
{
  return (float)rgba[1];
}

//---------------------------------------------------------------------------

float FColor::getBlue() const
{
  return (float)rgba[2];
}

//---------------------------------------------------------------------------

float FColor::getAlpha() const
{
  return (float)rgba[3];
}

//---------------------------------------------------------------------------

float FColor::getLuminance() const
{
  return (float)rgba[0]*0.30+(float)rgba[1]*0.59+(float)rgba[2]*0.11;
}

//---------------------------------------------------------------------------

void FColor::setRed(float red)
{
#ifndef NODEBUG
  if( checkRange(red) )
    {
      cout <<"ERROR: Illegal color value."<<endl; 
      cout <<"FColor::setRed(float)  value is "<<red<<endl;
      if(red>1.0) red=1.0;
      else red=0.0;
    }
#endif
  
  rgba[0] = (float)red;
}

//---------------------------------------------------------------------------

void FColor::setGreen(float green)
{
#ifndef NODEBUG
  if( checkRange(green) )
    {
      cout <<"ERROR: Illegal color value."<<endl; 
      cout <<"FColor::setGreen(float)  value is "<<green<<endl;
      if(green>1.0) green=1.0;
      else green=0.0;
    }
#endif
  
  rgba[1] = (float)green;
}

//---------------------------------------------------------------------------

void FColor::setBlue(float blue)
{
#ifndef NODEBUG
  if( checkRange(blue) )
    {
      cout <<"ERROR: Illegal color value."<<endl; 
      cout <<"FColor::setBlue(float)  value is "<<blue<<endl;
      if(blue>1.0) blue=1.0;
      else blue=0.0;
    }
#endif
  
  rgba[2] = (float)blue;
}

//---------------------------------------------------------------------------

void FColor::setAlpha(float alpha)
{
#ifndef NODEBUG
  if( checkRange(alpha) )
    {
      THROW_EXCEPTION(FException, "ERROR: Illegal opacity value.");
    }
#endif
  
  rgba[3] = (float)alpha;
}

//---------------------------------------------------------------------------

void FColor::setRGB( float r, float g, float b )
{
    rgba[0] = r;
    rgba[1] = g;
    rgba[2] = b;
}

//---------------------------------------------------------------------------

void FColor::setHSV( float h, float s, float v )
{
    if( s == 0.0 )
    {
	setRGB( v, v, v );
	return;
    }

    float i = floor(6.0*h);
    float f = 6.0*h - i;
    float p = v * (1.0-s);
    float q = v * (1.0-s*f);
    float t = v * (1.0-s*(1.0-f));

    int iq = (int)i % 6;

    switch( iq )
    {
    case 0:
    case 6:
	setRGB( v, t, p );
	return;
    case 1:
	setRGB( q, v, p );
	return;
    case 2:
	setRGB( p, v, t );
	return;
    case 3:
	setRGB( p, q, v );
	return;
    case 4:
	setRGB( t, p, v );
	return;
    case 5:
	setRGB( v, p, q );
	return;
    }
    // assert(false);
}

void FColor::getHSV( float &h, float &s, float &v) const
{
    float max = getRed();
    float min = getRed();
    if( max < getGreen()) max = getGreen();
    else min = getGreen();
    if( max < getBlue()) max = getBlue();
    if( min > getBlue()) min = getBlue();

    v = max;

    float delta = max -min;
    s = (max != 0.) ? (( delta )/max) : 0.0;

    if(s == 0.0) return; // h = undefined;

    if(getRed() == max)
    {
        h=(getGreen()-getBlue())/delta;
    }
    else if(getGreen() == max)
    {
        h=2.0+(getBlue()-getRed())/delta;
    }
    else //if(b==max)
    {
        h=4.0+(getRed()-getGreen())/delta;
    }
    h/=6.0;
    if(h<0.0)
        h+=1.;
}

void FColor::getHLS(float &h, float &l, float &s) const
{
    float max = getRed();
    float min = getRed();
    if( max < getGreen()) max = getGreen();
    else min = getGreen();
    if( max < getBlue()) max = getBlue();
    if( min > getBlue()) min = getBlue();
   
    l = (max+min)/2.;
    
    if( max == min )
    {
       s = 0;
      // h undefined;
    }
    else
    {
        float delta = max - min;
        s = (l<=0.5) ? ( delta /(max+min)): (delta/(2.0-(max+min)));
        if(getRed() == max)
            h = (getGreen()-getBlue())/delta;
        else if(getGreen() == max)
        {
            h=2.0+(getBlue()-getRed())/delta;
        }
        else //if(b==max)
        {
            h= 4.0+(getRed()-getGreen())/delta;
        }
        h*=60.0;
        if(h<0.) h+=360.;
    }
}

void FColor::setHLS(float h, float l, float s)
{
    float m1, m2;

    m2 =( l<=0.5 ) ? (l*(l+s)) : (l+s-l*s);
    m1 = 2.0*l-m2;
    if(s==0) { // h is undefined
        setRed(l);
        setGreen(l);
        setBlue(l);
    }
    else 
    {
        setRed(valueHLS(m1,m2,h+120.0));
        setGreen(valueHLS(m1,m2,h));
        setBlue(valueHLS(m1,m2,h-120.0));
    }
}
inline float FColor::valueHLS(float n1, float n2, float hue) const
{
    if(hue > 360.)
    {
        hue -=360.;
    }
    else if( hue < 0.)
    {
            hue += 360;
    }

    if(hue < 60.)
        return n1+(n2-n1)*hue/60.;
    if(hue < 180.)
        return n2;
    if(hue < 240.)
        return n1+(n2-n1)*(240.-hue)/60.;
    else
        return n1;
}
    

//---------------------------------------------------------------------------

bool FColor::checkRange(float value) const
{
  if( (value >= 0.0) && (value <= 1.0) )
    return false;
  else
    return true;
}

//---------------------------------------------------------------------------

const float* FColor::getOpenGLColor() const
{
  return rgba;
}

//---------------------------------------------------------------------------


FColor& FColor::operator*=( float f )
{
  for ( int i=0; i< 4; ++i )
  {
    rgba[ i ] *= f;
  }
  return *this;
}

//---------------------------------------------------------------------------

