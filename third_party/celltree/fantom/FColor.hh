//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FColor.hh,v $
// Language:  C++
// Date:      $Date: 2004/08/25 16:42:52 $
// Author:    $Author: hlawit $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#ifndef __FColor_hh
#define __FColor_hh

#include "FObject.hh"
#include <iosfwd>

//===========================================================================

/** The FColor class represents a color attribute for graphical primitives.
 *  Basically, it stores the RGBA information and provides methods for
 *  changing this color information. 
 */
class FColor: public FObject
{
public:

  /** default constructor
   */
  FColor();

  
  /** copy constructor
   *  \param
   *  color FColor to be copied.
   */
  FColor( const FColor& color );


  /** Constructs FColor object with given RGBA values.
   *  \param
   *  red red value in [0,1]
   *  \param
   *  blue blue value in [0,1]
   *  \param
   *  green green value in [0,1]
   *  \param
   *  alpha alpha value in [0,1]
   *  \exception
   *  FException is raised when one of the RGBA values is not in the
   *  interval [0,1].   
   */
  FColor(float red, float green, float blue, float alpha = 1.0 ); 


  /** destructor
   */
  virtual ~FColor();
  
  /** assignment operator.
   *  \param
   *  color FColor to be assigned
   */
  FColor operator= (FColor color);


  /** Returns the class name as string.
   *  \return class name   
   */
  virtual const FString& getClassName() const;


  /** Returns the contribution of red light to the current color.
   *  \return
   *  intensity of red in [0,1].
   */
  float getRed() const;


  /** Returns the contribution of green light to the current color.
   *  \return
   *  intensity of green in [0,1].
   */
  float getGreen() const;


  /** Returns the contribution of blue light to the current color.
   *  \return
   *  intensity of blue in [0,1].
   */
  float getBlue() const;


  /** Returns the alpha value (opacity).
   *  \return
   *  opacity alpha in [0,1].
   */  
  float getAlpha() const;

  /** Returns the luminance computed by
   * 0.3*red + 0.59*green +0.11*blue.
   */
  float getLuminance() const;

  /** Sets the contribution of red light to the current color.
   *  \param 
   *  red intensity of red in [0,1].
   *  \exception
   *  FException 
   */
  void setRed(float red);


  /** Sets the contribution of green light to the current color.
   *  \param 
   *  green intensity of green in [0,1].
   *  \exception
   *  FException is raised when the value is not in the
   *  interval [0,1].  
   */  
  void setGreen(float green);


  /** Sets the contribution of blue light to the current color.
   *  \param 
   *  blue intensity of blue [0...1].
   *  \exception
   *  FException is raised when the value is not in the
   *  interval [0,1]. 
   */    
  void setBlue(float blue);


  /** Sets the opacity alpha.
   *  \param 
   *  alpha opacity [0...1].
   *  \exception
   *  FException is raised when the value is not in the
   *  interval [0,1]. 
   */      
  void setAlpha(float alpha);


  /** Set RGB from HSV
   */
  void setHSV( float h, float s, float v );

  /** Get HSV value from color
   * h in [0,1), s and v in [0,1]
   * s == 0 => h undefined
   */
  void getHSV( float &h, float &s, float &v) const;

  /** Set RGB by HLS values
   * h,s and v in [0,1]
   */
  void setHLS( float h, float l, float s);
  
  /** Get HLS values
   * h,s and v in [0,1]
   */
  void getHLS( float &h, float &l, float &s) const;
  
  /** Set RGB
   */
  void setRGB( float r, float g, float b );

  /** Returns the RGBA values as a a pointer to a GLfloat[4] array.
   *  \return
   *  pointer to GLfloat array
   */      
  const float* getOpenGLColor() const;

  friend std::ostream& operator<< (std::ostream& os, const FColor& c);
  friend std::istream& operator>> (std::istream& is, FColor& c);


  friend FColor operator+(const FColor&, const FColor&);
  friend FColor operator*(const FColor&, double);
  friend FColor operator*(double, const FColor&);

  FColor& operator*=( float );

  float operator[]( int i ) const { return rgba[i]; }
  float& operator[]( int i ) { return rgba[i]; }
protected:

  /** Checks whether the argument is in the interval [0...1] or not.
   *  \param 
   *  value Value to be checked
   *  \return
   *  false if value is in [0...1], otherwise true.
   */
  bool checkRange( float value ) const;

private:
  // helper for setHLS
  float valueHLS(float n1, float n2, float hue) const;

  // data
  float rgba[4];
};

//===========================================================================
#endif // __FColor_hh
 
