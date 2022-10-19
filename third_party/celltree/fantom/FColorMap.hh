//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FColorMap.hh,v $
// Language:  C++
// Date:      $Date: 2004/08/26 12:11:00 $
// Author:    $Author: wiebel $
// Version:   $Revision: 1.6 $
//
//--------------------------------------------------------------------------- 

#ifndef __FColorMap_hh
#define __FColorMap_hh

#include "FColor.hh"
#include <math.h>
#include <vector>

/*
 * A bunch of colormaps. 
 */

/**
 * Base Class for all colormaps mapping scalar values to
 * colors.
 * \param value
 * value in [0,1] which will be mapped to a color
 * \returns
 * a color dependent on value and the type of the colormap
 */
class FColorMap
{
public:

    virtual ~FColorMap() {}
    
    //! Returns a color object for input double value (between 0 and 1)
    virtual FColor map( double value ) const = 0;
};

/////////////////////////////////////////////////////////////////////
//
// implementations
//
//

/**
 * Colormap returning a single fixed color.
 */
class FFixedColorMap: public FColorMap
{
  public:
    FFixedColorMap ( const FColor c ) : c(c){};
    virtual FColor map( double value ) const;
  private:
    FColor c;
    FFixedColorMap (){} // not implemented;
};

class FDefaultColorMap: public FColorMap
{
  public:
    virtual FColor map( double value ) const;
};
class FSegmentationColorMap: public FColorMap
{
  public:
    virtual FColor map( double value ) const;
};
class FDualColorMap: public FColorMap
{
  public:
    virtual FColor map( double value ) const;
};
class FHSVColorMap: public FColorMap
{
  public:
    virtual FColor map( double value ) const;
};

class FGrayScaleMap: public FColorMap
{
  public:
    virtual FColor map( double value ) const;
};

/** 
 * Taken from Thomas Lehmann: Bildverarbeitung fuer die Medizin.
 * A Colormap with constant intensity used for infrared images
 * from temperature sensitive cameras.
 */
class FConstantIntensityMap: public FColorMap
{
  public:
    virtual FColor map( double value ) const;
};

/**
 * Pseudo-coloring with monoton increasing intensity
 * Rotation of colormap specified by omega around 
 * the axis of a cube.
 */
class FMonotonIntensityMap: public FColorMap
{
  private:
    double frequency;
    double basecolor;
  protected:
    virtual double z(double) const;
    virtual double omega(double) const;
  public:
    //! tankes the number of spins as Frequency (0..M_PI) and the base color (0..M_PI)
    FMonotonIntensityMap( double freq = 4.*M_PI, double base = 0.0 )
    {
        frequency = freq;
        basecolor = base;
    }
    virtual FColor map( double value ) const;
    void setFrequency( double freq ){frequency = freq;}
    void setBaseColor( double base ){basecolor = base;}
};

#if USE_UGLY_ADAPTIVE_COLORMAP
class FAdaptiveColorMap: public FColorMap
{
    private:
    double tvals[5];
    bool adaptive;
    bool split;
    bool setColorMaxMin;
    double setMin;
    double setMax;
    public:
    FAdaptiveColorMap();
    virtual FColor map( double value ) const;

    /** this function takes all scalar that will be used
        and calculates a suitable colormap
	MUST BE CALLED at least once BEFORE EXECUTION of map(double)
    */
    virtual void calculateMap( std::vector<double> values );
    virtual void setQuartiles( const std::vector<double> &values );
    void setAdaptive(bool _adaptive=true){adaptive = _adaptive;}
    void setSplit(bool _split=true){split= _split;}
    void setColorMinMax(bool enabled, double _min=0, double _max=1){setColorMaxMin = enabled; setMin=_min; setMax=_max;}
};
#endif
//===========================================================================
#endif // __FColorMap_hh
 
