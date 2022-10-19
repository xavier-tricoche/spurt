//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FHexahedronCell.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:07 $
// Author:    $Author: garth $
// Version:   $Revision: 1.15 $
//
//--------------------------------------------------------------------------- 

#ifndef __FHexahedronCell_hh
#define __FHexahedronCell_hh

#include "FCell.hh"

//===========================================================================

/** 
 *The FHexahedronCell class defines an abstract hexahedron
 *in 3D space derived from the FCell abstract class.
 */

class FHexahedronCell : public FCell
{
protected :
  /** 
   *\par Description:
   * returns an hexahedron with indices set to invalid
   */
  FHexahedronCell();

  /** 
   *\par Description:
   * Constructor: returns a hexahedron which vertices are defined
   * by given indices
   *\param
   * vert: vector with indices of the tetrahedron vertices
   *\param
   * useVTKHexahedronEnumeration: 
   * is true if indexing of vertices should be done 
   *       like in the vtk hexahedron cell type
   * and false if indexing of vertices should be done 
   *       like in the vtk voxel cell type.
   * the default value is true.
   */
  FHexahedronCell(const vector<FIndex>&vert,bool useVTKHexahedronEnumeration=true);
  FHexahedronCell(const FIndex*vert);

public:

  /** 
   *\par Description:
   *implementation of FCell::sizeOfCellType()
   */
  virtual positive sizeOfCellType() const;


  /** 
   *\par Description:
   *Returns zeros lying in the cell.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned zeros.
   */
  void getZeros(list<FAMSingularPoint>& result) const;


  /**
   *\par Description:
   *changes enumeration from VTKvoxel
   *to VTKHexahedronEnumeration and vice versa
   *\pre v has at least 8 elements
   *\post none
   *\param v:
   *values in vertices of a cell,enumerated either 
   * like VTKvoxel or in VTKHexahedronEnumeration 
   */

  template<class T>
  static void fromToVoxelEnum(T*v)
  {
    T xch = v[3];  v[3] = v[2]; v[2] = xch;
    xch   = v[7];  v[7] = v[6]; v[6] = xch;
  }

  template<class T>
  static void fromToVoxelEnum(vector<T>&v)
  {
    T xch = v[3];  v[3] = v[2]; v[2] = xch;
    xch   = v[7];  v[7] = v[6]; v[6] = xch;
  }

  ///see hexahedron.fig
  static const geoInfo myGeoDescription;

  /**
   *\return
   * approximate size of cell in bytes
   */
  virtual positive memSize() const=0;

protected:

  
  FIndex myIndices[8];
  FRefArray myPositions[8];
  double myPositionData[24];
  FRefTensor myTensors[8];


private:

  ///functions for getZeros:
  typedef unsigned char byte;


  void splitX(const double * const t[8],//3*8 floats 
	      //(tensor values of edges)
	      float minb[3],//3 ints
	      float maxb[3],
	      int level,
	      byte signsx, //signs of 1st,2nd, 3rd tensor comp.
	      byte signsy,
	      byte signsz, 
	      list<FAMSingularPoint>& result
	      ) const;

  void splitY(const double * const t[8],//3*8 floats 
	      //(tensor values of edges)
	      float minb[3],//3 ints
	      float maxb[3],
	      int level,
	      byte signsx, //signs of 1st,2nd, 3rd tensor comp.
	      byte signsy,
	      byte signsz, 
	      list<FAMSingularPoint>& result
	      ) const ;

  void splitZ(const double * const t[8],//3*8 floats 
	      //(tensor values of edges)
	      float minb[3],//3 ints
	      float maxb[3],
	      int level,
	      byte signsx, //signs of 1st,2nd, 3rd tensor comp.
	      byte signsy,
	      byte signsz, 
	      list<FAMSingularPoint>& result
	      ) const ;
  
protected:

  //refines singularity in [minb,maxb] by newton-raphson
  virtual void computeSing(float minb[3],float maxb[3],list<FAMSingularPoint>& result) const =0;
  
};


#endif // __FHexahedronCell_hh
