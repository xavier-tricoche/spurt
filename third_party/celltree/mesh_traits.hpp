/*************************************************************************
* celltree: boundary interval hierarchy tree-based solution for point 
*           location and interpolation in unstructured meshes. This
*           technique was first proposed in:
*           "Fast, Memory-Efficient Cell location in Unstructured Grids
*            for Visualization"
*                by Christop Garth and Kenneth I Joy, VisWeek 2011.
* Author: Christoph Garth
* Copyright (C) 2011 Christoph Garth
* All rights reserved
**************************************************************************/

#ifndef __traits_hpp
#define __traits_hpp

namespace CellTree {

template<typename T>
struct variable_traits
{
};

template<typename T>
struct mesh_traits
{    
};

} // namespace CellTree

#endif // __traits_hpp
