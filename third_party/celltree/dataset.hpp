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

#ifndef __dataset_hpp
#define __dataset_hpp

#include <string>
#include "mesh.hpp"
#include "variable.hpp"

namespace CellTree {

class dataset
{
public:

    virtual ~dataset() {}

    virtual mesh*     read_mesh() const = 0;
    
    virtual variable* read_scalar_variable( unsigned int timestep, 
                                            const std::string& name = "" ) const = 0;
                                            
    virtual variable* read_vector_variable( unsigned int timestep, 
                                            const std::string& name = "" ) const = 0;
    
    virtual double       get_timestep( unsigned int timestep ) const = 0;
    virtual unsigned int get_num_timesteps() const = 0;

    virtual void load_from_file( const std::string& filename ) = 0;
};

} // namespace CellTree

#endif // __ds_reader_hpp
