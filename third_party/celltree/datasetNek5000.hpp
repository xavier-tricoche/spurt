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

#ifndef __datasetNek5000_hpp
#define __datasetNek5000_hpp

#include "avtNek5000FileFormat.h"
#include "dataset.hpp"
#include <string>

namespace CellTree {

/* datasetNek5000: import CFD mesh and associated velocity attribute
 *                 in Nek5000 format used e.g., by Paul Fisher at
 *                 Argonne National Laboratory, Argonne, Illinois, USA
 */ 
class datasetNek5000: public dataset
{
public:

    datasetNek5000( const std::string& metafile )
    {
        load_from_file(metafile);
    }
    
    void load_from_file(const std::string& metafile) {
        fmt(metafile.c_str());
    }

   virtual mesh* read_mesh() const
   {
       return fmt.GetMesh();
   };

   virtual variable* read_scalar_variable( unsigned int timestep, const std::string& name ) const
   {
       throw std::runtime_error( "not implemented" );
   }

   virtual variable* read_vector_variable( unsigned int timestep, const std::string& name ) const
   {
       if( name != "velocity" )
           throw std::runtime_error( "unknown variable " + name );

       return fmt.GetVectorVar( timestep );
   }

   virtual double get_timestep( unsigned int num ) const
   {
       std::vector<float> t;
       fmt.GetTimes( t );

       return t[num];
   }

   virtual unsigned int get_num_timesteps() const
   {
       return fmt.GetNTimesteps();
   }

protected:

   mutable avtNek5000FileFormat fmt;
};

} // namespace CellTree

#endif // __datasetNek5000_hpp
