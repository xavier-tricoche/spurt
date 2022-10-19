#ifndef __dataset_hpp
#define __dataset_hpp

#include <string>
#include <mesh.hpp>
#include <variable.hpp>

namespace celltree {

class dataset
{
public:

    ~dataset() {}

    virtual mesh*     read_mesh() const = 0;
    
    virtual variable* read_scalar_variable( unsigned int timestep, 
                                            const std::string& name = "" ) const = 0;
                                            
    virtual variable* read_vector_variable( unsigned int timestep, 
                                            const std::string& name = "" ) const = 0;
    
    virtual double       get_timestep( unsigned int timestep ) const = 0;
    virtual unsigned int get_num_timesteps() const = 0;

    static dataset* create( const std::string& filename );
};

} // namespace celltree

#endif // __ds_reader_hpp