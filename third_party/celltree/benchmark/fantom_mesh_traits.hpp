#ifndef __fantom_mesh_traits_hpp
#define __fantom_mesh_traits_hpp

#include <mesh_traits.hpp>
#include <FTensorField.hh>

template<>
struct mesh_traits< shared_ptr<FTensorField> >
{
    typedef std::size_t size_type;
    typedef double      coord_type;
    typedef double      value_type;

    static size_t memsize( shared_ptr<FTensorField> tf ) 
    {
        return tf->memSize();
    }
   
    template<typename V>
    static void extents( shared_ptr<FTensorField> tf, V* min, V* max )
    {
        FBoundingBox box = tf->getGrid()->getPositionSet()->getBoundingBox();  

        std::vector<double> range(6);
        box.getRange( range );

        min[0] = range[0];
        max[0] = range[1];
        min[1] = range[2];
        max[1] = range[3];
        min[2] = range[4];
        max[2] = range[5];
    }
   
    template<typename V>
    static void minmax( shared_ptr<FTensorField> tf, size_type index, V* min, V* max )
    {
        shared_ptr<FGrid> grid = tf->getGrid();
        
        std::vector<FIndex> ind;
        grid->getCellDefinitions()->getCellVerticesIndices( index, ind );

        std::vector<double> pos;
        grid->getPositionSet()->getPosition( pos, ind[0] );
        
        for( unsigned int d=0; d<3; ++d )
            min[d] = max[d] = pos[d];

        for( unsigned int n=1; n<ind.size(); ++n )
        {
            grid->getPositionSet()->getPosition( pos, ind[n] );

            for( unsigned int d=0; d<3; ++d )
            {
                if( pos[d] < min[d] )
                    min[d] = pos[d];

                if( pos[d] > max[d] )
                    max[d] = pos[d];
            }
        }            
    }

    static size_type points_size( shared_ptr<FTensorField> tf )
    {
        return tf->getGrid()->getPositionSet()->getNbPositions();
    }

    static size_type cells_size( shared_ptr<FTensorField> tf ) 
    {
        return tf->getGrid()->getCellDefinitions()->getNbCells();
    }
};

#endif // __fantom_mesh_traits_hpp