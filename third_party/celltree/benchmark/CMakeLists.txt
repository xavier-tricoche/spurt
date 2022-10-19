include_directories( ../fantom )

# native benchmark utility
add_executable( nbench 
                nbench.cpp
                memstat.cpp )
target_link_libraries( nbench celltree )

# fantom benchmark utility
# add_executable( fbench 
#                 fbench.cpp 
#                 memstat.cpp
#                 FDefaultInterpolator.cpp
#                 FCellTreeInterpolator.cpp )
# target_link_libraries( fbench celltree fantom )

# vtk benchmark utility
# set( VTK_FIND_VERSION 5.6 )
find_package( VTK )

if( VTK_FOUND )

    include_directories( ${VTK_INCLUDE_DIRS} )
    link_directories( ${VTK_LIBRARY_DIRS} )

    add_executable( vbench 
                     vbench.cpp
                    memstat.cpp
                    vtkInterpolator.cpp 
                    vtkCellTree.cpp 
                    vtkDLRReader.cpp )
    target_link_libraries( vbench vtkIO vtkCommon vtkFiltering vtkGraphics vtkNetCDF netcdf )

endif( VTK_FOUND )
