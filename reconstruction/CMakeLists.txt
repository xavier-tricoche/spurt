set( RECONSTRUCTION_SRCS
    resample
    resampleLS
    resampleRBF
    resample_celltree
)

if( VTK_FOUND AND NetCDF_FOUND )
    include( ${VTK_USE_FILE} )
    include_directories( ${CELLTREE_INCLUDE_DIR} )
    include_directories( ${STANN_INCLUDE_DIR} )
    foreach( name ${RECONSTRUCTION_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name}
                               ${Teem_LIBRARIES}
                               ${Boost_LIBRARIES}
                               ${CELLTREE_LIBRARIES}
                               ${NetCDF_LIBRARIES} 
                               ${VTK_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
    message( WARNING "Unable to build reconstruction executables" )
    if ( NOT NetCDF_FOUND )
        message( WARNING "NetCDF was not found" )
    else()
        message( WARNING "VTK was not found" )
    endif()
endif()
