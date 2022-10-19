set( FORMAT_SRC
    reorder_float_tagged_filenames
)

include_directories( ${NetCDF_INCLUDE_DIR} )
foreach( src ${FORMAT_SRC} )
    add_executable( ${src} ${src}.cpp )
    target_link_libraries( ${src} ${Teem_LIBRARIES})
    install( TARGETS ${src} RUNTIME DESTINATION bin )
endforeach( src ${FORMAT_SRC} )

if ( VTK_FOUND )
    set( CONVERTER_SRC
        vtk2nrrd
        nrrd2vtk
        vtk2xml
    )
    foreach( src ${CONVERTER_SRC} )
        include( ${VTK_USE_FILE} )
            add_executable( ${src} ${src}.cpp )
        target_link_libraries( ${src} ${Teem_LIBRARIES} ${VTK_LIBRARIES} ${Boost_LIBRARIES} )
        install( TARGETS ${src} RUNTIME DESTINATION bin )
    endforeach()
    
    if ( NetCDF_FOUND )
        set( DLR_CONVERTER_SRC
            dlr2vtk
        )
        foreach( src ${DLR_CONVERTER_SRC} )
            include( ${VTK_USE_FILE} )
            add_executable( ${src} ${src}.cpp )
            target_link_libraries( ${src} 
                                   ${Teem_LIBRARIES} 
                                   ${VTK_LIBRARIES} 
                                   ${Boost_LIBRARIES}
                                   ${VTK_NETCDF_LIBRARY} 
                                   ${NetCDF_LIBRARIES}
            )
            install( TARGETS ${src} RUNTIME DESTINATION bin )
        endforeach()
    endif()
endif()
