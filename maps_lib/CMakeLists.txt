set( MAPS_SRC
    # adaptive_map_sampling_03.31.2010.cpp
    definitions.cpp
    period.cpp
    newton.cpp
)

include_directories( ${TOKAMAK_INCLUDE_DIR} )
add_library( xavier_maps_lib ${LIBMODE} ${MAPS_SRC} )
target_link_libraries( xavier_maps_lib ${Teem_LIBRARIES}
                       xavier_poincare
                       # ${TOKAMAK_LIBRARIES}
                       # ${GL_BUNDLE_LIBRARIES}
                      )
install( TARGETS xavier_maps_lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

set( MAPS_LIB_EXEC_SRC_NEEDS_HDF5
    # ftle_tokamak_fast
    # ftle_tokamak_incremental
)

set( MAPS_LIB_EXEC_SRC_NEEDS_GLUT
    # analysis
    # map_analyzer
    # map_to_field
    # map_to_field_adaptive
    # map_to_vectors
    # parallel_analysis
    # parallel_analysis_stdmap
    pplot
    # safety_factor
    # show_map
    # stdmap_to_vectors
    # test_newton
    # test_period_analysis
    # visualize
)

set( MAPS_LIB_EXEC_SRC
    map_ftle_cr3bp
    map_ftle
    map_ftle_stdmap
    map_symplectic4D
)

if ( ${HDF_BUNDLE_FOUND} )
    include_directories( ${HDF_BUNDLE_INCLUDE_DIR} )
    foreach( target ${MAPS_LIB_EXEC_SRC_NEEDS_HDF5} )
        add_executable( ${target} ${target}.cpp )
        target_link_libraries( ${target} ${Teem_LIBRARIES} xavier_maps_lib
            ${HDF_BUNDLE_LIBRARIES}
            ${VTK_LIBRARIES}
        )
        install( TARGETS ${target} RUNTIME DESTINATION bin )
    endforeach()
endif()

if ( ${GL_BUNDLE_FOUND} )
    include_directories( ${GL_BUNDLE_INCLUDE_DIR} )
    foreach( target ${MAPS_LIB_EXEC_SRC_NEEDS_GLUT} )
        add_executable( ${target} ${target}.cpp )
        target_link_libraries( ${target} ${Teem_LIBRARIES} xavier_maps_lib
            ${GL_BUNDLE_LIBRARIES}
            ${VTK_LIBRARIES}
        )
        install( TARGETS ${target} RUNTIME DESTINATION bin )
    endforeach()
endif()


foreach( target ${MAPS_LIB_EXEC_SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} xavier_maps_lib
        ${TBB_LIBRARY_RELEASE} ${TBB_LIBRARIES}
        ${VTK_LIBRARIES}
    )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()
