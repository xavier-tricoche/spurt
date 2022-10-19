set(FLOW_EXEC_SRC
    best-scale
    color_ftle
    generate_seeds
    vector_ftle
    vector_ftle_euler
    catseye
    doublegyre
    abc
    flowmap_in_3D_nrrd
    #flow_map
    flow_map_new
    generate_flow
)

set(FLOW_EXEC_NEEDS_VTK
    tangent_lines
    lavd
    trace_particles
    view_flowmap
    view_flow
    test_upsample
    parse_lavd_header
    view_particles_path
    flow_map_unstructured
)

foreach( target ${FLOW_EXEC_SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} spurt_math  spurt_image ${Boost_LIBRARIES} )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()

if( ${VTK_FOUND} )
    include( ${VTK_USE_FILE} )
    foreach( name ${FLOW_EXEC_NEEDS_VTK} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name} ${VTK_LIBRARIES} ${Teem_LIBRARIES} spurt_math spurt_image ${Boost_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
    message( WARNING "VTK-dependent executables: " ${FTLE_EXEC_NEEDS_VTK}
             " cannot be built because VTK was not found" )
endif()
