
set(BROKEN_FLOW_EXEC
    best-scale
    color_ftle
    generate_seeds
    vector_ftle_euler
	flow_map_new
	generate_flow
    tangent_lines
    lavd_new
    lavd-vel_only
    lavd_intp
    lavd_multi
    vector_field
    vector_ftle
    tensor_ftle
)


set(FLOW_EXEC_SRC
    catseye
    doublegyre
    abc
    flowmap_in_3D_nrrd
	flow_map
    lavd
    trace_particles
    view_flowmap
    view_flow
    test_upsample
    parse_lavd_header
    #debug
    view_particles_path
    flow_map_unstructured
    flow_map_vtk
    time_dependent_flow_map
    time_dependent_flow_map_multi_step
    time_dependent_2d_flow_map_multi_step
    ftle
)

if( ${VTK_FOUND} )
    # include( ${VTK_USE_FILE} )
    foreach( name ${FLOW_EXEC_SRC} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
		target_link_libraries( ${name} ${VTK_LIBRARIES} ${Teem_LIBRARIES} BARG
			spurt_image spurt_flow ${Boost_LIBRARIES} 
            TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
            ${NetCDF_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
    message( WARNING "VTK-dependent executables: " ${FTLE_EXEC_WITH_VTK}
             " cannot be built because VTK was not found" )
endif()
