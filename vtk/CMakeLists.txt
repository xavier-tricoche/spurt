set( XAVIER_VTK_SRCS
	vectors_from_scalars
    view_grains
    view_grains_animation
    view_stress
    view_stress_static
    view_grains_static
    streamlines
    # grains_orientation_by_size
    new_grains_orientation_by_size
    # orientation_segmentation
    view_colormap
	color_scalar_image
    # new_grain_stress
    # new_streamlines
    values_on_grains
	volume_render
    volume_renderer
    streamlines_density
    streamlines_grains
    # test_lic
    view_lines
    test_jacobian
)

if( ${VTK_FOUND} )
    # include( ${VTK_USE_FILE} )
    foreach( name ${XAVIER_VTK_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name} ${VTK_LIBRARIES} ${Teem_LIBRARIES} ${Boost_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "VTK executables: " ${XAVIER_VTK_SRCS} " cannot be built."
        "\nVTK is not available in this installation."
    )
endif()
