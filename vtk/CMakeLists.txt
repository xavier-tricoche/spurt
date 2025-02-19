set( SPURT_VTK_SRCS
	vectors_from_scalars
    streamlines
    view_colormap
	color_scalar_image
	volume_render
    volume_renderer
    test_lic
    view_lines
    test_jacobian
    test_filter_ccs
)

if( ${VTK_FOUND} )
    # include( ${VTK_USE_FILE} )
    foreach( name ${SPURT_VTK_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name} ${VTK_LIBRARIES} ${Teem_LIBRARIES} ${Boost_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "VTK executables: " ${SPURT_VTK_SRCS} " cannot be built."
        "\nVTK is not available in this installation."
    )
endif()
