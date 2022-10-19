set( XAVIER_GRAPHICS_SRCS
	make_color_scale
)

foreach( name ${XAVIER_GRAPHICS_SRCS} )
    add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
    target_link_libraries( ${name} ${Boost_LIBRARIES} )
    install( TARGETS ${name} RUNTIME DESTINATION bin )
endforeach()
