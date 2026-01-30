if( ${VTK_FOUND} )
	include( ${VTK_USE_FILE} )
    set( SPURT_WAVES_SRCS
        view_front
    )
    foreach( name ${SPURT_WAVES_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name} ${Teem_LIBRARIES} ${VTK_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "wave executables: " ${SPURT_WAVE_SRCS} " cannot be built."
	    "\nVTK is not available in this installation."
	)
endif()
