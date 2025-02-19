file(GLOB SOURCES
    *.h
    *.cpp
	*.cxx
	)


if( ${VTK_FOUND} )
    #make this a library
	add_library( spurt_barg ${LIBMODE} ${SOURCES} )
	target_link_libraries( spurt_barg ${VTK_LIBRARIES} )
else()
    message( WARNING 
             "VTK-dependent executables for the boundary-aware-rectgrid project"
             " cannot be built because VTK was not found" )
endif()


