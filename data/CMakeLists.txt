set( DATA_EXEC_SRC
	test_locator
    test_image
    test_heap
)

include_directories( ${KDTREEPP_INCLUDE_DIR} )
foreach( target ${DATA_EXEC_SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries(${target} /usr/local/lib/libtbb.dylib /usr/local/lib/libtbbmalloc.dylib /usr/local/lib/libtbbmalloc_proxy.dylib ${Boost_LIBRARIES})
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()
