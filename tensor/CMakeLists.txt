set(TENSOR_SRC
	dpl
	tensor_lines
)

foreach( target ${TENSOR_SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} spurt_image ${Boost_LIBRARIES} )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()
