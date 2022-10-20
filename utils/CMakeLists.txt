set(SRC
    scalar2d
    scalar3d
)


foreach( target ${SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} ${MatIO_LIBRARIES} spurt_image ${Boost_LIBRARIES})
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()
