set(SRC
    scalar2d
    scalar3d
    abc_function
    sinc
    marschner_lobb
)

foreach( target ${SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} ${MatIO_LIBRARIES} spurt_image ${Boost_LIBRARIES} TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()
