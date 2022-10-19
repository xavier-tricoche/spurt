set( MATH_LIB_SRC
     matrix.cpp
)
add_library( spurt_math ${LIBMODE} ${MATH_LIB_SRC} )
target_link_libraries( spurt_math ${Teem_LIBRARIES} )
install( TARGETS spurt_math
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

set( MATH_EXEC
    test_rbfbasis
    voxelsubdiv
)

foreach( name ${MATH_EXEC} )
    add_executable( ${name} ${name}.cpp )
    target_link_libraries( ${name} spurt_math spurt_image ${Teem_LIBRARIES} )
endforeach()
install( TARGETS ${MATH_EXEC} RUNTIME DESTINATION bin )