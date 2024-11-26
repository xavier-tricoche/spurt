set( MATH_TEST
    test_exponentiation
    test_polynomial
    test_types
)
foreach( name ${MATH_TEST} )
    add_executable( ${name} ${name}.cpp )
endforeach()
install( TARGETS ${MATH_TEST} RUNTIME DESTINATION bin )

set( MATH_EXEC
    test_rbfbasis
    test_svd
    voxelsubdiv
)

foreach( name ${MATH_EXEC} )
    add_executable( ${name} ${name}.cpp )
    target_link_libraries( ${name} spurt_image ${Teem_LIBRARIES} )
endforeach()
install( TARGETS ${MATH_EXEC} RUNTIME DESTINATION bin )

set( MATH_EXEC_WITH_TOKAMAK
     # divergence_cleaning
     # divergence_cleaning_tokamak
     # divergence_cleaning_tokamak_neumann
     # naive_divergence_cleaning
     # test_poisson
)

foreach( name ${MATH_EXEC_WITH_TOKAMAK} )
    add_executable( ${name} ${name}.cpp )
    target_link_libraries( ${name} spurt_image ${Teem_LIBRARIES} ${TOKAMAK_LIBRARIES})
endforeach()
install( TARGETS ${MATH_EXEC_WITH_TOKAMAK} RUNTIME DESTINATION bin )
