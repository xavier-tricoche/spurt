set( DATA_EXEC_SRC
    test_locator
)

include_directories( ${KDTREEPP_INCLUDE_DIR} )
foreach( target ${DATA_EXEC_SRC} )
    add_executable( ${target} ${target}.cpp )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()
