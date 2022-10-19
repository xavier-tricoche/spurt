project( XAVIER_LEARNING )

set( LEARNING_SRCS
    test_isomap
)

include_directories( ${KDTREEPP_INCLUDE_DIR} )

if ( VTK_FOUND )
    # include( ${VTK_USE_FILE} )
    foreach( name ${LEARNING_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name}
                               ${Boost_LIBRARIES}
                               ${VTK_LIBRARIES}
                               ${Teem_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
endif()

set( LEARNING_NEED_TAPKEE_SRCS
  #  test_tapkee
)

if ( VTK_FOUND )
    # include( ${VTK_USE_FILE} )
    include_directories( ${Tapkee_INCLUDE_DIR} )
    foreach( name ${LEARNING_NEED_TAPKEE_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
        target_link_libraries( ${name}
                               ${Boost_LIBRARIES}
                               ${VTK_LIBRARIES}
                               ${Teem_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
endif()
