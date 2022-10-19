set( TOPOLOGY_EXEC_SRC
    scalar_topology
)

foreach( target ${TOPOLOGY_EXEC_SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} ${VTK_LIBRARIES} xavier_math xavier_image )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()