set(ORBITAL_EXEC_SRC
    map_ftle
    map_ftle_new
    draw_map
)

foreach( target ${ORBITAL_EXEC_SRC} )
    add_executable( orbital_${target} ${target}.cpp )
    target_link_libraries( orbital_${target}
                           ${Teem_LIBRARIES}
                           ${PNG_LIBRARIES}
                           ${TBB_LIBRARIES}
                           ${TBB_LIBRARY_RELEASE} )
    install( TARGETS orbital_${target} RUNTIME DESTINATION bin )
endforeach()

# include( ${VTK_USE_FILE} )
add_executable( orbital_draw_map_interactive draw_map_interactive.cpp )
target_link_libraries( orbital_draw_map_interactive
                       ${VTK_LIBRARIES} ${Teem_LIBRARIES} )
install( TARGETS orbital_draw_map_interactive RUNTIME DESTINATION bin )
