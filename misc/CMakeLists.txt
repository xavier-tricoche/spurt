set( MISC_SRC
    color_lines
    point2tris
    lmpd2vtk
    lmpd2vcl
)
foreach( src ${MISC_SRC} )
    add_executable( ${src} ${src}.cpp )
    target_link_libraries( ${src} ${Teem_LIBRARIES} )
    install( TARGETS ${src} RUNTIME DESTINATION bin )
endforeach()

add_executable( test_option_parse test_option_parse.cpp )
target_link_libraries( test_option_parse ${Boost_LIBRARIES} )
install( TARGETS test_option_parse RUNTIME DESTINATION bin )
