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
# /Users/xmt/code/boost-trunk/deb.shared/lib/libboost_program_options.dylib
# /Users/xmt/code/boost-trunk/deb.shared/lib/libboost_filesystem.dylib
# /Users/xmt/code/boost-trunk/deb.shared/lib/libboost_system.dylib )
install( TARGETS test_option_parse RUNTIME DESTINATION bin )
