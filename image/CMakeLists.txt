set(IMAGE_SRC
    probe.cpp
    integrate.cpp
    creaseline2d.cpp
)
add_library( spurt_image ${LIBMODE} ${IMAGE_SRC} )
target_link_libraries( spurt_image ${Teem_LIBRARIES} )
install( TARGETS spurt_image
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_executable( test_nrrd_manip test_nrrd_manip.cpp )
target_link_libraries( test_nrrd_manip ${Teem_LIBRARIES} )
install( TARGETS test_nrrd_manip RUNTIME DESTINATION bin )

add_executable( unux unux.cpp )
target_link_libraries( unux ${Teem_LIBRARIES} ${Boost_LIBRARIES} )
install( TARGETS unux RUNTIME DESTINATION bin )

add_executable( nrrd_stats nrrd_stats.cpp )
target_link_libraries( nrrd_stats ${Teem_LIBRARIES} ${Boost_LIBRARIES} )
install( TARGETS nrrd_stats RUNTIME DESTINATION bin )
