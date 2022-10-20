set( CREASE_SRC
    measure_wrapper.cpp
    extractor.cpp
    pvo.cpp
)
add_library( spurt_crease ${LIBMODE} ${CREASE_SRC} )
target_link_libraries( spurt_crease ${Teem_LIBRARIES} spurt_image spurt_math )
install( TARGETS spurt_crease
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

if ( VTK_FOUND )
    add_executable( peikert_sadlo ${EXECUTABLE_FLAG} peikert_sadlo.cpp )
    target_link_libraries( peikert_sadlo
                           ${Teem_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${CELLTREE_LIBRARIES}
                           ${NetCDF_LIBRARIES}
                           ${VTK_LIBRARIES}
                           BARG)
    install( TARGETS ${name} RUNTIME DESTINATION bin )
endif()
