set( CREASE_SRC
    measure_wrapper.cpp
    extractor.cpp
    pvo.cpp
)
add_library( xavier_crease ${LIBMODE} ${CREASE_SRC} )
target_link_libraries( xavier_crease ${Teem_LIBRARIES} xavier_image xavier_math )
install( TARGETS xavier_crease
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
