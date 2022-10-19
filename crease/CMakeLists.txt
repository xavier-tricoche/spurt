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