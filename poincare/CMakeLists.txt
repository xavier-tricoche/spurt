set(POINCARE_SRC
    ls.cpp
    newton.cpp
    map.cpp
)
add_library( spurt_poincare ${LIBMODE} ${POINCARE_SRC} )
target_link_libraries( spurt_poincare ${Boost_LIBRARIES} )
install( TARGETS spurt_poincare
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
