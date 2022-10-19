set(POINCARE_SRC
    ls.cpp
    newton.cpp
    map.cpp
)
add_library( xavier_poincare ${LIBMODE} ${POINCARE_SRC} )
target_link_libraries( xavier_poincare ${Boost_LIBRARIES} )
install( TARGETS xavier_poincare
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)