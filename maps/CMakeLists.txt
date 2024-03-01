set( MAPS_SRC
    # tokatopo.cpp
    # basic_definitions.cpp
)

set( MAPS_EXEC_SRC
    # q-based-fixed-points
    # q-based-fixedpoint-search
    # q-profile-fast
    # q-profile-period-specific
    # q-profile
    # q-version-1.1
    # sample_map
    # adaptive_map_sampling
)

if( ${CGARTH_TOKAMAK_FOUND} )
    add_library( xavier_maps ${LIBMODE} ${MAPS_SRC} )
    target_link_libraries(xavier_maps ${TEEM_LIBRARIES} ${CGARTH_TOKAMAK_LIBRARIES} xavier_poincare)
    install( TARGETS xavier_maps
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
    foreach(target ${MAPS_EXEC_SRC} )
        add_executable(${target} ${target}.cpp)
        target_link_libraries(${target} ${Teem_LIBRARIES} ${CGARTH_TOKAMAK_LIBRARIES} xavier_maps)
        install(TARGETS ${target} RUNTIME DESTINATION bin)
    endforeach()
else()
	message( WARNING "maps library and executables: " ${MAPS_EXEC_SRC} " cannot be built."
        "\nTokamak wrapper not available in this installation."
    )
endif()

add_executable(symplectic4d symplectic4d.cpp)
target_link_libraries(symplectic4d ${Teem_LIBRARIES} ${Boost_LIBRARIES} /usr/local/lib/libtbb.dylib /usr/local/lib/libtbbmalloc.dylib /usr/local/lib/libtbbmalloc_proxy.dylib)
install(TARGETS symplectic4d RUNTIME DESTINATION bin)


# add_executable(map4d map4d.cpp)
# target_link_libraries(map4d ${Teem_LIBRARIES} ${Boost_LIBRARIES} /usr/local/lib/libtbb.dylib /usr/local/lib/libtbbmalloc.dylib /usr/local/lib/libtbbmalloc_proxy.dylib)
# install(TARGETS map4d RUNTIME DESTINATION bin)
