# set( Spurt_CREASE_SRC
#     measure_wrapper.cpp
#     extractor.cpp
#     pvo.cpp
#     face.cpp
#     crease.cpp
# )
# add_library( spurt_crease ${LIBMODE} ${CREASE_SRC} )
# target_link_libraries( spurt_crease ${Teem_LIBRARIES} spurt_image )
# install( TARGETS spurt_crease
#   LIBRARY DESTINATION lib
#   ARCHIVE DESTINATION lib
# )

include_directories(
    .
    ${Spurt_SOURCE_DIR}
)

set( Spurt_CREASE_EXEC_SRC_WITH_VTK
    peikert_sadlo
    peikert
    multires_peikert
    schultz
    parallel_schultz
    schultz_Jan2025
    schultz_Feb2025
    filter_crease
)

set( Spurt_CREASE_EXEC_SRC_NO_VTK
    gH_sspace_selection
    gH_filtered_sspace_selection
)

if ( VTK_FOUND )
    foreach( exec ${Spurt_CREASE_EXEC_SRC_WITH_VTK} )
        add_executable( ${exec} ${EXECUTABLE_FLAG} ${exec}.cpp )
        target_link_libraries(${exec}
                           ${Teem_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${CELLTREE_LIBRARIES}
                           ${NetCDF_LIBRARIES}
                           ${VTK_LIBRARIES}
                           spurt_barg )
        install( TARGETS ${exec} RUNTIME DESTINATION bin )
    endforeach()
endif()

foreach( exec ${Spurt_CREASE_EXEC_SRC_NO_VTK} )
    add_executable( ${exec} ${EXECUTABLE_FLAG} ${exec}.cpp ) 
    target_link_libraries( ${exec}
                       ${Teem_LIBRARIES}
                       ${Boost_LIBRARIES}
                       TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy  )
    install( TARGETS gH_sspace_selection RUNTIME DESTINATION bin )
endforeach()