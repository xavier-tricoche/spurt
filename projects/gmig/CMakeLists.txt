project( MAARTEN_PROJECT )

set( MAARTEN_CORE_SRCS
    rbf_1station
    new_anisotropy
    # merge_files
    smooth_reconstruction
    single_source_field
)

set( MAARTEN_NEEDS_VTK_AND_KDTREEPP_SRCS
    travel_time
    view_data
    view_anisotropy
)

set( MAARTEN_NEEDS_GSL 
    anisotropy
    fitex
)

foreach( name ${MAARTEN_CORE_SRCS} )
    add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
    target_link_libraries( ${name} 
                           ${Teem_LIBRARIES}
                           ${Boost_LIBRARIES} )
    install( TARGETS ${name} RUNTIME DESTINATION bin )
endforeach()

if( VTK_FOUND AND KDTREEPP_FOUND )
    include( ${VTK_USE_FILE} )
    foreach( name ${MAARTEN_NEEDS_VTK_AND_KDTREEPP_SRCS} )
        add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
          target_link_libraries( ${name}
                                 ${Teem_LIBRARIES}
                                 ${VTK_LIBRARIES}
                                 ${Boost_LIBRARIES} )
          install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
endif()

if( GSL_FOUND )
    include_directories( ${GSL_INCLUDE_DIRS} )
    foreach( name ${MAARTEN_NEEDS_GSL} )
        add_executable( ${name} ${name}.cpp )
        target_link_libraries( ${name} 
                               ${GSL_LIBRARIES}
                               ${Teem_LIBRARIES}
                               ${Boost_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
endif()

