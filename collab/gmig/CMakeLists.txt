project( GMIG_PROJECT )

find_package( Qt4 )
find_package( Tapkee )

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_BINARY_DIR}
)

set( GMIG_LIB_SRCS
    data_IO
    anisotropy
)

set( GMIG_NEEDS_VTK_AND_KDTREEPP_SRCS
    travel_time_vis
#    view_data
#    view_anisotropy
)

set( GMIG_NEEDS_GSL 
)

if ( QT_FOUND AND Tapkee_FOUND)
    include( ${QT_USE_FILE} )
    include_directories( ${Tapkee_INCLUDE_DIR} )
    
    # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)

    set( QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "" )
    set( QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "" )
    set( QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "" )

    set( SRC mantle_vis.cpp mantle_vis_main.cpp )
    set( HDR mantle_vis_renderer.hpp mantle_vis_control.hpp )
    set( MOC_HDR mantle_vis_renderer.hpp mantle_vis_control.hpp )
    set( UI mantle_vis_renderer.ui mantle_vis_control.ui )
    qt4_wrap_ui( UI_SRC ${UI} )
    add_executable( mantle_vis ${EXECUTABLE_FLAG} ${SRC} ${UI_SRC} ${MOC_HDR} )
    target_link_libraries( mantle_vis ${Teem_LIBRARIES} ${Boost_LIBRARIES} 
                       ${VTK_LIBRARIES} ${QT_LIBRARIES} )
    install( TARGETS mantle_vis RUNTIME DESTINATION bin 
             BUNDLE DESTINATION bin )
else()
    message( WARNING "mantle_vis will not be built")
    if (QT_FOUND)
        message(WARNING "Qt was found")
    else()
        message(WARNING "Qt was not found")
    endif()
    if (Tapkee_FOUND)
        message(WARNING "Tapkee was found")
    else()
        message(WARNING "Tapkee was not found")
    endif()
endif()

if ( VTK_FOUND ) 
    # include( ${VTK_USE_FILE} )
    # add_executable( mantle_vis ${EXECUTABLE_FLAG} mantle_vis.cpp)
    # target_link_libraries( mantle_vis
    #                        ${Teem_LIBRARIES}
    #                        ${Boost_LIBRARIES}
    #                        ${VTK_LIBRARIES} )
    # install( TARGETS mantle_vis RUNTIME DESTINATION bin )


    # include( ${VTK_USE_FILE} )
    # add_executable( alt.mantle_vis ${EXECUTABLE_FLAG} alt.mantle_vis.cpp)
    # target_link_libraries( alt.mantle_vis
    #                        ${Teem_LIBRARIES}
    #                        ${Boost_LIBRARIES}
    #                        ${VTK_LIBRARIES} )
    # install( TARGETS alt.mantle_vis RUNTIME DESTINATION bin )
endif()

foreach( name ${GMIG_CORE_SRCS} )
    add_executable( ${name} ${EXECUTABLE_FLAG} ${name}.cpp )
    target_link_libraries( ${name} 
                           ${Teem_LIBRARIES}
                           ${Boost_LIBRARIES} )
    install( TARGETS ${name} RUNTIME DESTINATION bin )
endforeach()

if( VTK_FOUND AND KDTREEPP_FOUND AND VTK_DAG_FOUND )
    include( ${VTK_USE_FILE} )
    foreach( name ${GMIG_NEEDS_VTK_AND_KDTREEPP_SRCS} )
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
    foreach( name ${GMIG_NEEDS_GSL} )
        add_executable( ${name} ${name}.cpp )
        target_link_libraries( ${name} 
                               ${GSL_LIBRARIES}
                               ${Teem_LIBRARIES}
                               ${Boost_LIBRARIES} )
        install( TARGETS ${name} RUNTIME DESTINATION bin )
    endforeach()
endif()

