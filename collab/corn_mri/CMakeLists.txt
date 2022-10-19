project( CORN_MRI_PROJECT )

# find_package( Qt5 )

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_BINARY_DIR}
)

# set( CORN_MRI_LIB_SRCS
#     analysis
# )

# if ( VTK_FOUND AND QT_FOUND )
#     include( ${QT_USE_FILE} )
#
#     # Instruct CMake to run moc automatically when needed.
#     set(CMAKE_AUTOMOC ON)
#
#     set( QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "" )
#     set( QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "" )
#     set( QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "" )
#
#     set( SRC analysis.cpp main.cpp )
#     set( HDR interface.hpp )
#     set( MOC_HDR interface.hpp )
#     set( UI main_window.ui )
#     qt5_wrap_ui( UI_SRC ${UI} )
#     add_executable( mri_analysis ${EXECUTABLE_FLAG} ${SRC} ${UI_SRC} ${MOC_HDR} )
#     target_link_libraries( mri_analysis ${Teem_LIBRARIES} ${Boost_LIBRARIES}
#                        ${VTK_LIBRARIES} ${QT_LIBRARIES} )
#     install( TARGETS mri_analysis RUNTIME DESTINATION bin
#              BUNDLE DESTINATION bin )
# else()
#     message( WARNING "mri_analysis will not be built")
# endif()

set( CORN_MRI_SRCS
    measure
)

if( ${VTK_FOUND} )
	include( ${VTK_USE_FILE} )
    add_executable( measure ${EXECUTABLE_FLAGS} measure.cpp )
    target_link_libraries( measure ${Teem_LIBRARIES} ${Boost_LIBRARIES}
                          ${VTK_LIBRARIES} )
endif()
