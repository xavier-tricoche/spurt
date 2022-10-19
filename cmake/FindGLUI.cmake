# - Try to find GLUI
# Once done this will define
#
#
#   GLUI_FOUND - system has GLUI
#   GLUI_INCLUDE_DIR - GLUI include directory
#   GLUI_LIBRARIES - Link these to use GLUI

include(Utils)

if( GLUI_INCLUDE_DIR AND GLUI_LIBRARIES )
    set( GLUI_FIND_QUIETLY TRUE )
endif()

set( LIB_NAME glui )

list( APPEND GLUI_TENTATIVE_SRC_DIR
    $ENV{HOME}/code/glui
    $ENV{HOME}/code/glui/src
)

foreach( dir ${GLUI_TENTATIVE_SRC_DIR} )
    list( APPEND GLUI_TENTATIVE_INC_DIR ${dir} )
    list( APPEND GLUI_TENTATIVE_INC_DIR ${dir}/include )
    list( APPEND GLUI_TENTATIVE_INC_DIR ${dir}/GL )
    list( APPEND GLUI_TENTATIVE_INC_DIR ${dir}/include/GL )
    list( APPEND GLUI_TENTATIVE_LIB_DIR ${dir}/lib )
endforeach()

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND GLUI_TENTATIVE_INC_DIR ${dir}/glui )
    list( APPEND GLUI_TENTATIVE_INC_DIR ${dir}/GL)
endforeach()

find_path( GLUI_INCLUDE_DIR 
    NAMES glui.h 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${GLUI_TENTATIVE_INC_DIR}
)
find_library( GLUI_LIBRARIES
    NAME ${LIB_NAME}
    HINTS ${USER_LIKELY_LIBRARY_DIR} ${GLUI_TENTATIVE_LIB_DIR}
)

if( GLUI_INCLUDE_DIR AND GLUI_LIBRARIES )
   set( GLUI_FOUND TRUE )
   list( APPEND GLUI_INCLUDE_DIR ${GLUI_INCLUDE_DIR}/GL )
else()
   set( GLUI_FOUND FALSE )
endif()

if( GLUI_FOUND )
    if( NOT GLUI_FIND_QUIETLY )
        message( STATUS "Found GLUI :"
            "\nINC: " ${GLUI_INCLUDE_DIR} 
            "\nLIB: " ${GLUI_LIBRARIES} 
        )
    endif()
else()
    if( GLUI_FIND_REQUIRED )
        message( FATAL_ERROR "Could NOT find GLUI" )
    endif()
endif()

mark_as_advanced( GLUI_INCLUDE_DIR GLUI_LIBRARIES )

