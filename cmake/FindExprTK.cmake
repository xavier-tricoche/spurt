# - Try to find exprtk header
# Once done this will define
#
#   ExprTK_FOUND - system has exprtk.hpp
#   ExprTK_INCLUDE_DIR - exprtk.hpp include directory

include(Utils)

if( ExprTK_USE_MINE )
    set( ExprTK_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/fastmathparser )
    # set( ExprTK_FIND_QUIETLY TRUE )
    set( ExprTK_FOUND TRUE )
endif()

if( ExprTK_INCLUDE_DIR )
    set( ExprTK_FIND_QUIETLY TRUE )
    set( ExprTK_FOUND TRUE )
endif()

list( APPEND ExprTK_TENTATIVE_INC_DIR
    $ENV{HOME}/code/fastmathparser
)

find_path( ExprTK_INCLUDE_DIR 
    NAMES exprtk.hpp
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${ExprTK_TENTATIVE_INC_DIR}
)

if( ExprTK_INCLUDE_DIR )
   set( ExprTK_FOUND TRUE )
else()
   set( ExprTK_FOUND FALSE )
   if( ExprTK_FIND_REQUIRED )
       message( FATAL_ERROR "Could NOT find ExprTK" )
   endif()
endif()

mark_as_advanced( ExprTK_INCLUDE_DIR )
