# - Try to find libkdtree++
# Once done this will define
#
#
#   KDTREEPP_FOUND - system has libkdtree++
#   KDTREPP_INCLUDE_DIR - libkdtree++ include directory

include(Utils)

if( KDTREEPP_USE_MINE )
    set( KDTREEPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/libkdtree )
    # set( KDTREEPP_FIND_QUIETLY TRUE )
    set( KDTREEPP_FOUND TRUE )
endif()

if( KDTREEPP_INCLUDE_DIR )
    set( KDTREEPP_FIND_QUIETLY TRUE )
endif()

list( APPEND KDTREEPP_TENTATIVE_INC_DIR
    $ENV{HOME}/code
    $ENV{HOME}/code/libkdtree
    $ENV{HOME}/code/libkdtree++
)

find_path( KDTREEPP_INCLUDE_DIR 
    NAMES kdtree++/kdtree.hpp 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${KDTREEPP_TENTATIVE_INC_DIR}
)

if( KDTREEPP_INCLUDE_DIR )
   set( KDTREEPP_FOUND TRUE )
else()
   set( KDTREEPP_FOUND FALSE )
endif()

if( KDTREEPP_FOUND )
  if( NOT KDTREEPP_FIND_QUIETLY )
    message( STATUS "Found KDTREEPP : ${KDTREEPP_INCLUDE_DIR}" )
  endif()
else()
  if( KDTREEPP_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find KDTREEPP" )
  endif()
endif()

mark_as_advanced( KDTREEPP_INCLUDE_DIR )