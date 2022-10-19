# - Try to find NVIS
# Once done this will define
#
#
#   NVIS_FOUND - system has NVIS
#   NVIS_INCLUDE_DIR - NVIS include directory

include(Utils)

if( NVIS_USE_MINE )
    set( NVIS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/nvis )
    # set( NVIS_FIND_QUIETLY TRUE )
    set( NVIS_FOUND TRUE )
endif()

if( NVIS_INCLUDE_DIR )
    set( NVIS_FIND_QUIETLY TRUE )
endif()

list( APPEND NVIS_TENTATIVE_INC_DIR
    $ENV{HOME}/code/nvis
)

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND NVIS_TENTATIVE_INC_DIR ${dir}/nvis )
endforeach()

find_path( NVIS_INCLUDE_DIR 
    NAMES math/fixed_vector.hpp 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${NVIS_TENTATIVE_INC_DIR}
)

if( NVIS_INCLUDE_DIR )
   set( NVIS_FOUND TRUE )
else()
   set( NVIS_FOUND FALSE )
endif()

if( NVIS_FOUND )
  if( NOT NVIS_FIND_QUIETLY )
    message( STATUS "Found NVIS : ${NVIS_INCLUDE_DIR}" )
  endif()
else()
  if( NVIS_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find NVIS" )
  endif()
endif()

mark_as_advanced( NVIS_INCLUDE_DIR )