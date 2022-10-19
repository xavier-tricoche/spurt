# - Try to find STANN
# Once done this will define
#
#
#   STANN_FOUND - system has STANN
#   STANN_INCLUDE_DIR - STANN include directory

include(Utils)

if( STANN_USE_MINE )
    set( STANN_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/stann )
    # set( STANN_FIND_QUIETLY TRUE )
    set( STANN_FOUND TRUE )
endif()

if( STANN_INCLUDE_DIR )
    set( STANN_FIND_QUIETLY TRUE )
endif()

list( APPEND STANN_TENTATIVE_INC_DIR
    $ENV{HOME}/code/STANN/include
)

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND STANN_TENTATIVE_INC_DIR ${dir}/STANN )
endforeach()

find_path( STANN_INCLUDE_DIR 
    NAMES sfcnn.hpp
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${STANN_TENTATIVE_INC_DIR}
)

if( STANN_INCLUDE_DIR )
   set( STANN_FOUND TRUE )
else()
   set( STANN_FOUND FALSE )
endif()

if( STANN_FOUND )
  if( NOT STANN_FIND_QUIETLY )
    message( STATUS "Found STANN : ${STANN_INCLUDE_DIR}" )
  endif()
else()
  if( STANN_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find STANN" )
  endif()
endif()

mark_as_advanced( STANN_INCLUDE_DIR )