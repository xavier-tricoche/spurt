# - Try to find Tapkee
# Once done this will define
#
#
#   Tapkee_FOUND - system has Tapkee
#   Tapkee_INCLUDE_DIR - Tapkee include directory

include(Utils)

if( Tapkee_USE_MINE )
    set( Tapkee_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/tapkee )
    # set( Tapkee_FIND_QUIETLY TRUE )
    set( Tapkee_FOUND TRUE )
else()
	message( WARNING "Tapkee_USE_MINE is set to NO" )
endif()

if( Tapkee_INCLUDE_DIR )
    set( Tapkee_FIND_QUIETLY TRUE )
endif()

list( APPEND Tapkee_TENTATIVE_INC_DIR
    $ENV{HOME}/code/tapkee/include
)

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND Tapkee_TENTATIVE_INC_DIR ${dir}/tapkee )
endforeach()

find_path( Tapkee_INCLUDE_DIR 
    NAMES tapkee/tapkee.hpp 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${Tapkee_TENTATIVE_INC_DIR}
)

if( Tapkee_INCLUDE_DIR )
   set( Tapkee_FOUND TRUE )
else()
   set( Tapkee_FOUND FALSE )
endif()

if( Tapkee_FOUND )
  if( NOT Tapkee_FIND_QUIETLY )
    message( STATUS "Found Tapkee : ${Tapkee_INCLUDE_DIR}" )
  endif()
else()
  if( Tapkee_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find Tapkee" )
  endif()
endif()

mark_as_advanced( Tapkee_INCLUDE_DIR )