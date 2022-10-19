# - Try to find MatIO
# Once done this will define
#
#   MatIO_FOUND - system has MatIO
#   MatIO_INCLUDE_DIR - MatIO include directory
#   MatIO_LIBRARIES - Link these to use MatIO

include(Utils)

if( MatIO_INCLUDE_DIR AND MatIO_LIBRARIES )
    set( MatIO_FIND_QUIETLY TRUE )
endif()

set( LIB_NAME matio )

list( APPEND MatIO_TENTATIVE_INC_DIR
    $ENV{HOME}/code/MatIO
)

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND MatIO_TENTATIVE_INC_DIR ${dir}/MatIO )
endforeach()

find_path( MatIO_INCLUDE_DIR 
    NAMES matio.h
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${MatIO_TENTATIVE_INC_DIR}
)
find_library( MatIO_LIBRARIES
    NAME ${LIB_NAME}
    HINTS ${USER_LIKELY_LIBRARY_DIR} ${GLUI_TENTATIVE_LIB_DIR}
)

if( MatIO_INCLUDE_DIR AND MatIO_LIBRARIES )
   set( MatIO_FOUND TRUE )
else()
   set( MatIO_FOUND FALSE )
endif()

if( MatIO_FOUND )
  if( NOT MatIO_FIND_QUIETLY )
    message( STATUS "Found MatIO : ${MatIO_INCLUDE_DIR}" )
  endif()
else()
  if( MatIO_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find MatIO" )
  endif()
endif()

mark_as_advanced( MatIO_INCLUDE_DIR MatIO_LIBRARIES )