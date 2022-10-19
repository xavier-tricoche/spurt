# - Try to find GSL - The GNU Scientific Library
# Once done this will define
#
#
#   GSL_FOUND - system has GSL
#   GSL_INCLUDE_DIR - GSL include directory
#   GSL_LIBRARIES - Link these to use GSL

include(Utils)

set( LIB_NAME gsl )

if( GSL_INCLUDE_DIR AND GSL_LIBRARIES )
    set( GSL_FIND_QUIETLY TRUE )
endif()

find_path( GSL_INCLUDE_DIR 
    NAMES gsl/gsl_version.h 
    HINTS ${USER_LIKELY_INCLUDE_DIR}
)
find_library( GSL_LIBRARIES
    NAME ${LIB_NAME}
    HINTS ${USER_LIKELY_LIBRARY_DIR}
)

if( GSL_INCLUDE_DIR AND GSL_LIBRARIES )
    set( GSL_FOUND TRUE )
else()
    set( GSL_FOUND FALSE )
endif()

if( GSL_FOUND )
    if( NOT GSL_FIND_QUIETLY )
        message( STATUS "Found GSL :"
            "\nINC: " ${GSL_INCLUDE_DIR} 
            "\nLIB: " ${GSL_LIBRARIES} 
        )
    endif()
else()
    if( GSL_FIND_REQUIRED )
        message( FATAL_ERROR "Could NOT find GSL" )
    endif()
endif()

mark_as_advanced( GSL_INCLUDE_DIR GSL_LIBRARIES )

