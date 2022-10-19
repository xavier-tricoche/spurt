# - Try to find Teem
# Once done this will define
#
#
#   Teem_FOUND       - system has Teem
#   Teem_INCLUDE_DIR - Teem include directories
#   Teem_LIBRARIES   - libraries to link against to use Teem

include(Utils)

set( LIB_NAME teem )

if( NOT Teem_FOUND )
    find_package( Teem QUIET NO_MODULE PATHS ${USER_LIKELY_LIBRARY_DIR} )
endif()
if( Teem_FOUND )
    include( ${Teem_USE_FILE} )
    set( Teem_INCLUDE_DIR ${Teem_INCLUDE_DIRS} )
    # we are doing a search for the library although we know its name and
    # location just to recover its complete path with extension
    set( Teem_LIBRARIES )
    find_library( Teem_LIBRARIES NAME ${LIB_NAME} HINTS ${Teem_LIBRARY_DIRS} )
endif()

if( Teem_FOUND )
    set( Teem_FIND_QUIETLY TRUE )
endif()

if( NOT Teem_FOUND )
    set( Teem_TENTATIVE_INC_DIR $ENV{HOME}/code/teem/include )
    foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
        list( APPEND Teem_TENTATIVE_INC_DIR ${dir}/teem )
    endforeach()

    find_path( Teem_INCLUDE_DIR
        NAMES teem/nrrd.h 
        HINTS ${USER_LIKELY_INCLUDE_DIR} ${Teem_TENTATIVE_INC_DIR}
    )
    find_library( Teem_LIBRARIES
        NAME ${LIB_NAME}
        HINTS ${USER_LIKELY_LIBRARY_DIR} ${Teem_TENTATIVE_LIB_DIR}
    )

    if( Teem_INCLUDE_DIRS AND Teem_LIBRARIES )
        set( Teem_FOUND TRUE )
    else()
        set( Teem_FOUND FALSE )
    endif()
endif()

if( Teem_FOUND )
    if( NOT Teem_FIND_QUIETLY )
        message( STATUS "Found Teem :"
            "\nINC: " ${Teem_INCLUDE_DIRS} 
            "\nLIB: " ${Teem_LIBRARIES} 
        )
    endif()
else()
    if( Teem_FIND_REQUIRED )
        message( FATAL_ERROR "Could NOT find Teem" )
    endif()
endif()

mark_as_advanced( Teem_INCLUDE_DIR Teem_LIBRARIES )

