# - Try to find NetCDF
# Once done this will define
#
#
#   NetCDF_FOUND - system has NetCDF
#   NetCDF_INCLUDE_DIR - NetCDF include directory
#   NetCDF_LIBRARIES - Link these to use NetCDF

include(Utils)

set( LIB_NAME netcdf )

if( NetCDF_INCLUDE_DIR AND NetCDF_LIBRARIES )
    set( NetCDF_FIND_QUIETLY TRUE )
endif( NetCDF_INCLUDE_DIR AND NetCDF_LIBRARIES )

find_path( NetCDF_INCLUDE_DIR 
    NAMES netcdf.h 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${NetCDF_TENTATIVE_INC_DIR}
)
find_library( NetCDF_LIBRARIES
    NAME ${LIB_NAME}
    HINTS ${USER_LIKELY_LIBRARY_DIR} ${NetCDF_TENTATIVE_LIB_DIR}
)

if( NetCDF_INCLUDE_DIR AND NetCDF_LIBRARIES )
   set( NetCDF_FOUND TRUE )
else()
   set( NetCDF_FOUND FALSE )
endif()

if( NetCDF_FOUND )
    if( NOT NetCDF_FIND_QUIETLY )
        message( STATUS "Found NetCDF :"
            "\nINC: " ${NetCDF_INCLUDE_DIR} 
            "\nLIB: " ${NetCDF_LIBRARIES} 
        )
    endif()
else()
  if( NetCDF_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find NetCDF" )
  endif()
endif()

mark_as_advanced( NetCDF_INCLUDE_DIR NetCDF_LIBRARIES )

