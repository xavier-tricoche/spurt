# - Try to find GMM
# Once done this will define
#
#
#   GMM_FOUND - system has GMM
#   GMM_INCLUDE_DIR - GMM include directory

include(Utils)

if( GMM_INCLUDE_DIR )
    set( GMM_FIND_QUIETLY TRUE )
endif()

list( APPEND GMM_TENTATIVE_INC_DIR
    $ENV{HOME}/code/gmm/include
    $ENV{HOME}/code/getfem/include
    $ENV{HOME}/code/getfem-4.1/include
    $ENV{HOME}/code/getfem-4.1.1/include
)

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND GMM_TENTATIVE_INC_DIR ${dir}/gmm )
endforeach()

find_path( GMM_INCLUDE_DIR 
    NAMES gmm/gmm.h 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${GMM_TENTATIVE_INC_DIR}
)

if( GMM_INCLUDE_DIR )
   set( GMM_FOUND TRUE )
else()
   set( GMM_FOUND FALSE )
endif()

if( GMM_FOUND )
  if( NOT GMM_FIND_QUIETLY )
    message( STATUS "Found GMM : ${GMM_INCLUDE_DIR}" )
  endif()
else()
  if( GMM_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find GMM" )
  endif()
endif()

mark_as_advanced( GMM_INCLUDE_DIR )