# - Try to find Eigen
# Once done this will define
#
#   Eigen_FOUND - system has Eigen
#   Eigen_INCLUDE_DIR - Eigen include directory

include(Utils)

if( Eigen_INCLUDE_DIR )
    set( Eigen_FIND_QUIETLY TRUE )
endif()

list( APPEND Eigen_TENTATIVE_INC_DIR
    $ENV{HOME}/code/eigen
    $ENV{HOME}/code/eigen3
    $ENV{HOME}/code/eigen3.1.2
)

foreach( dir ${USER_LIKELY_INCLUDE_DIR} )
    list( APPEND Eigen_TENTATIVE_INC_DIR ${dir}/eigen ${dir}/eigen3 )
endforeach()

find_path( Eigen_INCLUDE_DIR 
    NAMES Eigen/Eigen 
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${Eigen_TENTATIVE_INC_DIR}
)

if( Eigen_INCLUDE_DIR )
   set( Eigen_FOUND TRUE )
else()
   set( Eigen_FOUND FALSE )
endif()

if( Eigen_FOUND )
  if( NOT Eigen_FIND_QUIETLY )
    message( STATUS "Found Eigen : ${Eigen_INCLUDE_DIR}" )
  endif()
else()
  if( Eigen_FIND_REQUIRED )
    message( FATAL_ERROR "Could NOT find Eigen" )
  endif()
endif()

mark_as_advanced( Eigen_INCLUDE_DIR )