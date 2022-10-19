# This is the root spurt CMakeLists file
cmake_minimum_required(VERSION 2.8)

# Run cmake --help-policy CMP<num> to see documentation.
if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
  cmake_policy( SET CMP0012 NEW )
endif( COMMAND cmake_policy )

project( spurt )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" )
set( THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/third_party )
set( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${THIRD_PARTY_DIR} )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release Public"
      FORCE)
endif()

option( BUILD_OPTION_OPENMP "use openmp parallelization" OFF )

set( LIBMODE SHARED )

#######################################################################

set( SPURT_SOURCE_DIR ${CMAKE_SOURCE_DIR} )

set( NVIS_INCLUDE_DIR     ${CMAKE_SOURCE_DIR}/third_party/nvis )
set( KDTREEPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/kdtree++ )
set( STANN_INCLUDE_DIR    ${CMAKE_SOURCE_DIR}/third_party/stann )
set( CELLTREE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/celltree/src )

# required packages
set( Boost_USE_STATIC_LIBS      OFF )
set( Boost_USE_MULTITHREADED    ON  )
set( Boost_USE_STATIC_RUNTIME   OFF )
find_package( Boost 1.55 REQUIRED
              COMPONENTS
              filesystem
              system
              program_options )  # C++ goodies
find_package( Eigen  REQUIRED )  # linear algebra
find_package( LAPACK REQUIRED )  # linear algebra
find_package( Teem   REQUIRED )  # I/O / image processing
find_package( PNG    REQUIRED )  # I/O


# optional packages

find_package( GLUI )      # graphics
find_package( GLUT )      # graphics
find_package( GSL )       # numerical analysis
find_package( HDF5 )      # I/O
find_package( MPI )       # I/O
find_package( OpenGL )    # graphics
find_package( NetCDF )    # I/O
find_package( VTK )       # graphics / visualization
find_package( ZLIB )      # I/O

# build configuration

# set( Tapkee_INCLUDE_DIR ${THIRD_PARTY_DIR}/tapkee )
# set( ExprTK_INCLUDE_DIR ${THIRD_PARTY_DIR}/fastmathparser )
set( KDTREEPP_INCLUDE_DIR ${THIRD_PARTY_DIR}/kdtree++ )
set( STANN_INCLUDE_DIR ${THIRD_PARTY_DIR}/stann )


# mandatory includes
include_directories( ${Teem_INCLUDE_DIRS} )
include_directories( ${XAVIER_SOURCE_DIR} )
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${Eigen_INCLUDE_DIR} )
include_directories( ${PNG_INCLUDE_DIRS} )
# include_directories( ${Tapkee_INCLUDE_DIR} )
# include_directories( ${ExprTK_INCLUDE_DIR} )
include_directories( ${KDTREEPP_INCLUDE_DIR} )
include_directories( ${STANN_INCLUDE_DIR} )
include_directories( ${THIRD_PARTY_DIR} )
add_definitions( ${PNG_DEFINITIONS} )
include_directories( . )

set( AVAILABLE_PACKAGES "" )
set( MISSING_PACKAGES "" )

if( VTK_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "VTK" )

    find_path( VTK_DAG
        NAMES vtkDirectedAcyclicGraph.h
        HINTS ${VTK_INCLUDE_DIRS}
    )
    if ( VTK_DAG )
        # message( WARNING "vtkDirectedAcyclicGraph.h found" )
        set( VTK_DAG_FOUND TRUE )
        add_definitions( -D__VTK_HAS_DAG__ )
        list( APPEND ${AVAILABLE_PACKAGES} "VTK_DAG" )
    else()
        message( WARNING "vtkDirectedAcyclicGraph.h NOT FOUND" )
        set( VTK_DAG_FOUND FALSE )
        list( APPEND ${MISSING_PACKAGES} "VTK_DAG" )
    endif()

    find_path( VTK_AXIS
        NAMES vtkAxisActor.h
        HINTS ${VTK_INCLUDE_DIRS}
    )
    if( VTK_AXIS )
        set( VTK_AXISACTOR_FOUND TRUE )
        add_definitions( -D__VTK_HAS_AXIS_ACTOR__ )
        list( APPEND ${AVAILABLE_PACKAGES} "VTK_AXIS" )
    else()
        message( WARNING "vtkAxisActor.h NOT FOUND" )
        set( VTK_AXISACTOR_FOUND FALSE )
        list( APPEND ${MISSING_PACKAGES} "VTK_AXIS" )
    endif()
else()
    list( APPEND ${MISSING_PACKAGES} "VTK" )
endif()

# GL* bundle
if( OPENGL_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "OPENGL" )
else()
    list( APPEND ${MISSING_PACKAGES} "OPENGL" )
endif()

if( GLUT_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "GLUT" )
else()
    list( APPEND ${MISSING_PACKAGES} "GLUT" )
endif()

if( GLUI_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "GLUI" )
else()
    list( APPEND ${MISSING_PACKAGES} "GLUI" )
endif()

if( OPENGL_FOUND AND GLUT_FOUND AND GLUI_FOUND )
    set( GL_BUNDLE_FOUND TRUE )
    list( APPEND GL_BUNDLE_INCLUDE_DIR
            ${OPENGL_INCLUDE_DIR}
            ${GLUT_INCLUDE_DIR}
            ${GLUI_INCLUDE_DIR}
        )
    list( APPEND GL_BUNDLE_LIBRARIES
            ${OPENGL_LIBRARIES}
            ${GLUT_LIBRARIES}
            ${GLUI_LIBRARIES}
        )
else()
    set( GL_BUNDLE_FOUND FALSE )
endif()

if( NetCDF_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "NetCDF" )
else()
    list( APPEND ${MISSING_PACKAGES} "NetCDF" )
endif()

# if both VTK and NetCDF are available we need to link against
# VTK_NETCDF to deal with cases where VTK's <netcdf.h> is included
# instead of the system's one (sigh...)
# if( VTK_FOUND AND NetCDF_FOUND )
#     find_library( VTK_NETCDF_LIBRARY
#         NAME "vtkNetCDF"
#         HINTS ${VTK_LIBRARY_DIRS}
#     )
# endif()

if( HDF5_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "HDF5" )
else()
    list( APPEND ${MISSING_PACKAGES} "HDF5" )
endif()

if( MPI_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "MPI" )
else()
    list( APPEND ${MISSING_PACKAGES} "MPI" )
endif()

if( ZLIB_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "ZLIB" )
else()
    list( APPEND ${MISSING_PACKAGES} "ZLIB" )
endif()

if( NOT HDF5_FOUND )
    if ( HDF5_hdf5_LIBRARY_RELEASE )
        list( APPEND HDF5_LIBRARIES ${HDF5_hdf5_LIBRARY_RELEASE} )
        set( HDF5_FOUND TRUE )
    endif()
endif()

# HDF5 bundle
if( HDF5_FOUND AND MPI_FOUND AND ZLIB_FOUND )
    set( HDF_BUNDLE_FOUND TRUE )
    list( APPEND HDF_BUNDLE_INCLUDE_DIR
            ${HDF5_INCLUDE_DIR}
            ${MPI_INCLUDE_DIR}
            ${ZLIB_INCLUDE_DIR}
        )
    list( APPEND HDF_BUNDLE_LIBRARIES
            ${HDF5_LIBRARIES}
            ${MPI_LIBRARIES}
            ${ZLIB_LIBRARIES}
        )
else()
    set( HDF_BUNDLE_FOUND FALSE )
endif()

if( GSL_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "GSL" )
else()
    list( APPEND ${MISSING_PACKAGES} "GSL" )
endif()

if( LAPACK_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "LAPACK" )
else()
    list( APPEND ${MISSING_PACKAGES} "LAPACK" )
endif()

# GSL bundle
if ( GSL_FOUND AND LAPACK_FOUND )
    list( APPEND GSL_LIBRARIES ${LAPACK_LIBRARIES} )
else()
    set( GSL_FOUND FALSE )
endif()

#######################################################################

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-result -Wno-deprecated-declarations" )

if( BUILD_OPTION_OPENMP )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
endif()

#######################################################################
# set build product destinations

set( LIBRARY_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/lib )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin )

link_directories( ${LIBRARY_OUTPUT_PATH} )

# individual libs
list( APPEND SPURT_LIBRARIES
    crease
    data
    format
    flow
    graphics
    image
    math
    misc
    reconstruction
)

foreach( lib ${SPURT_LIBRARIES})
  add_subdirectory( ${lib} )
endforeach()

######################################################################
# determine targets to be build and define their linking
######################################################################

IF( NOT SPURT_INSTALL_PRESET_IS_SET )
  SET( CMAKE_INSTALL_PREFIX "$ENV{HOME}/SPURT.default_install_dir" CACHE PATH "Install path prefix, prepended onto install directories." FORCE )
  SET( SPURT_INSTALL_PRESET_IS_SET 1 CACHE INTERNAL "Default install path written" FORCE )
ENDIF()

install( DIRECTORY bin DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY scripts DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY utils/GLSLShaderCode DESTINATION utils USE_SOURCE_PERMISSIONS)


######################################################################
# install all binaries
######################################################################
foreach( target ${SPURT_TARGETS} )
  install( TARGETS ${target}
        RUNTIME DESTINATION bin
  )
endforeach( target ${SPURT_TARGETS} )
