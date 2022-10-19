# This is the root Xavier CMakeLists file
cmake_minimum_required(VERSION 3.0)

# Run cmake --help-policy CMP<num> to see documentation.
if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
  cmake_policy( SET CMP0012 NEW )
endif( COMMAND cmake_policy )

project( XAVIER )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" )
set( THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/third_party )
set( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${THIRD_PARTY_DIR} )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release Public"
      FORCE)
endif()

option( BUILD_OPTION_OPENMP "use openmp parallelization" OFF )

option( BUILD_OPTION_TBB "use tbb parallelization" ON )

option( OLD_COMPILER "deactivate \"recent\" C++ features" OFF )

set( LIBMODE SHARED )

# use, i.e. don't skip the full RPATH for the build tree
#set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, use the install RPATH already
#set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" "/Applications/Xcode.app/Contents/Developer/Library/Frameworks/")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MACOSX_RPATH 1)
endif()

#######################################################################

set( XAVIER_SOURCE_DIR ${CMAKE_SOURCE_DIR} )

set( NVIS_INCLUDE_DIR     ${CMAKE_SOURCE_DIR}/third_party/nvis )
set( ExprTK_INCLUDE_DIR   ${CMAKE_SOURCE_DIR}/third_party/fastmathparser )
set( KDTREEPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/kdtree++ )
set( Tapkee_INCLUDE_DIR   ${CMAKE_SOURCE_DIR}/third_party/tapkee )
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
find_package( MatIO )     # I/O
find_package( MPI )       # I/O
find_package( OpenGL )    # graphics
find_package( NetCDF )    # I/O
find_package( VTK )       # graphics / visualization
find_package( ZLIB )      # I/O
find_package( TBB )       # parallelism

# build configuration

set( NVIS_INCLUDE_DIR ${THIRD_PARTY_DIR}/nvis )
set( Tapkee_INCLUDE_DIR ${THIRD_PARTY_DIR}/tapkee )
set( ExprTK_INCLUDE_DIR ${THIRD_PARTY_DIR}/fastmathparser )
set( KDTREEPP_INCLUDE_DIR ${THIRD_PARTY_DIR}/kdtree++ )
set( TOKAMAK_INCLUDE_DIR ${THIRD_PARTY_DIR}/tokamak )
set( TOKAMAK_LIBRARIES tokamak )
set( CELLTREE_INCLUDE_DIR ${THIRD_PARTY_DIR}/celltree/src )
set( CELLTREE_LIBRARIES celltree )
set( STANN_INCLUDE_DIR ${THIRD_PARTY_DIR}/stann )
set( GMM_INCLUDE_DIR ${THIRD_PARTY_DIR} )

# mandatory includes
include_directories( ${Teem_INCLUDE_DIRS} )
include_directories( ${XAVIER_SOURCE_DIR} )
include_directories( ${NVIS_INCLUDE_DIR}  )
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${Eigen_INCLUDE_DIR} )
include_directories( ${PNG_INCLUDE_DIRS} )
include_directories( ${Tapkee_INCLUDE_DIR} )
include_directories( ${ExprTK_INCLUDE_DIR} )
include_directories( ${KDTREEPP_INCLUDE_DIR} )
include_directories( ${STANN_INCLUDE_DIR} )
include_directories( ${CELLTREE_INCLUDE_DIR} )
include_directories( ${THIRD_PARTY_DIR} )
add_definitions( ${PNG_DEFINITIONS} )
include_directories( . )

set( AVAILABLE_PACKAGES "" )
set( MISSING_PACKAGES "" )

set( CGARTH_CELLTREE_LIBRARIES celltree )
set( CGARTH_TOKAMAK_LIBRARIES tokamak )


if( VTK_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "VTK" )
	if (VTK_MODULES_ENABLED MATCHES "LIC")
        set( VTK_LIC_FOUND TRUE )
        add_definitions( -D__VTK_HAS_LIC__ )
        list( APPEND ${AVAILABLE_PACKAGES} "VTK_LIC" )
    else()
        set( VTK_LIC_FOUND FALSE )
        list( APPEND ${MISSING_PACKAGES} "VTK_LIC" )
    endif()

    find_path( VTK_DAG
        NAMES vtkDirectedAcyclicGraph.h
        HINTS ${VTK_INCLUDE_DIRS}
    )
    if ( VTK_DAG )
        set( VTK_DAG_FOUND TRUE )
        add_definitions( -D__VTK_HAS_DAG__ )
        list( APPEND ${AVAILABLE_PACKAGES} "VTK_DAG" )
    else()
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
            ${GLUT_INCLUDE_DIR}/GL
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

if( HDF5_FOUND )
    # message( "HDF5 was found: INCLUDE_DIR=" ${HDF5_INCLUDE_DIR} ", LIBRARIES=" ${HDF5_LIBRARIES})
else()
    # message( "HDF5 was not found! INCLUDE_DIRS=" ${HDF5_INCLUDE_DIRS} ", LIBRARIES=" ${HDF5_LIBRARIES})
    if ( HDF5_hdf5_LIBRARY_RELEASE )
        # message( "Found HDF5_hdf5_LIBRARY_RELEASE=" ${HDF5_hdf5_LIBRARY_RELEASE} )
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

#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-result -Wno-deprecated-declarations" )

if( BUILD_OPTION_OPENMP )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
endif()

if( BUILD_OPTION_TBB )
endif()

if ( OLD_COMPILER )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__OLD_C11_COMPILER__")
endif()

#######################################################################
# set build product destinations

set( LIBRARY_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/lib )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin )

link_directories( ${LIBRARY_OUTPUT_PATH} )

# individual libs
list( APPEND XAVIER_LIBRARIES
	boundary-aware-rectgrid
    third_party
    crease
    data
    format
    flow
    geometry
    graphics
    image
    learning
    maps
    #maps_lib
    math
    misc
    # opencl
    #orbital
    #poincare
    reconstruction
    tensor
    topology
    utils
    #VTK

    # The libraries below are only relevant in the context of specific collaborations
    #collab
    #granular
    #waves
)

foreach( lib ${XAVIER_LIBRARIES})
  add_subdirectory( ${lib} )
endforeach()

######################################################################
# determine targets to be build and define their linking
######################################################################

IF( NOT XAVIER_INSTALL_PRESET_IS_SET )
  SET( CMAKE_INSTALL_PREFIX "$ENV{HOME}/XAVIER.default_install_dir" CACHE PATH "Install path prefix, prepended onto install directories." FORCE )
  SET( XAVIER_INSTALL_PRESET_IS_SET 1 CACHE INTERNAL "Default install path written" FORCE )
ENDIF( NOT XAVIER_INSTALL_PRESET_IS_SET )

install( DIRECTORY bin DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY scripts DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY utils/GLSLShaderCode DESTINATION utils USE_SOURCE_PERMISSIONS)


######################################################################
# install all binaries
######################################################################
foreach( target ${XAVIER_TARGETS} )
  install( TARGETS ${target}
        RUNTIME DESTINATION bin
  )
endforeach( target ${XAVIER_TARGETS} )
