# This is the root Spurt CMakeLists file
cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17) # Force CMake to use C++17 to compile

# Run cmake --help-policy CMP<num> to see documentation.
if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
  cmake_policy( SET CMP0012 NEW )
  cmake_policy( SET CMP0167 NEW )
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

string( TOLOWER "${CMAKE_BUILD_TYPE}" build_type_lower )
if( ${build_type_lower} MATCHES "debug" )
    message(WARNING "CMAKE_BUILD_TYPE is ${build_type_lower}")
    add_compile_definitions( SPURT_DEBUG )
else()
    message(WARNING "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
endif()

option( BUILD_OPTION_OPENMP "use openmp parallelization" ON )

option( BUILD_OPTION_TBB "use tbb parallelization" ON )

option( OLD_COMPILER "deactivate \"recent\" C++ features" OFF )

set( LIBMODE SHARED )

# Mac OS does weird stuff with dynamic libs
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MACOSX_RPATH 1)
endif()

#######################################################################

set( Spurt_SOURCE_DIR ${CMAKE_SOURCE_DIR} )

# required packages
set( Boost_USE_STATIC_LIBS      OFF )
set( Boost_USE_MULTITHREADED    ON  )
set( Boost_USE_STATIC_RUNTIME   OFF )
find_package( Boost 1.75 REQUIRED # something fairly recent
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

set( NVIS_INCLUDE_DIR     ${CMAKE_SOURCE_DIR}/third_party/nvis )
set( ExprTK_INCLUDE_DIR   ${CMAKE_SOURCE_DIR}/third_party/fastmathparser )
set( KDTREEPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/kdtree++ )
set( CELLTREE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/celltree/src )
set( CELLTREE_LIBRARIES celltree )
set( TOKAMAK_INCLUDE_DIR ${THIRD_PARTY_DIR}/tokamak )
set( TOKAMAK_LIBRARIES tokamak )
set( GMM_INCLUDE_DIR ${THIRD_PARTY_DIR} )

# mandatory includes
include_directories( 
    .
    ${Spurt_SOURCE_DIR} 
    ${Boost_INCLUDE_DIRS}
    ${Eigen_INCLUDE_DIR}
    ${PNG_INCLUDE_DIRS}
    ${Teem_INCLUDE_DIRS}
    ${NVIS_INCLUDE_DIR} 
    ${KDTREEPP_INCLUDE_DIR}
    ${CELLTREE_INCLUDE_DIR}
    ${THIRD_PARTY_DIR}
)
message( WARNING "include dirs are " ${INCLUDE_DIRECTORIES} ) 
message( WARNING "spurt source dir is " ${Spurt_SOURCE_DIR})
add_definitions( ${PNG_DEFINITIONS} )

set( AVAILABLE_PACKAGES "" )
set( MISSING_PACKAGES "" )

if( VTK_FOUND )
    list( APPEND ${AVAILABLE_PACKAGES} "VTK" )
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

if( BUILD_OPTION_OPENMP )
    if (APPLE)
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp" )
        set( OPENMP_LIBS "-lomp" )
    else()
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
    endif()
endif()

if( BUILD_OPTION_TBB )
    set( TBB_ALL_LIBS "TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy" )
endif()

#######################################################################
# set build product destinations

set( LIBRARY_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/lib )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin )

link_directories( ${LIBRARY_OUTPUT_PATH} )

# individual libs
list( APPEND Spurt_LIBRARIES
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
    vtk
)

foreach( lib ${Spurt_LIBRARIES})
  add_subdirectory( ${lib} )
endforeach()

######################################################################
# determine targets to be build and define their linking
######################################################################

IF( NOT Spurt_INSTALL_PRESET_IS_SET )
  SET( CMAKE_INSTALL_PREFIX "$ENV{HOME}/Spurt.default_install_dir" CACHE PATH "Install path prefix, prepended onto install directories." FORCE )
  SET( Spurt_INSTALL_PRESET_IS_SET 1 CACHE INTERNAL "Default install path written" FORCE )
ENDIF( NOT Spurt_INSTALL_PRESET_IS_SET )

install( DIRECTORY bin DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY scripts DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY utils/GLSLShaderCode DESTINATION utils USE_SOURCE_PERMISSIONS)


######################################################################
# install all binaries
######################################################################
foreach( target ${Spurt_TARGETS} )
  install( TARGETS ${target}
        RUNTIME DESTINATION bin
  )
endforeach( target ${Spurt_TARGETS} )
