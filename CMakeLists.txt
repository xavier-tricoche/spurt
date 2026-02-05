cmake_minimum_required(VERSION 3.10)

project(Flow)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" )
set( THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/third_party )
set( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${THIRD_PARTY_DIR} )

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
find_package( CGAL   REQUIRED )
find_package( GMP    REQUIRED )
find_package( MPFR   REQUIRED )
find_package( VTK    REQUIRED )       # graphics / visualization
find_package( ZLIB   REQUIRED )      # I/O
find_package( TBB    REQUIRED )       # parallelism


set( NVIS_INCLUDE_DIR     ${CMAKE_SOURCE_DIR}/third_party/nvis )
set( ExprTK_INCLUDE_DIR   ${CMAKE_SOURCE_DIR}/third_party/fastmathparser )
set( KDTREEPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/kdtree++ )
set( CELLTREE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/celltree/src )
set( CELLTREE_LIBRARIES celltree )

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

# include_directories(
# 	.
# 	./third_party/nvis
# 	../../eigen-3.4.0
# 	/mnt/wslg/distro/usr/include/vtk-9.1
# 	/mnt/wslg/distro/usr/include/teem
# 	)

set(LAVD
	flow/lavd.cpp
	image/probe.cpp
	)

# find_package(OpenMP REQUIRED)
# find_package(CGAL REQUIRED)
# find_package(GMP REQUIRED)
# find_package(MPFR REQUIRED)

add_executable(lavd ${LAVD})
target_link_libraries( lavd ${VTK_LIBRARIES} ${Teem_LIBRARIES}
			${Boost_LIBRARIES} 
            TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
            ${NetCDF_LIBRARIES} )
        install( TARGETS lavd RUNTIME DESTINATION bin )

# target_compile_options(lavd PRIVATE -frounding-math)

# link_directories(
# 	/mnt/wslg/distro/usr/lib/x86_64-linux-gnu
# 	/mnt/wslg/distro/usr/lib
# 	)

# find_package(Boost REQUIRED COMPONENTS program_options filesystem)


# target_link_libraries(lavd
# 	OpenMP::OpenMP_CXX
# 	vtkRenderingCore-9.1
# 	vtkCommonCore-9.1
# 	vtkCommonDataModel-9.1
# 	vtkInteractionWidgets-9.1
# 	vtkRenderingFreeType-9.1
# 	vtkRenderingAnnotation-9.1
# 	vtkIOXML-9.1
# 	vtkIOCore-9.1
# 	vtkIOLegacy-9.1
# 	tbb
# 	teem
# 	Boost::program_options
# 	Boost::filesystem
# 	CGAL
# 	gmp
# 	mpfr
# 	)

#set(VISUALIZER
#	flow/view_flow.cpp
#	)
#
#add_executable(visualizer ${VISUALIZER})
#
#target_link_libraries(visualizer
#	OpenMP::OpenMP_CXX
#	vtkRenderingCore-9.1
#	vtkCommonCore-9.1
#	vtkCommonDataModel-9.1
#	vtkInteractionWidgets-9.1
#	vtkRenderingFreeType-9.1
#	vtkRenderingAnnotation-9.1
#	vtkIOXML-9.1
#	vtkIOCore-9.1
#	vtkIOLegacy-9.1
#	vtkFiltersCore-9.1
#	vtkCommonExecutionModel-9.1
#	vtkImagingCore-9.1
#	vtkIOImage-9.1
#	tbb
#	teem
#	Boost::program_options
#	Boost::filesystem
#	)

set(NORMALIZER
	flow/field_normalizer.cpp
	image/probe.cpp
	)

message( WARNING "Boost_LIBRARIES are ${BOOST_LIBRARIES}")

add_executable(normalize ${NORMALIZER})
target_link_libraries( normalize ${VTK_LIBRARIES} ${Teem_LIBRARIES}
			Boost::filesystem Boost::program_options Boost::system 
            TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
            ${NetCDF_LIBRARIES} )
        install( TARGETS normalize RUNTIME DESTINATION bin )

# target_link_libraries(normalize
# 	OpenMP::OpenMP_CXX
#         vtkRenderingCore-9.1
#         vtkCommonCore-9.1
#         vtkCommonDataModel-9.1
#         vtkInteractionWidgets-9.1
#         vtkRenderingFreeType-9.1
#         vtkRenderingAnnotation-9.1
#         vtkIOXML-9.1
#         vtkIOCore-9.1
#         vtkIOLegacy-9.1
#         tbb
#         teem
#         Boost::program_options
#         Boost::filesystem
#         CGAL
#         gmp
#         mpfr
# 	)
