cmake_minimum_required(VERSION 3.10)

project(Flow)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)

include_directories(
	.
	./third_party/nvis
	../../eigen-3.4.0
	/mnt/wslg/distro/usr/include/vtk-9.1
	/mnt/wslg/distro/usr/include/teem
	)

set(LAVD
	flow/lavd.cpp
	image/probe.cpp
	)

find_package(OpenMP REQUIRED)
find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

add_executable(lavd ${LAVD})

target_compile_options(lavd PRIVATE -frounding-math)

link_directories(
	/mnt/wslg/distro/usr/lib/x86_64-linux-gnu
	/mnt/wslg/distro/usr/lib
	)

find_package(Boost REQUIRED COMPONENTS program_options filesystem)


target_link_libraries(lavd
	OpenMP::OpenMP_CXX
	vtkRenderingCore-9.1
	vtkCommonCore-9.1
	vtkCommonDataModel-9.1
	vtkInteractionWidgets-9.1
	vtkRenderingFreeType-9.1
	vtkRenderingAnnotation-9.1
	vtkIOXML-9.1
	vtkIOCore-9.1
	vtkIOLegacy-9.1
	tbb
	teem
	Boost::program_options
	Boost::filesystem
	CGAL
	gmp
	mpfr
	)

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

add_executable(normalize ${NORMALIZER})

target_link_libraries(normalize
	OpenMP::OpenMP_CXX
        vtkRenderingCore-9.1
        vtkCommonCore-9.1
        vtkCommonDataModel-9.1
        vtkInteractionWidgets-9.1
        vtkRenderingFreeType-9.1
        vtkRenderingAnnotation-9.1
        vtkIOXML-9.1
        vtkIOCore-9.1
        vtkIOLegacy-9.1
        tbb
        teem
        Boost::program_options
        Boost::filesystem
        CGAL
        gmp
        mpfr
	)
