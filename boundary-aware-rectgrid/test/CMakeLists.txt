cmake_minimum_required(VERSION 3.3)

project(BoundaryAwareRectilinearGrid)

find_package(VTK REQUIRED)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})


find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 
add_definitions(-DSTANDALONE_TEST)

set(EXEC_SRC
	test_BARG
)

file(GLOB SOURCES
    ../*.h
    ../*.cpp
	../*.cxx
	)

foreach( name ${EXEC_SRC} )
	add_executable( ${name} ${name}.cpp ${SOURCES} )
	target_include_directories ( ${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../..)
	target_link_libraries( ${name} ${VTK_LIBRARIES} Eigen3::Eigen )
	#install( TARGETS ${name} RUNTIME DESTINATION bin )
	
	vtk_module_autoinit(
		TARGETS ${name}
		MODULES ${VTK_LIBRARIES}
	)
endforeach()
	

