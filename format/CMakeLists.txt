set( FORMAT_NETCDF_SRC
    nc2vtk-triangles
    reorder_float_tagged_filenames
    nek5000toNRRD
    nek5000tonetcdf
    nek5000toVTK
)

set( FORMAT_VTK_NETCDF_SRC
    NCOM2nrrd
    DLR2VTK
)

set( FORMAT_VTK_SRC
    vtk2obj
    vtk2xml
    nrrd2vtk
)

if( NetCDF_FOUND )
    add_library( xavier_nek5000 ${LIBMODE} avtNek5000FileFormat.cpp )
    target_link_libraries( xavier_nek5000 ${NetCDF_LIBRARIES} )
    install( TARGETS xavier_nek5000
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )

	include_directories( ${NetCDF_INCLUDE_DIR} )
    foreach( src ${FORMAT_NETCDF_SRC} )
        add_executable( ${src} ${src}.cpp )
        target_link_libraries( ${src} ${NetCDF_LIBRARIES} ${Teem_LIBRARIES}
		                       ${VTK_NETCDF_LIBRARY} ${VTK_LIBRARIES}
                               xavier_nek5000)
        install( TARGETS ${src} RUNTIME DESTINATION bin )
    endforeach()
else()
    message( WARNING "format executables: " ${FORMAT_SRC} " cannot be built."
        "\nNetCDF not available in this installation."
    )
endif()

if ( VTK_FOUND )
	set( CONVERTER_SRC
		vtk2nrrd
        make_rectilinear_grid
        # nrrd2vtk
		# vtk2vti
		# vti2vtk
	)
	foreach( src ${CONVERTER_SRC} )
		add_executable( ${src} ${src}.cpp )
		target_link_libraries( ${src} ${Teem_LIBRARIES} ${VTK_LIBRARIES}
                               ${Boost_LIBRARIES} )
    	install( TARGETS ${src} RUNTIME DESTINATION bin )
	endforeach()
endif()

if( NetCDF_FOUND AND VTK_FOUND )
    include_directories( ${NetCDF_INCLUDE_DIR} )
    foreach( src ${FORMAT_VTK_NETCDF_SRC} )
        add_executable( ${src} ${src}.cpp )
        target_link_libraries( ${src} ${NetCDF_LIBRARIES} ${Teem_LIBRARIES}
                               ${Boost_LIBRARIES} ${VTK_LIBRARIES} )
        install( TARGETS ${src} RUNTIME DESTINATION bin )
    endforeach()
endif()

if( VTK_FOUND )
    foreach( src ${FORMAT_VTK_SRC} )
        add_executable( ${src} ${src}.cpp )
        target_link_libraries( ${src} ${Teem_LIBRARIES} ${Boost_LIBRARIES}
                               ${VTK_LIBRARIES} )
        install( TARGETS ${src} RUNTIME DESTINATION bin )
    endforeach()
endif()

if ( MatIO_FOUND )
    include_directories( ${MatIO_INCLUDE_DIR} )
    add_executable( matlab2nrrd matlab2nrrd.cpp )
    target_link_libraries( matlab2nrrd ${Boost_LIBRARIES} ${MatIO_LIBRARIES}
                           ${Teem_LIBRARIES} )
    install( TARGETS matlab2nrrd RUNTIME DESTINATION bin )
endif()

add_executable( read_sparse read_sparse.cpp )
target_link_libraries( read_sparse ${Teem_LIBRARIES} ${Boost_LIBRARIES} )
