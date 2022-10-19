set(GRAIN_SRC
    m3c_v1.1
    grain-stat
)
foreach( src ${GRAIN_SRC} )
    add_executable( ${src} ${src}.cpp )
    target_link_libraries( ${src} ${Teem_LIBRARIES} )
    install( TARGETS ${src} RUNTIME DESTINATION bin )
endforeach()

set(GEOMETRY_SRC
    boundary_extraction
    add_shear_stress
    refine_mesh
)
foreach( src ${GEOMETRY_SRC} )
    add_executable( ${src} ${src}.cpp )
    target_link_libraries( ${src} ${Teem_LIBRARIES} ${Boost_LIBRARIES}
                           ${VTK_LIBRARIES} ${VTK_NETCDF_LIBRARY}
                           ${NetCDF_LIBRARIES})
    install( TARGETS ${src} RUNTIME DESTINATION bin )
endforeach()
