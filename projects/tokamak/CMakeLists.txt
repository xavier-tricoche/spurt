include_directories( ${HDF5_INCLUDE_DIRS} ${NVIS_INCLUDE_DIR} )

add_library( tokamak SHARED
	tokamak_field.cpp
	tokamak_nimrod.cpp
	tokamak_nimrod_parametric.cpp
	poincare_map_rk56.cpp
	tokamak_nrrd.cpp
)

target_link_libraries( tokamak
                       ${Teem_LIBRARIES}
                       ${HDF5_LIBRARIES}
                       ${ZLIB_LIBRARIES} )

install( TARGETS tokamak
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib )
