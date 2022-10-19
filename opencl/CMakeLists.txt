#######################################################################

find_package( OpenCL )

if ( OpenCL_FOUND )
    find_path( VexCL_INC_DIR "vexcl/backend.hpp" )
    #find_path( MBA_INC_DIR "mba/mba.hpp" )
    # build configuration
    include_directories( ${OpenCL_INCLUDE_DIR} )
	# find_path( SIMD_INC_DIR "simd/vector_types.h" )

    set( OPENCL_SRC
        lorenz_ensemble
        abc_flowmap
        catseye_flowmap
        cr3bp_flowmap
        #multi_bspline_approx
		# unsteady_ftle
    )

    include_directories( ${VexCL_INC_DIR} )
	# include_directories( ${SIMD_INC_DIR}/simd )
    #include_directories( ${MBA_INC_DIR} )

    foreach( name ${OPENCL_SRC} )
        add_executable( ${name} ${name}.cpp )
        target_link_libraries( ${name} ${OpenCL_LIBRARY} ${Boost_LIBRARIES} ${TBB_LIBRARY_RELEASE} ${VTK_LIBRARIES} ${Teem_LIBRARIES} )
    endforeach()
    install( TARGETS ${OPENCL_SRC} RUNTIME DESTINATION bin )
else()
    message( WARNING "OpenCL not found. Disabling corresponding applications" )
endif()
