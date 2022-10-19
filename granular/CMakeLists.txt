set( GRANULAR_CORE_EXEC_SRC
    coord2col
    denis-approx
    denis-oa
    denis-orbits
	denis-orbitalaveraging
    discrete_holmes
    force_chain
    gamma
    gamma2vamp
    layers2col
    liggghts2nrrd
    lowest
    nameparser
    nrrd2txt
    pour
    txt2nrrd
    wave_location
)

# code with third party dependencies
set( GRANULAR_NEEDS_CGARTH_KDTREE_EXEC_SRC
    bump_height
    density2D
    test_locations
)

set( GRANULAR_NEEDS_GL_AND_CGARTH_TOKAMAK_EXEC_SRC
    show_particle_poincare
)

set( GRANULAR_NEEDS_KDTREEPP_EXEC_SRC
    MLSinterpolation
    particleFTLE
    slicer
    volumeFTLE
)

set( GRANULAR_NEEDS_STANN_EXEC_SRC
    contacts
)

set( GRANULAR_NEEDS_VTK_EXEC_SRC
    denis
    granular_visualize
    view_particles_ftle
    discrete_holmes_LIC
)

foreach( target ${GRANULAR_CORE_EXEC_SRC} )
    add_executable( ${target} ${target}.cpp )
    target_link_libraries( ${target} ${Teem_LIBRARIES} )
    install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach()

if( ${KDTREEPP_FOUND} )
    include_directories( ${KDTREEPP_INCLUDE_DIR} )
    foreach( target ${GRANULAR_NEEDS_KDTREEPP_EXEC_SRC} )
        add_executable( ${target} ${target}.cpp )
        target_link_libraries( ${target} ${Teem_LIBRARIES} xavier_math )
        install( TARGETS ${target} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "granular executables: " ${GRANULAR_NEEDS_KDTREEPP_EXEC_SRC} " cannot be built."
        "\nlibkdtree++ not available in this installation."
    )
endif()

if( GL_BUNDLE_FOUND AND CGARTH_TOKAMAK_FOUND )
	include_directories( ${GL_BUNDLE_INCLUDE_DIR} ${CGARTH_TOKAMAK_INCLUDE_DIR} )
    foreach( target ${GRANULAR_NEEDS_GL_AND_CGARTH_TOKAMAK_EXEC_SRC} )
        add_executable( ${target} ${target}.cpp )
        target_link_libraries( ${target} ${Teem_LIBRARIES} ${GL_BUNDLE_LIBRARIES} ${CGARTH_TOKAMAK_LIBRARIES} )
        install( TARGETS ${target} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "granular executables: " ${GRANULAR_NEEDS_GL_AND_CGARTH_TOKAMAK_EXEC_SRC} " cannot be built."
        "\nGL bundle or Tokamak library not available in this installation."
    )
endif()

if( ${STANN_FOUND} )
	include_directories( ${STANN_INCLUDE_DIR} )
    foreach( target ${GRANULAR_NEEDS_STANN_EXEC_SRC} )
        add_executable( ${target} ${target}.cpp )
        target_link_libraries( ${target} ${Teem_LIBRARIES} )
        install( TARGETS ${target} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "granular executables: " ${GRANULAR_NEEDS_STANN_EXEC_SRC} " cannot be built."
        "\nSTANN not available in this installation."
    )	
endif()

if( ${VTK_FOUND} )
	include( ${VTK_USE_FILE} )
    foreach( target ${GRANULAR_NEEDS_VTK_EXEC_SRC} )
        add_executable( ${target} ${target}.cpp )
        target_link_libraries( ${target} ${Teem_LIBRARIES} ${VTK_LIBRARIES}  ${Boost_LIBRARIES} )
        install( TARGETS ${target} RUNTIME DESTINATION bin )
    endforeach()
else()
	message( WARNING "granular executables: " ${GRANULAR_NEEDS_VTK_EXEC_SRC} " cannot be built."
        "\nVTK not available in this installation."
    )		
endif()
