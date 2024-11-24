# set( CREASE_SRC
#     measure_wrapper.cpp
#     extractor.cpp
#     pvo.cpp
#     face.cpp
#     crease.cpp
# )
# add_library( spurt_crease ${LIBMODE} ${CREASE_SRC} )
# target_link_libraries( spurt_crease ${Teem_LIBRARIES} spurt_image )
# install( TARGETS spurt_crease
#   LIBRARY DESTINATION lib
#   ARCHIVE DESTINATION lib
# )

if ( VTK_FOUND )
    add_executable( peikert_sadlo ${EXECUTABLE_FLAG} peikert_sadlo.cpp )
    target_link_libraries( peikert_sadlo
                           ${Teem_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${CELLTREE_LIBRARIES}
                           ${NetCDF_LIBRARIES}
                           ${VTK_LIBRARIES}
                           BARG)
    install( TARGETS peikert_sadlo RUNTIME DESTINATION bin )
endif()

add_executable( peikert ${EXECUTABLE_FLAG} peikert.cpp )
target_link_libraries( peikert
                       ${Teem_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${VTK_LIBRARIES} )
install( TARGETS peikert RUNTIME DESTINATION bin )

add_executable( schultz ${EXECUTABLE_FLAG} schultz.cpp )
target_link_libraries( schultz
                       ${Teem_LIBRARIES}
                       ${Boost_LIBRARIES} )
install( TARGETS schultz RUNTIME DESTINATION bin )

add_executable( multires_peikert ${EXECUTABLE_FLAG} multires_peikert.cpp )
target_link_libraries( multires_peikert
                       ${Teem_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${VTK_LIBRARIES}  )
install( TARGETS multires_peikert RUNTIME DESTINATION bin )

# add_executable( crease_surface ${EXECUTABLE_FLAG} crease_surface.cpp )
# target_link_libraries( crease_surface
#                        ${Teem_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${VTK_LIBRARIES} )
# install( TARGETS crease_surface RUNTIME DESTINATION bin )

# add_executable( creaseline ${EXECUTABLE_FLAG} creaseline.cpp )
# target_link_libraries( creaseline
#                        ${Teem_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${VTK_LIBRARIES}
#                        spurt_image
#                        spurt_crease)
# install( TARGETS creaseline RUNTIME DESTINATION bin )
#
# add_executable( creaseline2d ${EXECUTABLE_FLAG} creaseline2d.cpp )
# target_link_libraries( creaseline2d
#                        ${Teem_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${VTK_LIBRARIES}
#                        spurt_image
#                        spurt_crease)
# install( TARGETS creaseline2d RUNTIME DESTINATION bin )
#
# add_executable( creaseline3d ${EXECUTABLE_FLAG} creaseline3d.cpp )
# target_link_libraries( creaseline3d
#                        ${Teem_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${VTK_LIBRARIES}
#                        spurt_image
#                        spurt_crease)
# install( TARGETS creaseline3d RUNTIME DESTINATION bin )

# add_executable( crease_surface_v2 ${EXECUTABLE_FLAG} crease_surface_v2.cpp )
# target_link_libraries( crease_surface_v2
#                        ${Teem_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${VTK_LIBRARIES} )
# install( TARGETS crease_surface_v2 RUNTIME DESTINATION bin )