# - Try to find C. Garth's Tokamak's code
# Once done this will define
#
#
#   CGARTH_KDTREE_FOUND - system has C. Garth's implementation of kdtree
#   CGARTH_KDTREE_INCLUDE_DIR - C. Garth's kdtree include directory
#
#   CGARTH_TOKAMAK_FOUND - system has C. Garth's Tokamak code
#   CGARTH_TOKAMAK_INCLUDE_DIR - C. Garth's Tokamak code include directory
#   CGARTH_TOKAMAK_LIBRARIES - Link these to use C. Garth's Tokamak code
#
#   CGARTH_CELLTREE_FOUND - system has C. Garth's CellTree code
#   CGARTH_CELLTREE_INCUDE_DIR - C. Garth's CellTree code include directory
#   CGARTH_CELLTREE_LIBRARIES - Link these to use C. Garth's CellTree code

include(Utils)

set( TOKAMAK_LIB_NAME tokamak )
set( CELLTREE_LIB_NAME celltree )

list( APPEND CGARTH_TENTATIVE_SRC_DIR
    $ENV{HOME}/code/christoph
    $ENV{HOME}/code/chris
    $ENV{HOME}/code/cgarth
    $ENV{HOME}/code/garth
)

foreach( dir ${CGARTH_TENTATIVE_SRC_DIR} )
    list( APPEND CGARTH_TENTATIVE_INC_DIR ${dir} )
    list( APPEND CGARTH_TENTATIVE_INC_DIR ${dir}/include )
    list( APPEND CGARTH_TENTATIVE_LIB_DIR ${dir}/lib )
endforeach()

if( CGARTH_USE_THIRD_PARTY )
    set( CGARTH_CELLTREE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/celltree/src )
    set( CGARTH_TOKAMAK_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/tokamak )
    set( CGARTH_KDTREE_INCUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/kdtree )

    set( CGARTH_CELLTREE_LIBRARIES celltree )
    set( CGARTH_TOKAMAK_LIBRARIES tokamak )

    # set( CGARTH_CELLTREE_FIND_QUIETLY TRUE )
    # set( CGARTH_TOKAMAK_FIND_QUIETLY TRUE )
    # set( CGARTH_KDTREE_FIND_QUIETLY TRUE )

    set( CGARTH_KDTREE_FOUND TRUE )
    set( CGARTH_CELLTREE_FOUND TRUE )
    set( CGARTH_TOKAMAK_FOUND TRUE )

    if( CGARTH_TOKAMAK_INCLUDE_DIR AND CGARTH_TOKAMAK_LIBRARIES )
        set( CGARTH_TOKAMAK_FIND_QUIETLY TRUE )
        set( CGARTH_TOKAMAK_FOUND TRUE )
        message( WARNING "Tokamak found - done")
    else()
        # TOKAMAK
        find_path( CGARTH_TOKAMAK_INCLUDE_DIR
            NAMES tokamak/tokamak_nimrod.hpp
            HINTS ${USER_LIKELY_INCLUDE_DIR} ${CGARTH_TENTATIVE_INC_DIR}
        )
        find_library( CGARTH_TOKAMAK_LIBRARIES
            NAME ${TOKAMAK_LIB_NAME}
            HINTS ${USER_LIKELY_LIBRARY_DIR} ${CGARTH_TENTATIVE_LIB_DIR}
        )

        if( CGARTH_TOKAMAK_INCLUDE_DIR AND CGARTH_TOKAMAK_LIBRARIES AND HDF5_FOUND )
            set( CGARTH_TOKAMAK_FOUND TRUE )
            list( APPEND ${CGARTH_TOKAMAK_LIBRARIES} ${HDF5_LIBRARIES} )
        else()
            message( WARNING "Setting tokamak found to false" )
            message( WARNING "HDF5_FOUND: " ${HDF5_FOUND} )
            message( WARNING "HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS} )
            message( WARNING "HDF5_LIBRARIES: " ${HDF5_LIBRARIES})
            set( CGARTH_TOKAMAK_FOUND FALSE )
        endif()

        if( CGARTH_TOKAMAK_FOUND )
            if( TRUE )
                message( STATUS "Found CGARTH_TOKAMAK :"
                    "\nINC: " ${CGARTH_TOKAMAK_INCLUDE_DIR}
                    "\nLIB: " ${CGARTH_TOKAMAK_LIBRARIES}
                )
            endif()
        else()
            if( CGARTH_TOKAMAK_FIND_REQUIRED )
                message( FATAL_ERROR "Could NOT find CGARTH_TOKAMAK" )
            endif()
        endif()
    endif()
	# message( "CGARTH_TOKAMAK is found: " ${CGARTH_TOKAMAK_INCLUDE_DIR} )
endif()

if( CGARTH_KDTREE_INCLUDE_DIR AND CGARTH_KDTREE_LIBRARIES )
    set( CGARTH_KDTREE_FIND_QUIETLY TRUE )
    set( CGARTH_TOKAMAK_FOUND TRUE )
endif()

if( CGARTH_CELLTREE_INCLUDE_DIR AND CGARTH_CELLTREE_LIBRARIES )
    set( CGARTH_CELLTREE_FIND_QUIETLY TRUE )
    set( CGARTH_CELLTREE_FOUND TRUE )
endif()

# KDTREE
find_path( CGARTH_KDTREE_INCLUDE_DIR
    NAMES kdtree/kdtree.hpp
    HINTS ${USER_LIKELY_INCLUDE_DIR} ${CGARTH_TENTATIVE_INC_DIR}
)

if( CGARTH_KDTREE_INCLUDE_DIR )
   set( CGARTH_KDTREE_FOUND TRUE )
else()
   set( CGARTH_KDTREE_FOUND FALSE )
endif()

if( CGARTH_KDTREE_FOUND )
    if( NOT CGARTH_KDTREE_FIND_QUIETLY )
        message( STATUS "Found CGARTH_KDTREE :"
            "\nINC: " ${CGARTH_KDTREE_INCLUDE_DIR}
        )
    endif()
else()
    if( CGARTH_KDTREE_FIND_REQUIRED )
        message( FATAL_ERROR "Could NOT find CGARTH_KDTREE" )
    endif()
endif()

# CELLTREE
find_path( CGARTH_CELLTREE_INCLUDE_DIR
    NAMES celltree.hpp mesh.hpp interpolator.hpp
    HINTS ${USER_LIKELY_INCLUDE_DIR}
          ${CGARTH_TENTATIVE_INC_DIR}
          ${USER_LIKELY_INCLUDE_DIR}/celltree
          ${CGARTH_TENTATIVE_INC_DIR}/celltree
)
find_library( CGARTH_CELLTREE_LIBRARIES
    NAME ${CELLTREE_LIB_NAME}
    HINTS ${USER_LIKELY_LIBRARY_DIR} ${CGARTH_TENTATIVE_LIB_DIR}
)

if( CGARTH_CELLTREE_INCLUDE_DIR AND CGARTH_CELLTREE_LIBRARIES )
   set( CGARTH_CELLTREE_FOUND TRUE )
else()
   set( CGARTH_CELLTREE_FOUND FALSE )
endif()

if( CGARTH_CELLTREE_FOUND )
    if( NOT CGARTH_CELLTREE_FIND_QUIETLY )
        message( STATUS "Found CGARTH_CELLTREE :"
            "\nINC: " ${CGARTH_CELLTREE_INCLUDE_DIR}
            "\nLIB: " ${CGARTH_CELLTREE_LIBRARIES}
        )
    endif()
else()
    if( CGARTH_CELLTREE_FIND_REQUIRED )
        message( FATAL_ERROR "Could NOT find CGARTH_CELLTREE" )
    endif()
endif()

mark_as_advanced( CGARTH_TOKAMAK_INCLUDE_DIR
                  CGARTH_TOKAMAK_LIBRARIES
                  CGARTH_KDTREE_INCLUDE_DIR
                  CGARTH_CELLTREE_INCLUDE_DIR
                  CGARTG_CELLTREE_LIBRARIES )
