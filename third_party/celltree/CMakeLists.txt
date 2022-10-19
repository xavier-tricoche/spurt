include_directories( src )

add_subdirectory( src )
#add_subdirectory( benchmark )
#add_subdirectory( tests )

#find_package( CUDA )

if( CUDA_FOUND )
    add_subdirectory( cuda )
#   add_subdirectory( particles )
endif()
