add_library( celltree SHARED 
             cell.cpp
             dataset.cpp )
             
target_link_libraries( celltree ${NetCDF_LIBRARIES} xavier_nek5000 )

install( TARGETS celltree
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib )
