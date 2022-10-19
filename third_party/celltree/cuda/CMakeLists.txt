include_directories( . )

cuda_add_executable( gputest gputest.cu )
target_link_libraries( gputest celltree )

