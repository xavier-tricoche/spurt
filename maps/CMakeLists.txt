add_executable(map4d map4d.cpp)
target_link_libraries(map4d 
    ${Teem_LIBRARIES} 
    ${Boost_LIBRARIES} 
    TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
)
install(TARGETS map4d RUNTIME DESTINATION bin)
