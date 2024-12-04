add_executable(symplectic4d symplectic4d.cpp)
target_link_libraries(symplectic4d 
    ${Teem_LIBRARIES} 
    ${Boost_LIBRARIES} 
    TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
)
install(TARGETS symplectic4d RUNTIME DESTINATION bin)


add_executable(map4d map4d.cpp)
target_link_libraries(map4d 
    ${Teem_LIBRARIES} 
    ${Boost_LIBRARIES} 
    TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
)
install(TARGETS map4d RUNTIME DESTINATION bin)
