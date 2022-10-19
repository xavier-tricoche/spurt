set(MATH_SRC
	taylor.cpp
	matrix.cpp
)
add_library(xavier_math ${LIBMODE} 
	${MATH_SRC} 
)
target_link_libraries(xavier_math ${LAPACK_LIB} teem)

install( TARGETS xavier_math
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)