find_package( OpenGL REQUIRED )
list( APPEND OPENGL_INCLUDE_DIR ${OPENGL_INCLUDE_DIR}/GL )

find_package( GLUT REQUIRED )
list( APPEND GLUT_INCLUDE_DIR ${GLUT_INCLUDE_DIR}/GL )

find_package( GLEW REQUIRED )

include_directories( ../cuda ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR}
                     ${GLEW_INCLUDE_DIRS})

cuda_add_executable( particles
                     particles.cu
                     trackball.cpp
                     glm.cpp )

target_link_libraries( particles celltree ${OPENGL_gl_LIBRARY} ${GLUT_LIBRARY}
                       ${GLEW_LIBRARIES})
