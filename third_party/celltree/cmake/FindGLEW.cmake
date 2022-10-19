# - try to find GLEW library and include files
#  GLEW_INCLUDE_DIR, where to find GL/GLEW.h, etc.
#  GLEW_LIBRARIES, the libraries to link against
#  GLEW_FOUND, If false, do not try to use GLEW.

IF( WIN32 )
  FIND_PATH( GLEW_INCLUDE_DIR 
             NAMES GL/GLEW.h 
             PATHS ${GLEW_DIR}/include )

  FIND_LIBRARY( GLEW_LIBRARY 
                NAMES GLEW GLEW32
                PATHS ${OPENGL_LIBRARY_DIR}
                      ${GLEW_DIR}/Release )
ELSE( WIN32 )

  IF( APPLE )

    FIND_PATH( GLEW_INCLUDE_DIR GL/glew.h
               /System/Library/Frameworks/GLEW.framework/Versions/A/Headers
               ${OPENGL_INCLUDE_DIR} )
               
    FIND_LIBRARY( GLEW_LIBRARY 
                  NAMES GLEW glew
                  PATHS ${OPENGL_LIBRARY_DIR}
                        ${GLEW_ROOT_PATH}/Release
                        /sw/lib )
              
    SET( GLEW_LIBRARY "-framework GLEW" CACHE STRING "GLEW library for OSX") 

  ELSE( APPLE )
    
    FIND_PATH( GLEW_INCLUDE_DIR GL/glew.h
               /usr/include/GL
               /usr/openwin/share/include
               /usr/openwin/include
               /opt/graphics/OpenGL/include
               /opt/graphics/OpenGL/contrib/libGLEW )
  
    FIND_LIBRARY( GLEW_LIBRARY glew
                  /usr/openwin/lib )
        
  ENDIF( APPLE )
ENDIF( WIN32 )

MESSAGE( STATUS "Found GLEW headers at " ${GLEW_INCLUDE_DIR} )
MESSAGE( STATUS "Found GLEW library at " ${GLEW_LIBRARY} )

SET( GLEW_FOUND "NO" )

IF( GLEW_INCLUDE_DIR AND GLEW_LIBRARY )
    SET( GLEW_FOUND "YES" )
ENDIF( GLEW_INCLUDE_DIR AND GLEW_LIBRARY )
