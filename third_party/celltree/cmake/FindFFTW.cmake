# - Try to find fftw3
# Once done this will define
#
#  FFTW_FOUND - system has fftw
#  FFTW_INCLUDE_DIR - the BZip2 include directory
#  FFTW_LIBRARIES - Link these to use BZip2
#  FFTW_FLAGS

IF( FFTW_INCLUDE_DIR AND FFTW_LIBRARIES )
    SET( FFTW_FIND_QUIETLY TRUE )
ENDIF( FFTW_INCLUDE_DIR AND FFTW_LIBRARIES )

FIND_PATH( FFTW_INCLUDE_DIR fftw3.h )
FIND_LIBRARY( FFTW_LIBRARIES NAMES fftw3 )

IF( FFTW_INCLUDE_DIR AND FFTW_LIBRARIES )
   SET( FFTW_FOUND TRUE )
   #INCLUDE( CheckLibraryExists )
   #CHECK_LIBRARY_EXISTS( ${FFTW_LIBRARIES} fftw_kernel_malloc "" FFTW_FOUND )
ELSE( FFTW_INCLUDE_DIR AND FFTW_LIBRARIES )
   SET( FFTW_FOUND FALSE )
ENDIF( FFTW_INCLUDE_DIR AND FFTW_LIBRARIES )

IF( FFTW_FOUND )
  IF( NOT FFTW_FIND_QUIETLY )
    MESSAGE( STATUS "Found FFTW : ${FFTW_LIBRARIES}" )
  ENDIF( NOT FFTW_FIND_QUIETLY )
ELSE( FFTW_FOUND )
  IF( FFTW_FIND_REQUIRED )
    MESSAGE( FATAL_ERROR "Could NOT find FFTW" )
  ENDIF( FFTW_FIND_REQUIRED)
ENDIF( FFTW_FOUND )

MARK_AS_ADVANCED( FFTW_INCLUDE_DIR FFTW_LIBRARIES )
