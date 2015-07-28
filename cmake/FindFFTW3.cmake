INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(FFTW3_INCLUDE_DIR fftw3.h ${FFTW3_INCLUDE_DIRS})
FIND_LIBRARY(FFTW3_LIBRARY NAMES fftw3 ${FFTW3_NAMES} PATHS ${FFTW3_LIBRARIES})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3 REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR)