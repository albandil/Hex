INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(OPENCL_INCLUDE_DIR CL/cl.h PATHS ${OPENCL_INCLUDE_DIRS})
FIND_LIBRARY(OPENCL_LIBRARY NAMES OpenCL ${OPENCL_NAMES} PATHS ${OPENCL_LIBRARIES})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenCL REQUIRED_VARS OPENCL_LIBRARY OPENCL_INCLUDE_DIR)