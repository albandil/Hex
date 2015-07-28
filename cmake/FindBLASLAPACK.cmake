INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE(OpenBLAS)
IF(OPENBLAS_FOUND)
    SET(BLASLAPACK_INCLUDE_DIR ${OPENBLAS_INCLUDE_DIR})
    SET(BLASLAPACK_LIBRARIES ${OPENBLAS_LIBRARY_PATH})
ELSE(OPENBLAS_FOUND)
    FIND_PACKAGE(BLAS)
    FIND_PACKAGE(LAPACK)
    SET(BLASLAPACK_INCLUDE_DIR "")
    SET(BLASLAPACK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
ENDIF(OPENBLAS_FOUND)