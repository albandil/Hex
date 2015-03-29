INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(SQLITE3_INCLUDE_DIR sqlite3.h PATHS ${SQLITE3_INCLUDE_DIRS})
FIND_LIBRARY(SQLITE3_LIBRARY sqlite3 ${SQLITE3_LIBRARIES})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SQLITE3 REQUIRED_VARS SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR)
