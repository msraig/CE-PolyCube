# Add a library to the project using the specified source files. 
# In Linux/Unix, it will creates the libtet.a

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR}/Lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_SOURCE_DIR}/Lib)
set(CMAKE_DEBUG_POSTFIX "d")
add_library(TetGen STATIC tetgen.cxx predicates.cxx tetgen.h)
set_target_properties(TetGen PROPERTIES POSITION_INDEPENDENT_CODE ON)