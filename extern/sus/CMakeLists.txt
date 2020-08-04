include_directories("${PROJECT_SOURCE_DIR}")
FILE(GLOB_RECURSE project_headers *.h)
FILE(GLOB_RECURSE project_cpps *.cpp)
set(all_files ${project_headers} ${project_cpps})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR}/Lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_SOURCE_DIR}/Lib)
set(CMAKE_DEBUG_POSTFIX "d")
ADD_LIBRARY(SUS STATIC ${all_files})
target_link_libraries(SUS OpenMeshCore)
set_target_properties(SUS PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(SUS PROPERTIES "COMPILE_DEFINITIONS" _USE_MATH_DEFINES)
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)




