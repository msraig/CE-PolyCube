include_directories("${PROJECT_SOURCE_DIR}")
FILE(GLOB_RECURSE project_headers *.h)
FILE(GLOB_RECURSE project_cpps *.cpp)
set(all_files ${project_headers} ${project_cpps})


SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR}/Lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_SOURCE_DIR}/Lib)
set(CMAKE_DEBUG_POSTFIX "d")
ADD_LIBRARY(PolyCubeCut STATIC ${all_files})

set_target_properties(PolyCubeCut PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(PolyCubeCut ${CGAL_LIBRARIES} TetGen)

