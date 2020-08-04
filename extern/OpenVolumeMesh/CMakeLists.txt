cmake_minimum_required (VERSION 2.6)

if("${PROJECT_NAME}" STREQUAL "")
    message("Setting project name to OpenVolumeMesh")
    project (OpenVolumeMesh)
endif()

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# include our cmake files
include (ACGCommon)

acg_get_version ()

include (ACGOutput)

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# ========================================================================
# Call the subdirectories with the projects
# ========================================================================

add_subdirectory (src)

# Only add examples and doc to build target if OpenVolumeMesh
# is not built as external library only
if(${PROJECT_NAME} MATCHES "OpenVolumeMesh")
    add_subdirectory (examples)
    add_subdirectory (documentation)
endif()

# ========================================================================
# Bundle generation (Targets exist, now configure them)
# ========================================================================

# display results
acg_print_configure_header (OPENVOLUMEMESH "OpenVolumeMesh")
