# Set project variable such that OpenVolumeMesh does not
# build its unit tests and stuff
project(Examples)

# Add OpenVolumeMesh's source path to include path for compiler
include_directories(../src)

# Add target for first example
add_executable(simple_mesh EXCLUDE_FROM_ALL simple_mesh/simple_mesh.cc)

if(NOT WIN32)
  # Link against OpenVolumeMesh
  target_link_libraries(simple_mesh OpenVolumeMesh)
  
  # Set output directory to ${BINARY_DIR}/Examples
  set_target_properties(simple_mesh PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Examples)

else()
  # Link against OpenVolumeMesh
  target_link_libraries(simple_mesh OpenVolumeMesh) 
endif()



if(WIN32)
  # copy exe file to "Build" directory
  # Visual studio will create this file in a subdirectory so we can't use
  # RUNTIME_OUTPUT_DIRECTORY directly here
  add_custom_command (TARGET simple_mesh POST_BUILD
                      COMMAND ${CMAKE_COMMAND} -E
                              copy_if_different
                                ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/simple_mesh.exe
                                ${CMAKE_BINARY_DIR}/Examples/simple_mesh.exe)
endif()
