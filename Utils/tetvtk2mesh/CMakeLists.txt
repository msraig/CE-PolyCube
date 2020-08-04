# Set the project name
project (tetvtk2mesh)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    tetvtk2mesh.cpp
)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/Bin)

# Add an executable with the above sources
add_executable(tetvtk2mesh ${SOURCES})