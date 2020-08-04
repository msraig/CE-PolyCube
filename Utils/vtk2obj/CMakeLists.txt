# Set the project name
project (vtk2obj)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    ${CMAKE_SOURCE_DIR}/ScissorPoly/TetStructure.cpp
    vtk2obj.cpp
)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/Bin)

# Add an executable with the above sources
add_executable(vtk2obj ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(vtk2obj
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/ScissorPoly
		${CMAKE_SOURCE_DIR}/extern/algebra
)
