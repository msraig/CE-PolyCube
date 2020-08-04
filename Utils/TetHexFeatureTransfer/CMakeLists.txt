# Set the project name
project (TetHexFeatureTransfer)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    ${CMAKE_SOURCE_DIR}/ScissorPoly/Hex.cpp
    TetHexFeatureTransfer.cpp
    HexOperation.cpp
)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/Bin)

# Add an executable with the above sources
add_executable(TetHexFeatureTransfer ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(TetHexFeatureTransfer
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/ScissorPoly
		${CMAKE_SOURCE_DIR}/extern/algebra
)

# target_link_libraries(vtk2ovm PRIVATE optimized OpenVolumeMesh debug OpenVolumeMeshd)
add_compile_definitions(_USE_MATH_DEFINES)
add_compile_definitions(INCLUDE_TEMPLATES)
