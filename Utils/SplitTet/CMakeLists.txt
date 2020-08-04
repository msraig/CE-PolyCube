# Set the project name
project (SplitTet)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    ${CMAKE_SOURCE_DIR}/ScissorPoly/TetStructure.cpp
	${CMAKE_SOURCE_DIR}/ScissorPoly/TetOperation.cpp
    splittet.cpp
)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/Bin)

# Add an executable with the above sources
add_executable(SplitTet ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(SplitTet
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/ScissorPoly
		${CMAKE_SOURCE_DIR}/extern/algebra
)
