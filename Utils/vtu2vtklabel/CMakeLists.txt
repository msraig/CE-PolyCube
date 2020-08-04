# Set the project name
project (vtu2vtklabel)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    ${CMAKE_SOURCE_DIR}/ScissorPoly/TetStructure.cpp
    vtu2vtklabel.cpp
)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/Bin)

# Add an executable with the above sources
add_executable(vtu2vtklabel ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(vtu2vtklabel
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/ScissorPoly
		${CMAKE_SOURCE_DIR}/extern/algebra
		${CMAKE_SOURCE_DIR}/extern/rapidxml
)

target_link_libraries(vtu2vtklabel)
add_compile_definitions(_USE_MATH_DEFINES)
add_compile_definitions(INCLUDE_TEMPLATES)
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
