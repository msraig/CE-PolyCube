include (ACGCommon)

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(Boost 1.43.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
  set (USE_IN_SOURCE_TREE_BOOST true)
endif() 

# Add target for first example
add_executable(file_converter EXCLUDE_FROM_ALL ovm_converter.cc)
add_dependencies(file_converter OpenVolumeMesh)
	
if(NOT WIN32)
    
    # Link against all necessary libraries
    target_link_libraries(file_converter OpenVolumeMesh)
	  
    # Set output directory to ${BINARY_DIR}/FileConverter
    set_target_properties(file_converter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/FileConverter)
	
else()

    # Link against all necessary libraries
    target_link_libraries(file_converter OpenVolumeMesh) 

endif()
