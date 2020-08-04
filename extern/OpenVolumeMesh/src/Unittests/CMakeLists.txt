include (ACGCommon)

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Search for gtest headers and libraries
find_package(GoogleTest)

if(GTEST_FOUND)
   
    # Set correct include paths so that the compiler can find the headers
    include_directories(${GTEST_INCLUDE_DIRS})
    # set additional link directories
    link_directories(${GTEST_LIBRARY_DIR} )
    # Create new target named unittests_hexmeshing
    add_executable(unittests EXCLUDE_FROM_ALL unittests.cc unittests_common.cc)
	
	if ( NOT WIN32)
      # Link against all necessary libraries
      target_link_libraries(unittests OpenVolumeMesh ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)
	  
	  # Set output directory to ${BINARY_DIR}/Unittests
      set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Unittests)
	
    else()
	  # Link against all necessary libraries
      target_link_libraries(unittests OpenVolumeMesh ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES}) 
	endif()
	

    
	if ( NOT WIN32 )
	  # Set compiler flags
      set_target_properties(unittests PROPERTIES COMPILE_FLAGS "-Wall -Wextra -Wuninitialized -Wunused -g -pedantic -ansi -Wno-long-long")
	else()
	  # Set compiler flags
      set_target_properties(unittests PROPERTIES COMPILE_FLAGS "" )

	  # copy exe file to "Build" directory
      # Visual studio will create this file in a subdirectory so we can't use
      # RUNTIME_OUTPUT_DIRECTORY directly here
      add_custom_command (TARGET unittests POST_BUILD
                          COMMAND ${CMAKE_COMMAND} -E
                                  copy_if_different
                                    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/unittests.exe
                                    ${CMAKE_BINARY_DIR}/Unittests/unittests.exe)
	endif()
	
    acg_copy_after_build(unittests ${CMAKE_CURRENT_SOURCE_DIR}/TestFiles ${CMAKE_BINARY_DIR}/Unittests/)

else(GTEST_FOUND)
    message("Google testing framework was not found!")
endif(GTEST_FOUND)
