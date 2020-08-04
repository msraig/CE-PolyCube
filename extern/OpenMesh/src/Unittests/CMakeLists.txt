include (ACGCommon)

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if ( NOT DEFINED OPENMESH_BUILD_UNIT_TESTS)
    set( OPENMESH_BUILD_UNIT_TESTS false CACHE BOOL "Enable or disable unit test builds in OpenMesh." )
endif()

if ( OPENMESH_BUILD_UNIT_TESTS )
  # Search for gtest headers and libraries
  find_package(GoogleTest)

  if(GTEST_FOUND)

    enable_testing()
   
    # Set correct include paths so that the compiler can find the headers
    include_directories(${GTEST_INCLUDE_DIRS})
    # set additional link directories
    link_directories(${GTEST_LIBRARY_DIR} )

    if ( CMAKE_GENERATOR MATCHES "^Visual Studio 11.*" ) 
      add_definitions( /D _VARIADIC_MAX=10 )
    endif()

    # Create new target named unittests_hexmeshing
    FILE(GLOB UNITTEST_SRC *.cc)
          # Create unittest executable
          acg_add_executable(unittests ${UNITTEST_SRC})
 
          # For the unittest we don't want the install rpath as set by acg_add_executable 
          set_target_properties ( unittests PROPERTIES  BUILD_WITH_INSTALL_RPATH 0 )

   # Set output directory to ${BINARY_DIR}/Unittests
    set (OUTPUT_DIR "${CMAKE_BINARY_DIR}/Unittests")
      set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}) 
      foreach(CONFIG ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${CONFIG} UPCONFIG)
        set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${UPCONFIG} ${OUTPUT_DIR})
    endforeach()

	
	if ( NOT WIN32)
           # Link against all necessary libraries
           target_link_libraries(unittests OpenMeshCore OpenMeshTools ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)
	  
	  
	
       else()
	    # Link against all necessary libraries
		if ( OPENMESH_BUILD_SHARED )
			add_definitions( -DOPENMESHDLL )

		endif()
      target_link_libraries(unittests OpenMeshCore OpenMeshTools ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES}) 
    endif()
	

    
	if ( NOT WIN32 )
	    # Set compiler flags
      set_target_properties(unittests PROPERTIES COMPILE_FLAGS "-g -pedantic -Wno-long-long")
	else()
	    # Set compiler flags
      set_target_properties(unittests PROPERTIES COMPILE_FLAGS "" )

    endif()
    
    if ( OPENMESH_BUILD_SHARED )
	
      #-------- copy dlls to unittests --------
      SET(OPENMESH_TARGETS "OpenMeshTools" "OpenMeshCore")
    
       foreach(TAR ${OPENMESH_TARGETS} )
         add_custom_command(TARGET unittests POST_BUILD 
                            COMMAND "${CMAKE_COMMAND}" -E copy 
                            "$<TARGET_FILE:${TAR}>"
                            "${CMAKE_BINARY_DIR}/Unittests/$<TARGET_FILE_NAME:${TAR}>" 
                            COMMENT "Copying OpenMesh targets to unittests directory")
       endforeach(TAR)
     endif()
	
    acg_copy_after_build(unittests ${CMAKE_CURRENT_SOURCE_DIR}/TestFiles ${CMAKE_BINARY_DIR}/Unittests/)

    add_test(NAME AllTestsIn_OpenMesh_tests WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Unittests" COMMAND "${CMAKE_BINARY_DIR}/Unittests/unittests")

  else(GTEST_FOUND)
      message("Google testing framework was not found!")
  endif(GTEST_FOUND)
endif()
