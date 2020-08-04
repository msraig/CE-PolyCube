set (CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)

macro(source_group_by_dir source_files)
 if(MSVC)
  set(sgbd_cur_dir ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(sgbd_file ${${source_files}})
   string(REGEX REPLACE ${sgbd_cur_dir}/\(.*\) \\1 sgbd_fpath ${sgbd_file})
   string(REGEX REPLACE "\(.*\)/.*" \\1 sgbd_group_name ${sgbd_fpath})
   string(COMPARE EQUAL ${sgbd_fpath} ${sgbd_group_name} sgbd_nogroup)
   string(REPLACE "/" "\\" sgbd_group_name ${sgbd_group_name})
   if(sgbd_nogroup)
    set(sgbd_group_name "\\")
   endif(sgbd_nogroup)
   source_group(${sgbd_group_name} FILES ${sgbd_file})
  endforeach(sgbd_file)
 endif(MSVC)
endmacro(source_group_by_dir)


if(NOT CMAKE_BUILD_TYPE) 
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

AUX_SOURCE_DIRECTORY(MeshStructure SOURCES)
AUX_SOURCE_DIRECTORY(../extern/ANN SOURCES)
SET(SOURCES ${SOURCES})

FILE(GLOB_RECURSE project_headers MeshStructure/*.h ../extern/ANN/*.h)
FILE(GLOB_RECURSE project_cpps MeshStructure/*.cpp ../extern/ANN/*.cpp)
set(all_files ${project_headers} ${project_cpps})
source_group_by_dir(all_files)

INCLUDE_DIRECTORIES(MeshStructure ../extern/ANN/)

if(${MSVC}) 
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR}/Lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_SOURCE_DIR}/Lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_SOURCE_DIR}/Bin)
elseif(${CYGWIN})
SET(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Bin)
SET(LIBRARY_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Lib)
elseif(${MINGW})
SET(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Bin)
SET(LIBRARY_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Lib)
elseif(${UNIX})
SET(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Bin)
SET(LIBRARY_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Lib)
endif(${MSVC})


set(CMAKE_DEBUG_POSTFIX "d")

ADD_LIBRARY(HexOpt STATIC ${all_files})
set_target_properties(HexOpt PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(HexOpt PRIVATE OpenMP::OpenMP_CXX)

ADD_EXECUTABLE(HexaConvert Driver/HexaConvert.cpp)
target_link_libraries(HexaConvert PRIVATE OpenMP::OpenMP_CXX HexOpt)

ADD_EXECUTABLE(HexaRefiner Driver/HexaRefiner.cpp)
target_link_libraries(HexaRefiner PRIVATE OpenMP::OpenMP_CXX HexOpt)
