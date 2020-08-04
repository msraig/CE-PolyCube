cmake_minimum_required(VERSION 3.12)

project(ScissorPolyCube)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Please install mosek 8, CGAL, Boost
# Please set BOOST_ROOT CGAL_DIR in your system path if they are not set automatically.


find_package(CGAL REQUIRED)
find_package(MOSEK REQUIRED)
include(${CGAL_USE_FILE})

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${BOOST_LIBRARYDIR})
endif()
if (CGAL_FOUND)
    include_directories(${CGAL_DIR})
endif()

include_directories(${CMAKE_SOURCE_DIR}/extern)
include_directories(${CMAKE_SOURCE_DIR}/extern/OpenMesh/src)
include_directories(${CMAKE_SOURCE_DIR}/extern/OpenVolumeMesh/src)
include_directories(${CMAKE_SOURCE_DIR}/extern/eigen3)

add_subdirectory(${CMAKE_SOURCE_DIR}/extern/gflags)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/OpenMesh)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/OpenVolumeMesh)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/tetgen)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/eigen3)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/SUS)
add_subdirectory(${CMAKE_SOURCE_DIR}/PolyCubeCut)
add_subdirectory(${CMAKE_SOURCE_DIR}/ScissorPoly)
add_subdirectory(${CMAKE_SOURCE_DIR}/HexOptimizer)
add_subdirectory(${CMAKE_SOURCE_DIR}/Utils)
