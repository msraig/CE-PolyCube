find_package(OpenMP REQUIRED)
find_package(Boost 1.60 REQUIRED COMPONENTS filesystem)

include_directories(${PROJECT_SOURCE_DIR} ${MOSEK_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(../extern ../extern/algebra ../extern/TetGen  ../PolyCubeCut ../extern/sus ../extern/MathSuite/SPARSE ../extern/MathSuite/CHOLMOD ../extern/MathSuite/UMFPACK ../extern/MathSuite/HLBFGS ../extern/ANN)
LINK_DIRECTORIES(../Lib ../../PolyCut/ScissorPolyCube/BasicMath/MathSuite/lib64)

FILE(GLOB_RECURSE project_headers *.h ../extern/MathSuite/SPARSE/Sparse_Matrix.h ../extern/MathSuite/SPARSE/Iterative_Solver.h ../extern/MathSuite/SPARSE/Sparse_Solver.h)
FILE(GLOB_RECURSE project_cpps *.cpp ../extern/ANN/*.cpp ../extern/MathSuite/HLBFGS/*.cpp)
set(all_files ${project_headers} ${project_cpps})


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${PROJECT_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${PROJECT_SOURCE_DIR}/Bin)
set(CMAKE_DEBUG_POSTFIX "d")
ADD_EXECUTABLE(ScissorPoly ${all_files})
target_link_libraries(ScissorPoly PRIVATE ${MOSEK_LIBRARIES} gflags::gflags OpenMP::OpenMP_CXX TetGen SUS PolyCubeCut shlwapi optimized OpenMeshCore optimized OpenVolumeMesh optimized SUS_polycube optimized SUS_standard optimized SparseSolver debug OpenMeshCored debug OpenVolumeMeshd debug SUS_polycubed debug SparseSolverd debug SUS_standardd Boost::filesystem)
add_compile_definitions(_USE_MATH_DEFINES)
add_compile_definitions(INCLUDE_TEMPLATES)
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
add_compile_options(/wd4267)
