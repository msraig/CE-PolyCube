cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)


PROJECT (IO)

find_package(VTK REQUIRED)
if (vtkIOExport_LOADED AND vtkIOImport_LOADED)
vtk_module_config(VTK
  vtkCommonCore
  vtkIOLegacy
  vtkIOXML
  vtkIOGeometry
  #vtkTestingRendering
  vtkInteractionStyle
  vtkIOImport
  vtkIOExport
  vtksys
)
else()
vtk_module_config(VTK
  vtkCommonCore
  vtkIOLegacy
  vtkIOXML
  vtkIOGeometry
  #vtkTestingRendering
  vtkInteractionStyle
  vtkIOImport
  vtksys
)
endif()
include(${VTK_USE_FILE})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/Bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/Bin)

add_executable(vtudecode MACOSX_BUNDLE vtudecode.cxx)
target_link_libraries(vtudecode ${VTK_LIBRARIES})

