project(SurfRemesh)
cmake_minimum_required(VERSION 2.8)

find_package(SlicerExecutionModel REQUIRED)

include_directories(/System/Library/Frameworks/Accelerate.framework/Headers)
include_directories(Mesh)
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
include_directories(/System/Library/Frameworks/Accelerate.framework/Headers)

generateCLP(SRemeshCLP SRemesh.xml)

set (Remesh_SOURCES SurfaceRemeshing.h SurfaceRemeshing.cpp)

set (Mesh_SOURCES Mesh/Mesh.cpp Mesh/MeshIO.cpp Mesh/MeshIO.h Mesh/ObjIO.cpp Mesh/ObjIO.h Mesh/VtkIO.cpp Mesh/VtkIO.h)

add_executable(SurfRemesh ${Remesh_SOURCES} ${Mesh_SOURCES} main.cpp ${SRemeshCLP})

