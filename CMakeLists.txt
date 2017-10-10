project(SurfRemesh)
cmake_minimum_required(VERSION 2.8)

find_package(SlicerExecutionModel REQUIRED)

include_directories(/System/Library/Frameworks/Accelerate.framework/Headers)
include_directories(Mesh)
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
include_directories(/System/Library/Frameworks/Accelerate.framework/Headers)

generateCLP(SRemeshCLP wrapper/SRemesh.xml)

set (Remesh_SOURCES src/SurfaceRemeshing.cpp)

# MESH library
find_package(MeshLib REQUIRED)
include_directories(${Mesh_INCLUDE_DIRS} ${Mesh_INCLUDE_DIRS}/Util)
link_directories(${Mesh_LIBRARY_DIRS})

include_directories(src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /home-nfs/masi-shared-home/home/local/VANDERBILT/lyui/bin)
add_executable(SurfRemesh ${Remesh_SOURCES} wrapper/main.cpp ${SRemeshCLP})
target_link_libraries(SurfRemesh Mesh)
