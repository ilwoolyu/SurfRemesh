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
set(MESH_LIB_DIR /NIRAL/work/ilwoolyu/MeshLib)
include_directories(${MESH_LIB_DIR}/include ${MESH_LIB_DIR}/include/Util)
#add_library(Mesh SHARED IMPORTED)
#set_target_properties(Mesh PROPERTIES IMPORTED_LOCATION ${MESH_LIB_DIR}/res/libMesh.so.1.0)
link_directories(${MESH_LIB_DIR}/res)

include_directories(src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable(SurfRemesh ${Remesh_SOURCES} wrapper/main.cpp ${SRemeshCLP})
target_link_libraries(SurfRemesh Mesh)
