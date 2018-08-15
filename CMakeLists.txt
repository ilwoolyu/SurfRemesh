cmake_minimum_required(VERSION 3.5)

project(SurfRemesh)

#-----------------------------------------------------------------------------
set(MODULE_NAME SRemesh)

#-----------------------------------------------------------------------------
# MeshLib
find_package(MeshLib REQUIRED)

# SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  )

set(MODULE_SRCS
  SurfaceRemeshing.cpp
  )

set(MODULE_TARGET_LIBRARIES
  Mesh
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
#  add_subdirectory(Testing)
endif()
