CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(poreflow)

# The version number.
set (POREFLOW_VERSION_MAJOR 1)
set (POREFLOW_VERSION_MINOR 0)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/Modules ${CMAKE_MODULE_PATH})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_language(CXX)

find_program (INTEL_COMPILER icpc)
if (INTEL_COMPILER)
  include (CMakeForceCompiler)
  cmake_force_cxx_compiler (icpc "Intel C++ Compiler")
endif (INTEL_COMPILER)

set (POREFLOW_LIBRARIES)

FIND_PACKAGE(VTK REQUIRED NO_MODULE)
if(VTK_FOUND)
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")

  add_definitions(-DHAVE_VTK)
  include(${VTK_USE_FILE})
  set (POREFLOW_LIBRARIES ${VTK_LIBRARIES} ${POREFLOW_LIBRARIES})
endif()
 
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

FIND_PACKAGE(CGAL REQUIRED)
if(CGAL_FOUND)
  include(${CGAL_USE_FILE})
  set (POREFLOW_LIBRARIES CGAL CGAL_ImageIO mpfr boost_thread pthread ${POREFLOW_LIBRARIES})
endif()

FIND_PACKAGE(Boost REQUIRED)
if(Boost_FOUND)
  include(${Boost_INCLUDE_DIRS})
  set (POREFLOW_LIBRARIES boost_system boost_filesystem ${POREFLOW_LIBRARIES})
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-deprecated -std=c++0x")

include_directories(include)

file(GLOB CXX_SOURCES src/CTImage.cpp src/writers.cpp src/mesh_conversion.cpp)

ADD_EXECUTABLE(convert_microct src/convert_microct.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(convert_microct ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(create_hourglass ./src/create_hourglass.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(create_hourglass ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(mesh_microct ./src/mesh_microct.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(mesh_microct ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(tarantula2gmsh ./src/tarantula2gmsh.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(tarantula2gmsh ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(vtk2gmsh ./src/vtk2gmsh.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(vtk2gmsh ${POREFLOW_LIBRARIES})

