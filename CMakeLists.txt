CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(poreflow)

# The version number.
set (POREFLOW_VERSION_MAJOR 1)
set (POREFLOW_VERSION_MINOR 0)

set (CMAKE_MODULE_PATH CMake/Modules ${CMAKE_MODULE_PATH})

enable_language(CXX)

set (POREFLOW_LIBRARIES)

FIND_PACKAGE(VTK REQUIRED)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
  set (POREFLOW_LIBRARIES vtkIO vtkGraphics vtkFiltering vtkCommon vtksys ${POREFLOW_LIBRARIES})
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

file(GLOB CXX_SOURCES src/CTImage.cpp src/writers.cpp)

ADD_EXECUTABLE(./bin/convert_microct src/convert_microct.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(./bin/convert_microct ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(./bin/create_hourglass ./src/create_hourglass.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(./bin/create_hourglass ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(./bin/mesh_microct ./src/mesh_microct.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(./bin/mesh_microct ${POREFLOW_LIBRARIES})

ADD_EXECUTABLE(./bin/tarantula2gmsh ./src/tarantula2gmsh.cpp ${CXX_SOURCES})
TARGET_LINK_LIBRARIES(./bin/tarantula2gmsh ${POREFLOW_LIBRARIES})

