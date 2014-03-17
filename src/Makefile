CXX = g++

CPPFLAGS = -I/usr/include/vtk-5.8
CXXFLAGS = $(CPPFLAGS) -g -O2 -frounding-math -Wno-deprecated

CGAL_LIBS = -lCGAL -lCGAL_ImageIO -lgmp -lmpfr -lboost_thread -lpthread
VTK_LIBS = -lvtkIO -lvtkHybrid -lvtkGraphics -lvtkFiltering -lvtkCommon

BOOST_LIBS =  -lboost_filesystem -lboost_system
LIBS = $(CGAL_LIBS) $(VTK_LIBS) $(BOOST_LIBS)

binaries: ../bin/mesh_raw ../bin/tarantula2gmsh

default: $(binaries)

../bin/mesh_raw:mesh_raw.cpp
	$(CXX) $(CXXFLAGS) -o ../bin/mesh_raw mesh_raw.cpp $(LIBS)

../bin/tarantula2gmsh: tarantula2gmsh.cpp
	$(CXX) $(CXXFLAGS) -o ../bin/tarantula2gmsh tarantula2gmsh.cpp $(LIBS)
