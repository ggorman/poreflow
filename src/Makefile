CXX = g++

CPPFLAGS = -I/usr/include/vtk-5.8
CXXFLAGS = $(CPPFLAGS) -g -O2 -frounding-math -Wno-deprecated

CGAL_LIBS = -lCGAL -lCGAL_ImageIO -lgmp -lmpfr -lboost_thread -lpthread
VTK_LIBS = -lvtkIO -lvtkHybrid -lvtkGraphics -lvtkFiltering -lvtkCommon

LIBS = $(CGAL_LIBS) $(VTK_LIBS)

default:
	$(CXX) $(CXXFLAGS) -o ../bin/tarantula2gmsh tarantula2gmsh.cpp $(LIBS)
	$(CXX) $(CXXFLAGS) -o ../bin/mesh_2d_tiff mesh_2d_tiff.cpp $(LIBS)
	$(CXX) $(CXXFLAGS) -o ../bin/mesh_3d_tiff mesh_3d_tiff.cpp $(LIBS)
