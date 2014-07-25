CXX = g++

CPPFLAGS = -I/usr/include/vtk-5.8
CXXFLAGS = $(CPPFLAGS) -g -O2 -frounding-math -Wno-deprecated

CGAL_LIBS = -lCGAL -lCGAL_ImageIO -lgmp -lmpfr -lboost_thread -lpthread
VTK_LIBS = -lvtkIO -lvtkHybrid -lvtkGraphics -lvtkFiltering -lvtkCommon

BOOST_LIBS =  -lboost_filesystem -lboost_system
LIBS = $(CGAL_LIBS) $(VTK_LIBS) $(BOOST_LIBS)

objs = writers.o CTImage.o

binaries: ../bin/convert_microct ../bin/mesh_microct ../bin/tarantula2gmsh ../bin/create_hourglass

VPATH = ./:./src

.SUFFIXES: .cpp .F90 .o .lo

.cpp.o:
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

default: $(binaries)

../bin/mesh_microct: $(objs) mesh_microct.cpp
	$(CXX) $(CXXFLAGS) -o ../bin/mesh_microct mesh_microct.cpp $(objs) $(LIBS)

../bin/convert_microct: $(objs) convert_microct.cpp
	$(CXX) $(CXXFLAGS) -o ../bin/convert_microct convert_microct.cpp $(objs) $(LIBS)

../bin/tarantula2gmsh: $(objs) tarantula2gmsh.cpp
	$(CXX) $(CXXFLAGS) -o ../bin/tarantula2gmsh tarantula2gmsh.cpp $(objs) $(LIBS)

../bin/create_hourglass: $(objs) create_hourglass.cpp
	$(CXX) $(CXXFLAGS) -o ../bin/create_hourglass create_hourglass.cpp $(objs) $(LIBS)

clean:
	rm -rf $(binaries) $(objs)
